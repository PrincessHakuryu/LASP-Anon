#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nnp_vasp_like_ifc2_ifc3.py

VASP-like finite-displacement workflow using an ASE/NNP calculator.

Given a POSCAR and an ASE calculator factory, this script generates:
  - FORCE_CONSTANTS_2ND   via phonopy VASP-style displacements + fake vasprun.xml forces
  - NNP-relaxed POSCAR    optional ASE optimization before IFC generation
  - FORCE_CONSTANTS_3RD   via thirdorder_vasp.py sow/reap + fake vasprun.xml forces
  - NNP_BTE_REPORT.yaml/txt with FC2 imaginary-frequency/soft-mode diagnostics

This intentionally keeps the traditional VASP/thirdorder/ShengBTE workflow:
  phonopy / thirdorder generate displaced POSCARs;
  your ASE calculator replaces VASP single-point force calculations;
  fake vasprun.xml files contain the NNP forces in VASP-compatible units (eV/Angstrom);
  phonopy / thirdorder_vasp.py parse those XML files and write IFCs.

Requirements:
  pip/conda packages: ase, numpy
  executables in PATH: phonopy, thirdorder_vasp.py

Example:
  1. Edit the USER CALCULATOR CONFIG section near the top of this file.
  2. Run, for example:

  python nnp_vasp_like_batch_ifc_shengbte_singlefile.py \
      --batch-root . \
      --fc2-dim 2 2 1 \
      --fc3-dim 2 2 1 \
      --fc3-cutoff -3 \
      --qgrid 12 12 12 \
      --overwrite

The script has an internal make_calc() function, so no external calc_factory.py
is required. If desired, --calc-factory path/to/file.py:make_calc can still be
used as an override.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import contextlib
import multiprocessing as mp
import importlib.util
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np
import yaml
from ase.io import read, write
import datetime

# =============================================================================
# USER CALCULATOR CONFIG
# =============================================================================
# Edit this section when moving the script or changing models.
#
# Default assumes the compiled NequIP model is placed in the SAME directory as
# this script. You can also override the model path from command line:
#   --model-path /path/to/NequIP-OAM-L-ASE.nequip.pt2
#
# Supported CALCULATOR_KIND values in make_calc():
#   "nequip_compiled"
#   "nequix"
#   "custom"  -> edit make_calc() manually.

CALCULATOR_KIND = "nequix"
NNP_MODEL_PATH = "NequIP-OAM-L-ASE.nequip.pt2"
NNP_DEVICE = "cuda"

# For Nequix, edit these if CALCULATOR_KIND = "nequix".
NEQUIX_MODEL_NAME = "nequix-oam-1-pft"
NEQUIX_BACKEND = "torch"
NEQUIX_USE_KERNEL = True


def make_calc():
    """
    Return an ASE-compatible calculator.

    This is the only function you usually need to edit. It replaces the old
    external calc_factory.py:make_calc file.
    """
    kind = str(CALCULATOR_KIND).lower().strip()

    if kind == "nequip":
        from nequip.integrations.ase import NequIPCalculator
        import torch._inductor.codecache
        return NequIPCalculator.from_compiled_model(
            str(NNP_MODEL_PATH),
            device=str(NNP_DEVICE),
        )

    if kind == "nequix":
        from nequix.calculator import NequixCalculator

        return NequixCalculator(
            NEQUIX_MODEL_NAME,
            model_path="nequix-oam-1-pft.nqx",
            backend=NEQUIX_BACKEND,
            use_kernel=NEQUIX_USE_KERNEL,
        )

    if kind == "custom":
        # Put your custom ASE calculator here, for example:
        #
        # from my_package import MyCalculator
        # return MyCalculator(...)
        #
        raise NotImplementedError("Edit make_calc() for CALCULATOR_KIND='custom'.")

    raise ValueError(f"Unknown CALCULATOR_KIND: {CALCULATOR_KIND!r}")


# =============================================================================
# Utilities
# =============================================================================

def natural_int_from_name(path: Path, default: int = 0) -> int:
    """Return the last integer appearing in a filename."""
    nums = re.findall(r"\d+", path.name)
    return int(nums[-1]) if nums else default


def natural_sort(paths: Iterable[Path]) -> list[Path]:
    return sorted(paths, key=lambda p: (natural_int_from_name(p), str(p)))


def run_checked(
    cmd: Sequence[str],
    cwd: Path,
    input_text: str | None = None,
    log_path: Path | None = None,
) -> subprocess.CompletedProcess:
    """Run a command and raise with useful log tail if it fails."""
    cwd = Path(cwd)
    if log_path is not None:
        log_path = Path(log_path)
        with log_path.open("w", encoding="utf-8") as log:
            proc = subprocess.run(
                list(cmd),
                cwd=cwd,
                input=input_text,
                text=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
    else:
        proc = subprocess.run(
            list(cmd),
            cwd=cwd,
            input=input_text,
            text=True,
            check=False,
        )

    if proc.returncode != 0:
        tail = ""
        if log_path is not None and log_path.is_file():
            lines = log_path.read_text(encoding="utf-8", errors="ignore").splitlines()
            tail = "\n".join(lines[-120:])
        raise RuntimeError(
            "Command failed:\n"
            f"  cwd: {cwd}\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  returncode: {proc.returncode}\n"
            + (f"\nLast log lines from {log_path}:\n{tail}\n" if tail else "")
        )

    return proc


def load_calc_factory(spec: str) -> Callable[[], object]:
    """
    Load a calculator factory from:
      path/to/calc_factory.py:make_calc
    or:
      importable.module:make_calc
    """
    if ":" not in spec:
        raise ValueError(
            "--calc-factory must look like 'calc_factory.py:make_calc' "
            "or 'some.module:make_calc'"
        )

    module_spec, func_name = spec.split(":", 1)

    if module_spec.endswith(".py") or Path(module_spec).is_file():
        module_path = Path(module_spec).resolve()
        if not module_path.is_file():
            raise FileNotFoundError(module_path)
        name = module_path.stem
        spec_obj = importlib.util.spec_from_file_location(name, module_path)
        if spec_obj is None or spec_obj.loader is None:
            raise ImportError(f"Could not import {module_path}")
        module = importlib.util.module_from_spec(spec_obj)
        sys.modules[name] = module
        spec_obj.loader.exec_module(module)
    else:
        module = __import__(module_spec, fromlist=[func_name])

    factory = getattr(module, func_name)
    if not callable(factory):
        raise TypeError(f"{spec} did not resolve to a callable")
    return factory


def write_minimal_vasprun_forces(filename: Path, forces: np.ndarray) -> None:
    """
    Write minimal vasprun.xml containing only:
      <modeling><calculation><varray name="forces">...</varray></calculation></modeling>

    This is sufficient for thirdorder_vasp.py read_forces and commonly works
    for phonopy force extraction from VASP XML.
    """
    forces = np.asarray(forces, dtype=float)
    if forces.ndim != 2 or forces.shape[1] != 3:
        raise ValueError(f"forces must have shape (natoms, 3), got {forces.shape}")

    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)

    with filename.open("w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write("<modeling>\n")
        fh.write("  <calculation>\n")
        fh.write('    <varray name="forces">\n')
        for fx, fy, fz in forces:
            fh.write(f"      <v> {fx:.16e} {fy:.16e} {fz:.16e} </v>\n")
        fh.write("    </varray>\n")
        fh.write("  </calculation>\n")
        fh.write("</modeling>\n")


def evaluate_poscar_force(
    poscar: Path,
    calc,
    xml_out: Path,
    copy_poscar_to: Path | None = None,
    forces_out: Path | None = None,
) -> np.ndarray:
    """
    Read POSCAR, evaluate forces using ASE calculator, and write vasprun.xml.

    Important: do not sort, relax, symmetrize, or wrap atom order here.
    thirdorder/phonopy need the exact displaced-supercell atom order.
    """
    atoms = read(str(poscar), format="vasp")
    atoms.pbc = True
    atoms.calc = calc
    forces = np.asarray(atoms.get_forces(), dtype=float)

    if copy_poscar_to is not None:
        copy_poscar_to.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(poscar, copy_poscar_to)

    if forces_out is not None:
        forces_out.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(forces_out, forces, fmt="%.16e")

    write_minimal_vasprun_forces(xml_out, forces)
    return forces


def ensure_clean_workdir(workdir: Path, overwrite: bool) -> None:
    workdir = Path(workdir)
    if workdir.exists():
        if not overwrite:
            raise FileExistsError(
                f"Workdir already exists: {workdir}\n"
                "Use --overwrite to remove it first, or choose a new --workdir."
            )
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)



def _extract_displacements_from_phonopy_yaml(path: Path, natoms_supercell: int) -> list[tuple[int, np.ndarray]]:
    """
    Read displacement atom index and displacement vector from phonopy_disp.yaml.

    Supports common phonopy YAML layouts:
      phonopy_dataset.first_atoms[].number        # usually 0-based
      displacement_dataset.first_atoms[].number   # usually 0-based
      displacements[].atom or displacements[].number

    Returns:
      list of (atom_number_1based, displacement_vector)
    """
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    raw: list[tuple[int, np.ndarray]] = []

    dataset = (
        data.get("phonopy_dataset")
        or data.get("displacement_dataset")
        or data.get("dataset")
        or {}
    )

    first_atoms = dataset.get("first_atoms") or data.get("first_atoms")
    if first_atoms is not None:
        for entry in first_atoms:
            if "number" not in entry:
                raise RuntimeError(f"Missing 'number' in phonopy displacement entry: {entry}")
            if "displacement" not in entry:
                raise RuntimeError(f"Missing 'displacement' in phonopy displacement entry: {entry}")
            raw.append((int(entry["number"]), np.asarray(entry["displacement"], dtype=float)))

    elif "displacements" in data:
        for entry in data["displacements"]:
            if "atom" in entry:
                atom_num = int(entry["atom"])
            elif "number" in entry:
                atom_num = int(entry["number"])
            else:
                raise RuntimeError(f"Missing atom/number in phonopy displacement entry: {entry}")
            if "displacement" not in entry:
                raise RuntimeError(f"Missing 'displacement' in phonopy displacement entry: {entry}")
            raw.append((atom_num, np.asarray(entry["displacement"], dtype=float)))

    else:
        raise RuntimeError(
            f"Could not find displacements in {path}. "
            "Please inspect phonopy_disp.yaml format."
        )

    if not raw:
        raise RuntimeError(f"No displacements found in {path}")

    nums = [x[0] for x in raw]
    min_num, max_num = min(nums), max(nums)

    # Normalize to FORCE_SETS 1-based atom numbering.
    # phonopy_dataset.first_atoms[].number is normally 0-based.
    if min_num == 0:
        normalized = [(num + 1, disp) for num, disp in raw]
    elif min_num >= 1 and max_num <= natoms_supercell:
        normalized = raw
    elif max_num <= natoms_supercell - 1:
        normalized = [(num + 1, disp) for num, disp in raw]
    else:
        raise RuntimeError(
            f"Cannot infer atom-index base from {path}: min={min_num}, max={max_num}, "
            f"natoms_supercell={natoms_supercell}"
        )

    for atom_num, disp in normalized:
        if atom_num < 1 or atom_num > natoms_supercell:
            raise RuntimeError(
                f"Invalid atom number {atom_num} in displacement yaml "
                f"for natoms_supercell={natoms_supercell}"
            )
        if disp.shape != (3,):
            raise RuntimeError(f"Invalid displacement vector shape: {disp.shape}")

    return normalized


def _write_force_sets_type1(
    filename: Path,
    displacements: list[tuple[int, np.ndarray]],
    force_blocks: list[np.ndarray],
    subtract_drift: bool = True,
) -> None:
    """
    Write phonopy FORCE_SETS type-1 format directly.

    This avoids relying on phonopy's VASP vasprun.xml parser for fake XML files.
    """
    filename = Path(filename)
    if not force_blocks:
        raise RuntimeError("No force blocks to write to FORCE_SETS")

    natoms = int(np.asarray(force_blocks[0]).shape[0])
    if len(displacements) != len(force_blocks):
        raise RuntimeError(
            f"Number of displacements ({len(displacements)}) does not match "
            f"number of force blocks ({len(force_blocks)})"
        )

    with filename.open("w", encoding="utf-8") as fh:
        fh.write(f"{natoms}\n")
        fh.write(f"{len(force_blocks)}\n")

        for (atom_num, disp), forces in zip(displacements, force_blocks):
            forces = np.asarray(forces, dtype=float)
            if forces.shape != (natoms, 3):
                raise RuntimeError(f"Bad force block shape {forces.shape}; expected ({natoms}, 3)")

            if subtract_drift:
                forces = forces - forces.mean(axis=0, keepdims=True)

            fh.write("\n")
            fh.write(f"{int(atom_num)}\n")
            fh.write(f"  {disp[0]: .16e} {disp[1]: .16e} {disp[2]: .16e}\n")
            for fx, fy, fz in forces:
                fh.write(f"  {fx: .16e} {fy: .16e} {fz: .16e}\n")



# =============================================================================
# ASE/NNP relaxation and harmonic-stability diagnostics
# =============================================================================

def _atoms_basic_summary(atoms) -> dict:
    """Small JSON/YAML-friendly summary for an ASE Atoms object."""
    atoms.pbc = True
    cell = np.asarray(atoms.cell.array, dtype=float)
    lengths = atoms.cell.lengths()
    angles = atoms.cell.angles()
    return {
        "formula": atoms.get_chemical_formula(mode="hill"),
        "natoms": int(len(atoms)),
        "volume_A3": float(atoms.get_volume()),
        "cell_lengths_A": [float(x) for x in lengths],
        "cell_angles_deg": [float(x) for x in angles],
        "cell_matrix_A": cell.tolist(),
    }


def _safe_energy_forces_stress(atoms) -> dict:
    """Evaluate calculator energy/forces/stress when available."""
    info: dict = {}
    try:
        e = float(atoms.get_potential_energy())
        info["energy_eV"] = e
        info["energy_per_atom_eV"] = e / max(len(atoms), 1)
    except Exception as exc:
        info["energy_error"] = repr(exc)

    try:
        forces = np.asarray(atoms.get_forces(), dtype=float)
        fnorm = np.linalg.norm(forces, axis=1)
        info["max_force_eV_per_A"] = float(fnorm.max()) if len(fnorm) else 0.0
        info["mean_force_eV_per_A"] = float(fnorm.mean()) if len(fnorm) else 0.0
        info["force_drift_eV_per_A"] = forces.mean(axis=0).tolist() if len(forces) else [0.0, 0.0, 0.0]
    except Exception as exc:
        info["forces_error"] = repr(exc)

    try:
        # ASE stress order is xx yy zz yz xz xy, normally in eV/A^3.
        stress = np.asarray(atoms.get_stress(voigt=True), dtype=float)
        info["stress_voigt_eV_per_A3"] = stress.tolist()
        # 1 eV/A^3 = 1602.1766208 kbar.
        info["stress_voigt_kbar"] = (stress * 1602.1766208).tolist()
        info["max_abs_stress_kbar"] = float(np.max(np.abs(stress * 1602.1766208)))
    except Exception as exc:
        info["stress_error"] = repr(exc)

    return info


def _get_optimizer_class(name: str):
    """Return ASE optimizer class by simple name."""
    lname = str(name).lower().strip()
    if lname == "fire":
        from ase.optimize import FIRE
        return FIRE
    if lname == "bfgs":
        from ase.optimize import BFGS
        return BFGS
    if lname == "lbfgs":
        from ase.optimize import LBFGS
        return LBFGS
    raise ValueError(f"Unknown relaxation optimizer: {name!r}. Use fire, bfgs, or lbfgs.")


def relax_structure_with_ase(
    *,
    input_poscar: Path,
    output_poscar: Path,
    calc,
    fmax: float,
    steps: int,
    optimizer_name: str,
    relax_cell: bool,
    constant_volume: bool,
    verbose: bool,
) -> dict:
    """
    Relax a structure on the NNP potential energy surface using ASE.

    If relax_cell=True but the calculator cannot provide stress, this function falls
    back to atom-only relaxation and records the fallback in the returned report.
    """
    input_poscar = Path(input_poscar)
    output_poscar = Path(output_poscar)
    output_poscar.parent.mkdir(parents=True, exist_ok=True)

    atoms = read(str(input_poscar), format="vasp")
    atoms.pbc = True
    atoms.calc = calc

    report: dict = {
        "enabled": True,
        "input_poscar": str(input_poscar),
        "output_poscar": str(output_poscar),
        "requested_relax_cell": bool(relax_cell),
        "constant_volume": bool(constant_volume),
        "optimizer": str(optimizer_name),
        "fmax_target_eV_per_A": float(fmax),
        "max_steps": int(steps),
        "initial_structure": _atoms_basic_summary(atoms),
        "initial_calculator_eval": _safe_energy_forces_stress(atoms),
    }

    dyn_target = atoms
    used_cell_filter = False
    fallback_reason = None

    if relax_cell:
        # Check stress once before constructing a cell filter. Some NNP calculators
        # provide forces but not stress.
        try:
            _ = atoms.get_stress(voigt=True)
            try:
                from ase.filters import FrechetCellFilter
                dyn_target = FrechetCellFilter(atoms, constant_volume=bool(constant_volume))
                used_cell_filter = True
            except Exception:
                try:
                    from ase.filters import ExpCellFilter
                    dyn_target = ExpCellFilter(atoms, constant_volume=bool(constant_volume))
                    used_cell_filter = True
                except Exception:
                    from ase.filters import UnitCellFilter
                    dyn_target = UnitCellFilter(atoms, constant_volume=bool(constant_volume))
                    used_cell_filter = True
        except Exception as exc:
            fallback_reason = (
                "Calculator did not provide stress; falling back to atom-only relaxation: "
                + repr(exc)
            )
            dyn_target = atoms
            used_cell_filter = False

    report["used_cell_filter"] = bool(used_cell_filter)
    if fallback_reason:
        report["fallback_reason"] = fallback_reason

    log_path = output_poscar.parent / "relax.log"
    traj_path = output_poscar.parent / "relax.traj"
    opt_cls = _get_optimizer_class(optimizer_name)

    if verbose:
        mode = "cell+ions" if used_cell_filter else "ions only"
        print(f"\n=== NNP relaxation: {mode}, optimizer={optimizer_name}, fmax={fmax:g}, steps={steps} ===")

    opt = opt_cls(dyn_target, logfile=str(log_path), trajectory=str(traj_path))
    converged = bool(opt.run(fmax=float(fmax), steps=int(steps)))

    atoms.wrap()
    write(str(output_poscar), atoms, format="vasp", vasp5=True, direct=True, sort=False)

    report.update(
        {
            "converged": bool(converged),
            "optimizer_steps_taken": int(getattr(opt, "nsteps", -1)),
            "relax_log": str(log_path),
            "relax_trajectory": str(traj_path),
            "final_structure": _atoms_basic_summary(atoms),
            "final_calculator_eval": _safe_energy_forces_stress(atoms),
        }
    )

    if verbose:
        ff = report["final_calculator_eval"].get("max_force_eV_per_A", None)
        ep = report["final_calculator_eval"].get("energy_per_atom_eV", None)
        print(f"NNP relaxation converged: {converged}")
        if ff is not None:
            print(f"Final max force: {ff:.6g} eV/A")
        if ep is not None:
            print(f"Final energy/atom: {ep:.10g} eV")

    return report


def no_relax_structure_copy(*, input_poscar: Path, output_poscar: Path, calc) -> dict:
    """Copy input POSCAR to output POSCAR and evaluate the NNP at the copied structure."""
    input_poscar = Path(input_poscar)
    output_poscar = Path(output_poscar)
    output_poscar.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_poscar, output_poscar)
    atoms = read(str(output_poscar), format="vasp")
    atoms.pbc = True
    atoms.calc = calc
    return {
        "enabled": False,
        "input_poscar": str(input_poscar),
        "output_poscar": str(output_poscar),
        "structure": _atoms_basic_summary(atoms),
        "calculator_eval": _safe_energy_forces_stress(atoms),
    }


def _write_simple_conf(path: Path, text: str) -> None:
    path = Path(path)
    path.write_text(text.strip() + "\n", encoding="utf-8")


def _phonopy_gamma_modes_from_qpoints(qpoints_yaml: Path, poscar: Path, top_n: int = 8) -> list[dict]:
    """Read Gamma frequencies/eigenvectors and classify translation-like vs optical-like modes."""
    with Path(qpoints_yaml).open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    atoms = read(str(poscar), format="vasp")
    symbols = atoms.get_chemical_symbols()
    masses = np.asarray(atoms.get_masses(), dtype=float)

    phonon0 = data["phonon"][0]
    bands = phonon0["band"]
    nat_eig = len(bands[0].get("eigenvector", []))
    if nat_eig and nat_eig != len(masses):
        # Fall back to phonopy YAML masses when present.
        pts = data.get("points") or data.get("atoms") or []
        if len(pts) == nat_eig:
            symbols = [p.get("symbol", p.get("name", f"atom{i+1}")) for i, p in enumerate(pts)]
            masses = np.asarray([float(p.get("mass", p.get("atomic_mass", 1.0))) for p in pts], dtype=float)
        else:
            # Last-resort neutral masses to keep reporting alive.
            symbols = [f"atom{i+1}" for i in range(nat_eig)]
            masses = np.ones(nat_eig, dtype=float)

    modes: list[dict] = []
    for ib, band in enumerate(bands, start=1):
        freq = float(band["frequency"])
        mode: dict = {"band": int(ib), "frequency_THz": freq}

        eig = band.get("eigenvector")
        if eig is not None:
            ev = np.array(
                [[complex(c[0], c[1]) for c in atom] for atom in eig],
                dtype=complex,
            )
            if ev.shape[0] == len(masses):
                u = ev / np.sqrt(masses)[:, None]
                avg = u.mean(axis=0)
                denom = max(float(np.linalg.norm(u)), 1e-30)
                ratio = float(np.linalg.norm(u - avg) / denom)
                contrib = np.linalg.norm(u, axis=1) ** 2
                if float(contrib.sum()) > 0:
                    contrib = contrib / contrib.sum()
                top = np.argsort(contrib)[-int(top_n):][::-1]
                mode["translation_ratio"] = ratio
                mode["top_atom_contributions"] = [
                    {
                        "atom_index_1based": int(i + 1),
                        "symbol": str(symbols[i]) if i < len(symbols) else f"atom{i+1}",
                        "contribution_fraction": float(contrib[i]),
                    }
                    for i in top
                ]

        modes.append(mode)

    return modes


def _parse_mesh_yaml(mesh_yaml: Path, imag_tol: float) -> dict:
    """Parse phonopy mesh.yaml and count imaginary modes."""
    with Path(mesh_yaml).open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    phonons = data["phonon"]
    nq_irred = len(phonons)
    nband = len(phonons[0]["band"]) if phonons else 0
    total_weight = int(sum(int(p.get("weight", 1)) for p in phonons))
    total_weighted_modes = int(total_weight * nband)

    neg_raw: list[tuple[float, int, int, int, list[float]]] = []
    neg_sig: list[tuple[float, int, int, int, list[float]]] = []
    by_band_raw: dict[int, int] = {}
    by_band_sig: dict[int, int] = {}
    weighted_by_band_sig: dict[int, int] = {}

    min_freq = None
    min_entry = None

    for iq, p in enumerate(phonons):
        q = [float(x) for x in p.get("q-position", [0.0, 0.0, 0.0])]
        w = int(p.get("weight", 1))
        for ib, band in enumerate(p["band"], start=1):
            f = float(band["frequency"])
            if min_freq is None or f < min_freq:
                min_freq = f
                min_entry = {"frequency_THz": f, "band": int(ib), "iq": int(iq), "weight": int(w), "q_position": q}
            if f < -1.0e-5:
                neg_raw.append((f, ib, iq, w, q))
                by_band_raw[ib] = by_band_raw.get(ib, 0) + 1
            if f < -abs(float(imag_tol)):
                neg_sig.append((f, ib, iq, w, q))
                by_band_sig[ib] = by_band_sig.get(ib, 0) + 1
                weighted_by_band_sig[ib] = weighted_by_band_sig.get(ib, 0) + w

    neg_raw_sorted = sorted(neg_raw, key=lambda x: x[0])
    neg_sig_sorted = sorted(neg_sig, key=lambda x: x[0])

    weighted_neg_raw = int(sum(x[3] for x in neg_raw))
    weighted_neg_sig = int(sum(x[3] for x in neg_sig))
    max_branch_weighted_sig = max(weighted_by_band_sig.values()) if weighted_by_band_sig else 0

    return {
        "nq_irreducible": int(nq_irred),
        "nbands": int(nband),
        "total_q_weight": int(total_weight),
        "total_weighted_modes": int(total_weighted_modes),
        "min_frequency_THz": float(min_freq) if min_freq is not None else None,
        "min_frequency_entry": min_entry,
        "imag_tol_THz": float(abs(imag_tol)),
        "negative_entries_raw_lt_minus_1e_5": int(len(neg_raw)),
        "negative_entries_significant": int(len(neg_sig)),
        "weighted_negative_entries_raw": int(weighted_neg_raw),
        "weighted_negative_entries_significant": int(weighted_neg_sig),
        "weighted_significant_fraction_of_all_modes": (
            float(weighted_neg_sig / total_weighted_modes) if total_weighted_modes else 0.0
        ),
        "max_single_branch_weighted_significant_fraction_of_qmesh": (
            float(max_branch_weighted_sig / total_weight) if total_weight else 0.0
        ),
        "negative_count_by_band_raw": {int(k): int(v) for k, v in sorted(by_band_raw.items())},
        "negative_count_by_band_significant": {int(k): int(v) for k, v in sorted(by_band_sig.items())},
        "weighted_negative_count_by_band_significant": {
            int(k): int(v) for k, v in sorted(weighted_by_band_sig.items())
        },
        "most_negative_20_raw": [
            {
                "frequency_THz": float(f),
                "band": int(ib),
                "iq": int(iq),
                "weight": int(w),
                "q_position": [float(x) for x in q],
            }
            for f, ib, iq, w, q in neg_raw_sorted[:20]
        ],
        "most_negative_20_significant": [
            {
                "frequency_THz": float(f),
                "band": int(ib),
                "iq": int(iq),
                "weight": int(w),
                "q_position": [float(x) for x in q],
            }
            for f, ib, iq, w, q in neg_sig_sorted[:20]
        ],
    }


def classify_harmonic_stability(
    *,
    gamma_modes: list[dict],
    mesh: dict,
    imag_tol: float,
    soft_mode_threshold: float,
    optical_ratio_threshold: float,
    soft_branch_qfraction: float,
    many_imag_fraction: float,
) -> dict:
    """
    Classify NNP harmonic stability. This is a screening label only.

    The workflow intentionally continues to FC3/CONTROL regardless of this label.
    """
    imag_tol = abs(float(imag_tol))
    soft_thr = abs(float(soft_mode_threshold))
    optical_thr = float(optical_ratio_threshold)

    gamma_soft_optical = []
    gamma_tiny_acoustic = []
    for m in gamma_modes:
        f = float(m.get("frequency_THz", 0.0))
        ratio = m.get("translation_ratio")
        if f < -soft_thr and ratio is not None and float(ratio) >= optical_thr:
            gamma_soft_optical.append(m)
        if f < -1e-5 and f >= -imag_tol and ratio is not None and float(ratio) < optical_thr:
            gamma_tiny_acoustic.append(m)

    min_freq = mesh.get("min_frequency_THz")
    neg_sig = int(mesh.get("negative_entries_significant", 0))
    neg_raw = int(mesh.get("negative_entries_raw_lt_minus_1e_5", 0))
    frac_sig = float(mesh.get("weighted_significant_fraction_of_all_modes", 0.0))
    branch_frac = float(mesh.get("max_single_branch_weighted_significant_fraction_of_qmesh", 0.0))

    if gamma_soft_optical:
        label = "SOFT_GAMMA_OPTICAL"
        severity = "red"
        reason = (
            "At least one Gamma mode is strongly imaginary and optical-like "
            f"(frequency < -{soft_thr:g} THz and translation_ratio >= {optical_thr:g})."
        )
    elif (min_freq is not None and float(min_freq) < -soft_thr) or branch_frac >= float(soft_branch_qfraction) or frac_sig >= float(many_imag_fraction):
        label = "SOFT_BRANCH_OR_MANY_IMAG"
        severity = "red"
        reason = (
            "Significant imaginary modes form a soft branch or occupy a non-negligible "
            "part of the q mesh."
        )
    elif neg_sig > 0:
        label = "BORDERLINE_SIGNIFICANT_IMAG"
        severity = "yellow"
        reason = f"There are {neg_sig} mesh entries below -{imag_tol:g} THz, but no strong Gamma optical flag."
    elif neg_raw > 0:
        # Collapse former "green-yellow" into green. Tiny imaginary modes
        # above the significant threshold are treated as numerical/acoustic drift.
        label = "STABLE_NNP_HARMONIC"
        severity = "green"
        reason = (
            "Only tiny imaginary frequencies were detected below -1e-5 THz but above "
            f"the significant threshold -{imag_tol:g} THz; classified as green."
        )
    else:
        label = "STABLE_NNP_HARMONIC"
        severity = "green"
        reason = f"No imaginary mesh entries below -{imag_tol:g} THz."

    return {
        "label": label,
        "severity": severity,
        "reason": reason,
        "will_continue_to_fc3_and_bte_inputs": True,
        "thresholds": {
            "imag_tol_THz": imag_tol,
            "soft_mode_threshold_THz": soft_thr,
            "optical_translation_ratio_threshold": optical_thr,
            "soft_branch_qfraction_threshold": float(soft_branch_qfraction),
            "many_imag_fraction_threshold": float(many_imag_fraction),
        },
        "gamma_soft_optical_modes": gamma_soft_optical,
        "gamma_tiny_acoustic_like_negative_modes": gamma_tiny_acoustic,
    }


def analyze_fc2_harmonic_stability(
    *,
    workdir: Path,
    phonopy_exe: str,
    fc2_dim: Sequence[int],
    mesh_dims: Sequence[int],
    imag_tol: float,
    soft_mode_threshold: float,
    optical_ratio_threshold: float,
    soft_branch_qfraction: float,
    many_imag_fraction: float,
    verbose: bool,
) -> dict:
    """Run phonopy Gamma+mesh checks and return a report dictionary."""
    workdir = Path(workdir)
    dim_str = " ".join(str(int(x)) for x in fc2_dim)
    mesh_str = " ".join(str(int(x)) for x in mesh_dims)

    report: dict = {
        "analysis_status": "not_started",
        "fc2_dim": [int(x) for x in fc2_dim],
        "mesh": [int(x) for x in mesh_dims],
        "frequency_unit": "THz",
    }

    if verbose:
        print("\n=== FC2 harmonic stability / soft-mode screen ===")
        print(f"Stability mesh: {mesh_str}; imag_tol={imag_tol:g} THz")

    try:
        gamma_conf = workdir / "q_gamma_check.conf"
        _write_simple_conf(
            gamma_conf,
            f"""
DIM = {dim_str}
QPOINTS = 0 0 0
EIGENVECTORS = .TRUE.
""",
        )
        run_checked(
            [phonopy_exe, "-c", "POSCAR", "--readfc", str(gamma_conf.name)],
            cwd=workdir,
            log_path=workdir / "phonopy_gamma_check.log",
        )
        qpoints_yaml = workdir / "qpoints.yaml"
        qpoints_copy = workdir / "qpoints_gamma_check.yaml"
        if qpoints_yaml.is_file():
            shutil.copy2(qpoints_yaml, qpoints_copy)
        else:
            raise RuntimeError("phonopy Gamma check did not create qpoints.yaml")

        gamma_modes = _phonopy_gamma_modes_from_qpoints(
            qpoints_yaml=qpoints_copy,
            poscar=workdir / "POSCAR",
            top_n=8,
        )

        mesh_conf = workdir / "mesh_stability_check.conf"
        _write_simple_conf(
            mesh_conf,
            f"""
DIM = {dim_str}
MESH = {mesh_str}
GAMMA_CENTER = .TRUE.
WRITE_MESH = .TRUE.
MESH_FORMAT = YAML
""",
        )
        run_checked(
            [phonopy_exe, "-c", "POSCAR", "--readfc", str(mesh_conf.name)],
            cwd=workdir,
            log_path=workdir / "phonopy_mesh_stability_check.log",
        )
        mesh_yaml = workdir / "mesh.yaml"
        mesh_copy = workdir / "mesh_stability_check.yaml"
        if mesh_yaml.is_file():
            shutil.copy2(mesh_yaml, mesh_copy)
        else:
            raise RuntimeError("phonopy mesh check did not create mesh.yaml")

        mesh_report = _parse_mesh_yaml(mesh_copy, imag_tol=imag_tol)
        classification = classify_harmonic_stability(
            gamma_modes=gamma_modes,
            mesh=mesh_report,
            imag_tol=imag_tol,
            soft_mode_threshold=soft_mode_threshold,
            optical_ratio_threshold=optical_ratio_threshold,
            soft_branch_qfraction=soft_branch_qfraction,
            many_imag_fraction=many_imag_fraction,
        )

        report.update(
            {
                "analysis_status": "ok",
                "gamma_modes": gamma_modes,
                "mesh_summary": mesh_report,
                "classification": classification,
            }
        )

        if verbose:
            label = classification["label"]
            minf = mesh_report.get("min_frequency_THz")
            neg = mesh_report.get("negative_entries_significant")
            print(f"Harmonic-stability label: {label}")
            print(f"Minimum mesh frequency: {minf:.6g} THz; significant negative entries: {neg}")

    except Exception as exc:
        report.update(
            {
                "analysis_status": "failed",
                "error": repr(exc),
                "classification": {
                    "label": "ANALYSIS_FAILED_CONTINUING",
                    "severity": "unknown",
                    "reason": "FC2 was generated, but the phonopy stability-analysis step failed.",
                    "will_continue_to_fc3_and_bte_inputs": True,
                },
            }
        )
        if verbose:
            print(f"WARNING: FC2 stability analysis failed but workflow will continue: {exc}")

    return report


def write_case_reports(workdir: Path, report: dict) -> tuple[Path, Path]:
    """Write YAML and text reports in each *-BTE workdir."""
    workdir = Path(workdir)
    yaml_path = workdir / "NNP_BTE_REPORT.yaml"
    txt_path = workdir / "NNP_BTE_REPORT.txt"

    with yaml_path.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(report, fh, sort_keys=False, allow_unicode=True)

    lines: list[str] = []
    lines.append("NNP/Anon BTE generation report")
    lines.append("=" * 34)
    lines.append(f"case_name: {report.get('case_name', '')}")
    lines.append(f"source_poscar: {report.get('source_poscar', '')}")
    lines.append(f"created_at: {report.get('created_at', '')}")
    lines.append(f"workflow_status: {report.get('workflow_status', '')}")
    lines.append("")

    relax = report.get("relaxation", {}) or {}
    lines.append("[NNP relaxation]")
    lines.append(f"enabled: {relax.get('enabled')}")
    if relax.get("enabled"):
        lines.append(f"requested cell relaxation: {relax.get('requested_relax_cell')}")
        lines.append(f"used cell filter: {relax.get('used_cell_filter')}")
        lines.append(f"converged: {relax.get('converged')}")
        if relax.get("fallback_reason"):
            lines.append(f"fallback_reason: {relax.get('fallback_reason')}")
        fini = relax.get("final_calculator_eval", {}) or {}
        if "energy_per_atom_eV" in fini:
            lines.append(f"final_energy_per_atom_eV: {fini['energy_per_atom_eV']:.10g}")
        if "max_force_eV_per_A" in fini:
            lines.append(f"final_max_force_eV_per_A: {fini['max_force_eV_per_A']:.6g}")
        if "max_abs_stress_kbar" in fini:
            lines.append(f"final_max_abs_stress_kbar: {fini['max_abs_stress_kbar']:.6g}")
    else:
        ev = relax.get("calculator_eval", {}) or {}
        if "max_force_eV_per_A" in ev:
            lines.append(f"input_max_force_eV_per_A: {ev['max_force_eV_per_A']:.6g}")
    lines.append("")

    stab = report.get("harmonic_stability", {}) or {}
    cls = stab.get("classification", {}) or {}
    mesh = stab.get("mesh_summary", {}) or {}
    lines.append("[FC2 harmonic stability]")
    lines.append(f"analysis_status: {stab.get('analysis_status')}")
    lines.append(f"label: {cls.get('label')}")
    lines.append(f"severity: {cls.get('severity')}")
    lines.append(f"reason: {cls.get('reason')}")
    if mesh:
        lines.append(f"min_frequency_THz: {mesh.get('min_frequency_THz')}")
        lines.append(f"negative_entries_raw_lt_minus_1e_5: {mesh.get('negative_entries_raw_lt_minus_1e_5')}")
        lines.append(f"negative_entries_significant: {mesh.get('negative_entries_significant')}")
        lines.append(f"weighted_negative_entries_significant: {mesh.get('weighted_negative_entries_significant')}")
        lines.append(f"weighted_significant_fraction_of_all_modes: {mesh.get('weighted_significant_fraction_of_all_modes')}")
        lines.append(f"negative_count_by_band_raw: {mesh.get('negative_count_by_band_raw')}")
        lines.append("most_negative_5_raw:")
        for row in (mesh.get("most_negative_20_raw") or [])[:5]:
            lines.append(
                f"  f={row['frequency_THz']:.6g} THz band={row['band']} "
                f"q={row['q_position']} weight={row['weight']}"
            )
    lines.append("")

    lines.append("[Generated final files]")
    for key in ("FORCE_CONSTANTS_2ND", "FORCE_CONSTANTS_3RD", "CONTROL"):
        p = report.get("outputs", {}).get(key)
        lines.append(f"{key}: {p}")
    lines.append("")
    lines.append("Note: FC3/CONTROL are generated regardless of the stability label; use the label only for later filtering.")

    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return yaml_path, txt_path

# =============================================================================
# IFC2: phonopy VASP-like route
# =============================================================================

def generate_ifc2_with_phonopy_vasp_like(
    *,
    workdir: Path,
    calc,
    phonopy_exe: str,
    fc2_dim: Sequence[int],
    fc2_distance: float,
    verbose: bool,
) -> Path:
    """
    Generate FORCE_CONSTANTS_2ND using:
      phonopy -d --dim="a b c" --amplitude distance -c POSCAR
      NNP forces on POSCAR-*
      phonopy -f fc2/disp-*/vasprun.xml
      phonopy --writefc --dim="a b c" -c POSCAR
      cp FORCE_CONSTANTS FORCE_CONSTANTS_2ND
    """
    workdir = Path(workdir)
    dim_str = " ".join(str(int(x)) for x in fc2_dim)

    if verbose:
        print("\n=== IFC2: phonopy VASP-like displacement workflow ===")
        print(f"FC2 supercell dim : {dim_str}")
        print(f"FC2 displacement  : {fc2_distance:g} Å")

    # Generate displacements.
    cmd = [
        phonopy_exe,
        "-d",
        f"--dim={dim_str}",
        f"--amplitude={float(fc2_distance)}",
        "-c",
        "POSCAR",
    ]
    run_checked(cmd, cwd=workdir, log_path=workdir / "phonopy_fc2_displace.log")

    poscars = [
        p for p in workdir.glob("POSCAR-*")
        if re.fullmatch(r"POSCAR-\d+", p.name)
    ]
    poscars = natural_sort(poscars)
    if not poscars:
        raise RuntimeError("phonopy did not generate POSCAR-### displacement files.")

    if verbose:
        print(f"FC2 displaced POSCARs: {len(poscars)}")

    # Evaluate FC2 forces. We still write fake vasprun.xml files for bookkeeping,
    # but we DO NOT ask phonopy to parse them. Some modern phonopy versions require
    # much more complete vasprun.xml files. Instead, write FORCE_SETS directly.
    force_blocks: list[np.ndarray] = []
    for count, poscar in enumerate(poscars, start=1):
        idx = natural_int_from_name(poscar)
        disp_dir = workdir / "fc2" / f"disp-{idx:05d}"
        xml_path = disp_dir / "vasprun.xml"
        forces = evaluate_poscar_force(
            poscar=poscar,
            calc=calc,
            xml_out=xml_path,
            copy_poscar_to=disp_dir / "POSCAR",
            forces_out=disp_dir / "forces.dat",
        )
        force_blocks.append(np.asarray(forces, dtype=float))
        if verbose:
            print(f"[FC2] {count:5d}/{len(poscars):5d} {poscar.name} -> {disp_dir / 'forces.dat'}")

    # Create FORCE_SETS directly from phonopy_disp.yaml + NNP force blocks.
    natoms_supercell = int(force_blocks[0].shape[0])
    displacements = _extract_displacements_from_phonopy_yaml(
        workdir / "phonopy_disp.yaml",
        natoms_supercell=natoms_supercell,
    )
    if len(displacements) != len(force_blocks):
        raise RuntimeError(
            f"phonopy_disp.yaml has {len(displacements)} displacements, "
            f"but evaluated {len(force_blocks)} force blocks."
        )

    _write_force_sets_type1(
        workdir / "FORCE_SETS",
        displacements=displacements,
        force_blocks=force_blocks,
        subtract_drift=True,
    )

    if verbose:
        print(f"FC2 FORCE_SETS written directly: {workdir / 'FORCE_SETS'}")

    # Create FORCE_CONSTANTS from FORCE_SETS.
    cmd = [
        phonopy_exe,
        f"--dim={dim_str}",
        "--writefc",
        "--full-fc",
        "-c",
        "POSCAR",
    ]
    run_checked(cmd, cwd=workdir, log_path=workdir / "phonopy_fc2_writefc.log")

    fc2_src = workdir / "FORCE_CONSTANTS"
    if not fc2_src.is_file():
        raise RuntimeError("phonopy --writefc finished but FORCE_CONSTANTS was not created.")

    fc2_dst = workdir / "FORCE_CONSTANTS_2ND"
    shutil.copy2(fc2_src, fc2_dst)

    _validate_fc2_text_header_for_shengbte(
        fc2_file=fc2_dst,
        poscar=workdir / "POSCAR",
        fc2_dim=fc2_dim,
    )

    if verbose:
        print(f"IFC2 written: {fc2_dst}")
        print("IFC2 full-supercell header check: OK")

    return fc2_dst


# =============================================================================
# IFC3: thirdorder_vasp.py VASP-like route
# =============================================================================

def generate_ifc3_with_thirdorder_vasp_like(
    *,
    workdir: Path,
    calc,
    thirdorder_exe: str,
    fc3_dim: Sequence[int],
    fc3_cutoff: str,
    verbose: bool,
) -> Path:
    """
    Generate FORCE_CONSTANTS_3RD using:
      thirdorder_vasp.py sow na nb nc cutoff
      NNP forces on 3RD.POSCAR.*
      find/list XMLs | thirdorder_vasp.py reap na nb nc cutoff
    """
    workdir = Path(workdir)
    dims = [str(int(x)) for x in fc3_dim]

    if verbose:
        print("\n=== IFC3: thirdorder_vasp.py VASP-like workflow ===")
        print(f"FC3 supercell dim : {' '.join(dims)}")
        print(f"FC3 cutoff        : {fc3_cutoff}  (e.g. -3 = third nearest neighbor)")

    # Generate thirdorder displacements.
    cmd = [thirdorder_exe, "sow"] + dims + [str(fc3_cutoff)]
    run_checked(cmd, cwd=workdir, log_path=workdir / "thirdorder_sow.log")

    # Typical names are 3RD.POSCAR.00001, 3RD.POSCAR.00002, ...
    poscars = [
        p for p in workdir.glob("3RD.POSCAR.*")
        if p.is_file() and not p.name.endswith(".xml")
    ]
    poscars = natural_sort(poscars)
    if not poscars:
        raise RuntimeError(
            "thirdorder_vasp.py sow did not generate any 3RD.POSCAR.* files."
        )

    if verbose:
        print(f"FC3 displaced POSCARs: {len(poscars)}")

    xml_paths_abs: list[Path] = []
    for count, poscar in enumerate(poscars, start=1):
        idx = natural_int_from_name(poscar)
        scf_dir = workdir / "fc3" / f"SCF-{idx:05d}"
        xml_path = scf_dir / "vasprun.xml"
        evaluate_poscar_force(
            poscar=poscar,
            calc=calc,
            xml_out=xml_path,
            copy_poscar_to=scf_dir / "POSCAR",
            forces_out=scf_dir / "forces.dat",
        )
        xml_paths_abs.append(xml_path)
        if verbose:
            print(f"[FC3] {count:5d}/{len(poscars):5d} {poscar.name} -> {xml_path}")

    # thirdorder_vasp.py reap reads a list of vasprun.xml paths from stdin.
    # Keep paths relative to workdir for portability.
    xml_paths_rel = [str(p.relative_to(workdir)) for p in natural_sort(xml_paths_abs)]
    input_text = "\n".join(xml_paths_rel) + "\n"

    cmd = [thirdorder_exe, "reap"] + dims + [str(fc3_cutoff)]
    run_checked(
        cmd,
        cwd=workdir,
        input_text=input_text,
        log_path=workdir / "thirdorder_reap.log",
    )

    fc3 = workdir / "FORCE_CONSTANTS_3RD"
    if not fc3.is_file():
        raise RuntimeError(
            "thirdorder_vasp.py reap finished but FORCE_CONSTANTS_3RD was not created."
        )

    if verbose:
        print(f"IFC3 written: {fc3}")

    return fc3




def _validate_fc2_text_header_for_shengbte(
    fc2_file: Path,
    poscar: Path,
    fc2_dim: Sequence[int],
) -> None:
    """Check FORCE_CONSTANTS_2ND first line matches natoms(unitcell)*prod(scell)."""
    fc2_file = Path(fc2_file)
    poscar = Path(poscar)
    atoms = read(str(poscar), format="vasp")
    natoms_unit = len(atoms)
    expected = natoms_unit * int(np.prod([int(x) for x in fc2_dim]))

    if not fc2_file.is_file():
        raise RuntimeError(f"Missing FC2 file: {fc2_file}")

    first = fc2_file.read_text(encoding="utf-8", errors="ignore").splitlines()[0].split()
    if not first:
        raise RuntimeError(f"{fc2_file} has an empty first line")

    try:
        got = int(first[0])
    except Exception as exc:
        raise RuntimeError(f"Cannot parse first line of {fc2_file!s}: {first}") from exc

    if got != expected:
        raise RuntimeError(
            f"{fc2_file} first line is {got}, but ShengBTE CONTROL expects "
            f"natoms*scell = {natoms_unit}*{tuple(int(x) for x in fc2_dim)} = {expected}.\n"
            "This usually means phonopy wrote compact FC2. Re-run with phonopy --full-fc."
        )



# =============================================================================
# OUTCAR parser for NAC
# =============================================================================

_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?")


def _extract_floats(line: str) -> list[float]:
    return [float(x) for x in _FLOAT_RE.findall(line)]


def _find_last_index_contains(lines: list[str], keyword: str) -> int:
    for i in range(len(lines) - 1, -1, -1):
        if keyword in lines[i]:
            return i
    return -1


def parse_epsilon_from_outcar(outcar: Path) -> np.ndarray:
    lines = Path(outcar).read_text(encoding="utf-8", errors="ignore").splitlines()
    idx = _find_last_index_contains(lines, "MACROSCOPIC STATIC DIELECTRIC TENSOR")
    if idx < 0:
        raise RuntimeError(f"Could not find dielectric tensor in OUTCAR: {outcar}")

    rows = []
    for line in lines[idx + 1 : idx + 30]:
        vals = _extract_floats(line)
        if len(vals) >= 3:
            rows.append(vals[:3])
            if len(rows) == 3:
                break

    eps = np.asarray(rows, dtype=float)
    if eps.shape != (3, 3):
        raise RuntimeError(f"Failed to parse dielectric tensor from {outcar}, got {eps.shape}")
    return eps


def parse_born_from_outcar(outcar: Path, natoms: int) -> np.ndarray:
    lines = Path(outcar).read_text(encoding="utf-8", errors="ignore").splitlines()
    idx = _find_last_index_contains(lines, "BORN EFFECTIVE CHARGES")
    if idx < 0:
        raise RuntimeError(f"Could not find Born effective charges in OUTCAR: {outcar}")

    born_tensors = []
    i = idx + 1
    while i < len(lines):
        s = lines[i].strip().lower()
        if s.startswith("ion"):
            rows = []
            j = i + 1
            while j < len(lines) and len(rows) < 3:
                vals = _extract_floats(lines[j])
                if len(vals) >= 4:
                    row_index = int(round(vals[0]))
                    if row_index in (1, 2, 3):
                        rows.append(vals[-3:])
                elif rows:
                    break
                j += 1

            if len(rows) == 3:
                born_tensors.append(rows)
                if len(born_tensors) == natoms:
                    break
                i = j
                continue
        i += 1

    born = np.asarray(born_tensors, dtype=float)
    if born.shape != (natoms, 3, 3):
        raise RuntimeError(
            f"Failed to parse Born charges from {outcar}; got {born.shape}, "
            f"expected ({natoms}, 3, 3)."
        )
    return born


def read_nac_from_outcar(outcar: Path, natoms: int) -> tuple[np.ndarray, np.ndarray]:
    outcar = Path(outcar)
    if not outcar.is_file():
        raise FileNotFoundError(f"NAC requested but OUTCAR not found: {outcar}")
    epsilon = parse_epsilon_from_outcar(outcar)
    born = parse_born_from_outcar(outcar, natoms=natoms)
    return epsilon, born


# =============================================================================
# ShengBTE CONTROL writer
# =============================================================================

def write_shengbte_control(
    filename: Path,
    poscar: Path,
    qgrid: Sequence[int],
    fc2_dim: Sequence[int],
    tmin: float,
    tmax: float,
    tstep: float,
    scalebroad: float = 0.1,
    iterative: bool = True,
    diis: bool = True,
    itergpu: bool = True,
    nanowires: bool = False,
    lfactor_nm: float = 0.1,
    nonanalytic: bool = False,
    epsilon: np.ndarray | None = None,
    born: np.ndarray | None = None,
) -> Path:
    """
    Write ShengBTE CONTROL file from POSCAR.

    lfactor=0.1 means lattice vectors are supplied in Angstrom.
    iterative=True maps to convergence=.TRUE.
    scell(:) is the FC2 supercell used to generate FORCE_CONSTANTS_2ND.
    """
    atoms = read(str(poscar), format="vasp")
    atoms.pbc = True

    qgrid = [int(x) for x in qgrid]
    scell = [int(x) for x in fc2_dim]
    if len(qgrid) != 3 or any(x <= 0 for x in qgrid):
        raise ValueError(f"qgrid must be 3 positive integers, got {qgrid}")
    if len(scell) != 3 or any(x <= 0 for x in scell):
        raise ValueError(f"fc2_dim/scell must be 3 positive integers, got {scell}")

    symbols = atoms.get_chemical_symbols()
    unique_elems: list[str] = []
    for sym in symbols:
        if sym not in unique_elems:
            unique_elems.append(sym)
    elem_to_type = {sym: i + 1 for i, sym in enumerate(unique_elems)}
    types = [elem_to_type[sym] for sym in symbols]

    frac = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
    frac %= 1.0
    frac[np.isclose(frac, 1.0)] = 0.0
    cell_ang = np.asarray(atoms.cell.array, dtype=float)

    eps_arr = None
    born_arr = None
    if nonanalytic:
        if epsilon is None or born is None:
            raise ValueError("nonanalytic=True requires epsilon and born.")
        eps_arr = np.asarray(epsilon, dtype=float)
        born_arr = np.asarray(born, dtype=float)
        if eps_arr.shape != (3, 3):
            raise ValueError(f"epsilon must be (3,3), got {eps_arr.shape}")
        if born_arr.shape != (len(atoms), 3, 3):
            raise ValueError(
                f"born must be ({len(atoms)},3,3), got {born_arr.shape}"
            )

    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)

    with filename.open("w", encoding="utf-8") as fh:
        fh.write("&allocations\n")
        fh.write(f"    nelements={len(unique_elems)}\n")
        fh.write(f"    natoms={len(atoms)}\n")
        fh.write(f"    ngrid(:)={qgrid[0]} {qgrid[1]} {qgrid[2]}\n")
        fh.write("&end\n\n")

        fh.write("&crystal\n")
        fh.write(f"    lfactor={lfactor_nm:.10f}\n")
        for i in range(3):
            fh.write(
                f"    lattvec(:,{i + 1})="
                f"{cell_ang[i, 0]:.12f} {cell_ang[i, 1]:.12f} {cell_ang[i, 2]:.12f}\n"
            )

        elements_str = " ".join(f'"{elem}"' for elem in unique_elems)
        fh.write(f"    elements={elements_str}\n")
        fh.write("    types=" + " ".join(str(t) for t in types) + "\n")

        for i, pos in enumerate(frac, start=1):
            fh.write(
                f"    positions(:,{i})="
                f"{pos[0]:.12f} {pos[1]:.12f} {pos[2]:.12f}\n"
            )

        if eps_arr is not None:
            for i in range(3):
                fh.write(
                    f"    epsilon(:,{i + 1})="
                    f"{eps_arr[0, i]:.12f} {eps_arr[1, i]:.12f} {eps_arr[2, i]:.12f}\n"
                )

        if born_arr is not None:
            for atom_idx in range(len(atoms)):
                for i in range(3):
                    fh.write(
                        f"    born(:,{i + 1},{atom_idx + 1})="
                        f"{born_arr[atom_idx, 0, i]:.12f} "
                        f"{born_arr[atom_idx, 1, i]:.12f} "
                        f"{born_arr[atom_idx, 2, i]:.12f}\n"
                    )

        fh.write(f"    scell(:)={scell[0]} {scell[1]} {scell[2]}\n")
        fh.write("&end\n\n")

        fh.write("&parameters\n")
        fh.write(f"    T_min={float(tmin):.10g}\n")
        fh.write(f"    T_max={float(tmax):.10g}\n")
        fh.write(f"    T_step={float(tstep):.10g}\n")
        fh.write(f"    scalebroad={float(scalebroad):.10g}\n")
        fh.write("&end\n\n")

        fh.write("&flags\n")
        fh.write(f"    nonanalytic={'.TRUE.' if nonanalytic else '.FALSE.'}\n")
        fh.write(f"    nanowires={'.TRUE.' if nanowires else '.FALSE.'}\n")
        fh.write(f"    convergence={'.TRUE.' if iterative else '.FALSE.'}\n")
        fh.write(f"    itergpu={'.TRUE.' if itergpu else '.FALSE.'}\n")
        fh.write(f"    DIIS={'.TRUE.' if diis else '.FALSE.'}\n")
        fh.write("&end\n")

    return filename


def _validate_fc2_text_header_for_shengbte(
    fc2_file: Path,
    poscar: Path,
    fc2_dim: Sequence[int],
) -> None:
    """Check FORCE_CONSTANTS_2ND first line matches natoms(unitcell)*prod(scell)."""
    fc2_file = Path(fc2_file)
    poscar = Path(poscar)
    atoms = read(str(poscar), format="vasp")
    natoms_unit = len(atoms)
    expected = natoms_unit * int(np.prod([int(x) for x in fc2_dim]))

    if not fc2_file.is_file():
        raise RuntimeError(f"Missing FC2 file: {fc2_file}")

    first = fc2_file.read_text(encoding="utf-8", errors="ignore").splitlines()[0].split()
    if not first:
        raise RuntimeError(f"{fc2_file} has an empty first line")
    got = int(first[0])
    if got != expected:
        raise RuntimeError(
            f"{fc2_file} first line is {got}, but ShengBTE CONTROL expects "
            f"natoms*scell = {natoms_unit}*{tuple(int(x) for x in fc2_dim)} = {expected}.\n"
            "This usually means phonopy wrote compact FC2. Re-run with phonopy --full-fc."
        )


# =============================================================================
# Batch / cleanup helpers
# =============================================================================

KEEP_FINAL_FILES = {"CONTROL", "FORCE_CONSTANTS_2ND", "FORCE_CONSTANTS_3RD", "NNP_BTE_REPORT.yaml", "NNP_BTE_REPORT.txt"}


def cleanup_workdir_keep_only_final(workdir: Path, verbose: bool = True) -> None:
    """Delete intermediates and keep only CONTROL/2ND/3RD plus NNP_BTE_REPORT."""
    workdir = Path(workdir)
    for item in list(workdir.iterdir()):
        if item.name in KEEP_FINAL_FILES:
            continue
        if item.is_dir():
            shutil.rmtree(item)
        else:
            item.unlink()
    if verbose:
        kept = ", ".join(sorted(KEEP_FINAL_FILES))
        print(f"=> Cleaned intermediates in {workdir}; kept only: {kept}")


def find_case_directories(
    root: Path,
    recursive: bool = False,
    poscar_name: str = "POSCAR",
    allow_contcar: bool = True,
) -> list[Path]:
    root = Path(root)
    pattern_dirs = root.rglob("*") if recursive else root.iterdir()
    cases: list[Path] = []

    for p in pattern_dirs:
        if not p.is_dir():
            continue
        name = p.name
        if name.endswith("-BTE") or name in {"__pycache__", "ifc_nnp", "BTE"}:
            continue
        if (p / poscar_name).is_file() or (allow_contcar and (p / "CONTCAR").is_file()):
            cases.append(p)

    return sorted(cases, key=lambda x: str(x))


def choose_structure_file(case_dir: Path, poscar_name: str, allow_contcar: bool) -> Path:
    """Choose the structure file for a case.

    POSCAR is preferred. If POSCAR is missing and CONTCAR fallback is allowed,
    emit a visible warning before using CONTCAR. This avoids silent benchmark
    mismatches when a raw POSCAR was expected.
    """
    p = case_dir / poscar_name
    if p.is_file():
        return p

    c = case_dir / "CONTCAR"
    if allow_contcar and c.is_file():
        print(
            f"WARNING: {case_dir} has no {poscar_name}; falling back to CONTCAR: {c}",
            file=sys.stderr,
            flush=True,
        )
        return c

    raise FileNotFoundError(f"No {poscar_name} or CONTCAR found in {case_dir}")


def case_output_dir(case_dir: Path, output_style: str, outdir_name: str) -> Path:
    if output_style == "sibling-bte":
        return case_dir.parent / f"{case_dir.name}-BTE"
    if output_style == "inside":
        return case_dir / outdir_name
    raise ValueError(f"Unknown output_style: {output_style}")


def generate_single_case(
    *,
    poscar: Path,
    workdir: Path,
    calc,
    args,
    verbose: bool,
) -> list[Path]:
    """
    Generate relaxed structure -> FC2 -> FC2 soft-mode report -> FC3 -> CONTROL.

    Important design choice for high-throughput screening:
    harmonic-stability labels are written to NNP_BTE_REPORT, but the workflow
    continues to FC3 and CONTROL regardless of imaginary frequencies.
    """
    poscar = Path(poscar).resolve()
    workdir = Path(workdir).resolve()

    if workdir.exists():
        if args.overwrite:
            shutil.rmtree(workdir)
        else:
            raise FileExistsError(
                f"Output directory already exists: {workdir}. Use --overwrite."
            )
    workdir.mkdir(parents=True, exist_ok=True)

    report: dict = {
        "schema_version": 1,
        "workflow_status": "started",
        "created_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "case_name": workdir.name,
        "workdir": str(workdir),
        "source_poscar": str(poscar),
        "parameters": {
            "fc2_dim": [int(x) for x in args.fc2_dim],
            "fc3_dim": [int(x) for x in args.fc3_dim],
            "fc2_distance_A": float(args.fc2_distance),
            "fc3_cutoff": str(args.fc3_cutoff),
            "qgrid": [int(x) for x in args.qgrid],
            "nac": bool(args.nac),
            "clean_intermediates": bool(args.clean_intermediates),
        },
        "outputs": {},
    }

    outputs: list[Path] = []

    try:
        # Keep a copy until cleanup; CONTROL contains the final relaxed structure.
        shutil.copy2(poscar, workdir / "POSCAR.input")

        if args.relax:
            report["relaxation"] = relax_structure_with_ase(
                input_poscar=workdir / "POSCAR.input",
                output_poscar=workdir / "POSCAR",
                calc=calc,
                fmax=args.relax_fmax,
                steps=args.relax_steps,
                optimizer_name=args.relax_optimizer,
                relax_cell=args.relax_cell,
                constant_volume=args.relax_constant_volume,
                verbose=verbose,
            )
        else:
            report["relaxation"] = no_relax_structure_copy(
                input_poscar=workdir / "POSCAR.input",
                output_poscar=workdir / "POSCAR",
                calc=calc,
            )

        write_case_reports(workdir, report)

        fc2_path = generate_ifc2_with_phonopy_vasp_like(
            workdir=workdir,
            calc=calc,
            phonopy_exe=args.phonopy,
            fc2_dim=args.fc2_dim,
            fc2_distance=args.fc2_distance,
            verbose=verbose,
        )
        outputs.append(fc2_path)
        report["outputs"]["FORCE_CONSTANTS_2ND"] = str(fc2_path)
        report["fc2_generation"] = {"status": "ok", "path": str(fc2_path)}
        write_case_reports(workdir, report)

        report["harmonic_stability"] = analyze_fc2_harmonic_stability(
            workdir=workdir,
            phonopy_exe=args.phonopy,
            fc2_dim=args.fc2_dim,
            mesh_dims=args.stability_mesh,
            imag_tol=args.imag_tol,
            soft_mode_threshold=args.soft_mode_threshold,
            optical_ratio_threshold=args.optical_ratio_threshold,
            soft_branch_qfraction=args.soft_branch_qfraction,
            many_imag_fraction=args.many_imag_fraction,
            verbose=verbose,
        )
        write_case_reports(workdir, report)

        # Always generate FC3 even if FC2 has imaginary modes.
        fc3_path = generate_ifc3_with_thirdorder_vasp_like(
            workdir=workdir,
            calc=calc,
            thirdorder_exe=args.thirdorder,
            fc3_dim=args.fc3_dim,
            fc3_cutoff=str(args.fc3_cutoff),
            verbose=verbose,
        )
        outputs.append(fc3_path)
        report["outputs"]["FORCE_CONSTANTS_3RD"] = str(fc3_path)
        report["fc3_generation"] = {"status": "ok", "path": str(fc3_path)}
        write_case_reports(workdir, report)

        nac_enabled = bool(args.nac)
        epsilon = None
        born = None
        if nac_enabled:
            atoms = read(str(workdir / "POSCAR"), format="vasp")
            outcar = poscar.parent / args.outcar_name
            epsilon, born = read_nac_from_outcar(outcar, natoms=len(atoms))
            if verbose:
                print(f"=> NAC data read from {outcar}")

        control_path = write_shengbte_control(
            filename=workdir / "CONTROL",
            poscar=workdir / "POSCAR",
            qgrid=args.qgrid,
            fc2_dim=args.fc2_dim,
            tmin=args.tmin,
            tmax=args.tmax,
            tstep=args.tstep,
            scalebroad=args.scalebroad,
            iterative=not args.rta,
            diis=not args.no_diis,
            itergpu=not args.no_itergpu,
            nanowires=args.nanowires,
            nonanalytic=nac_enabled,
            epsilon=epsilon,
            born=born,
        )
        outputs.append(control_path)
        report["outputs"]["CONTROL"] = str(control_path)
        report["control_generation"] = {"status": "ok", "path": str(control_path)}

        # Final sanity check for FC2 after CONTROL is written.
        _validate_fc2_text_header_for_shengbte(
            fc2_file=workdir / "FORCE_CONSTANTS_2ND",
            poscar=workdir / "POSCAR",
            fc2_dim=args.fc2_dim,
        )

        report["workflow_status"] = "success"
        report["finished_at"] = datetime.datetime.now().isoformat(timespec="seconds")
        write_case_reports(workdir, report)

        if args.clean_intermediates:
            cleanup_workdir_keep_only_final(workdir, verbose=verbose)

        return outputs

    except Exception as exc:
        report["workflow_status"] = "failed"
        report["error"] = repr(exc)
        report["finished_at"] = datetime.datetime.now().isoformat(timespec="seconds")
        try:
            write_case_reports(workdir, report)
        except Exception:
            pass
        raise


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Batch/single-case VASP-like IFC2/IFC3 generation using ASE/NNP forces, "
            "phonopy, thirdorder_vasp.py, and ShengBTE CONTROL writer."
        )
    )

    p.add_argument("--poscar", default="POSCAR", help="Single-case POSCAR path.")
    p.add_argument(
        "--calc-factory",
        default=None,
        help=(
            "Optional override, e.g. old-style calc_factory.py:make_calc. "
            "If omitted, the internal make_calc() near the top of this script is used."
        ),
    )
    p.add_argument(
        "--model-path",
        default=None,
        help="Optional override for NNP_MODEL_PATH used by internal make_calc().",
    )
    p.add_argument(
        "--device",
        default=None,
        help="Optional override for NNP_DEVICE used by internal make_calc(), e.g. cuda or cpu.",
    )
    p.add_argument("--workdir", default="ifc_nnp", help="Single-case output directory.")

    p.add_argument("--batch-root", default=None, help="Parent directory containing many case folders.")
    p.add_argument("--recursive", action="store_true", help="Recursively search batch-root.")
    p.add_argument("--poscar-name", default="POSCAR", help="Structure filename to look for in each case.")
    p.add_argument("--no-contcar", action="store_true", help="Do not fall back to CONTCAR.")
    p.add_argument(
        "--output-style",
        choices=["sibling-bte", "inside"],
        default="sibling-bte",
        help="Batch output location: sibling-bte -> case-BTE beside case; inside -> case/outdir-name.",
    )
    p.add_argument("--outdir-name", default="ifc_nnp", help="Used when --output-style inside.")
    p.add_argument("--continue-on-error", action="store_true", help="Continue batch after a case fails.")
    p.add_argument("--parallel-jobs", type=int, default=1,
                   help="Number of independent batch workers. Use 2 for two GPUs. Default: 1.")
    p.add_argument("--gpu-devices", default=None,
                   help=(
                       "Comma-separated GPU ids assigned to workers, e.g. '0,1'. "
                       "If omitted, uses CUDA_VISIBLE_DEVICES when set, otherwise 0..parallel_jobs-1. "
                       "Each worker gets one CUDA_VISIBLE_DEVICES value and processes its own serial queue."
                   ))
    p.add_argument("--mp-start-method", choices=["spawn", "forkserver", "fork"], default="spawn",
                   help="Multiprocessing start method for parallel batch mode. Default: spawn.")

    p.add_argument("--relax", action=argparse.BooleanOptionalAction, default=True,
                   help="NNP-relax each structure before IFC generation. Default: true.")
    p.add_argument("--relax-cell", action=argparse.BooleanOptionalAction, default=True,
                   help="Relax cell+ions when stress is available; otherwise fallback to ions-only. Default: true.")
    p.add_argument("--relax-constant-volume", action="store_true",
                   help="When relaxing cell, keep volume approximately constant through the ASE cell filter.")
    p.add_argument("--relax-fmax", type=float, default=1.0e-3,
                   help="ASE relaxation fmax in eV/A. Default: 1e-3.")
    p.add_argument("--relax-steps", type=int, default=500,
                   help="Maximum ASE relaxation steps. Default: 500.")
    p.add_argument("--relax-optimizer", choices=["fire", "bfgs", "lbfgs"], default="fire",
                   help="ASE optimizer for NNP relaxation. Default: fire.")

    p.add_argument("--fc2-dim", nargs=3, type=int, default=[2, 2, 1])
    p.add_argument("--fc3-dim", nargs=3, type=int, default=[2, 2, 1])
    p.add_argument("--fc2-distance", type=float, default=0.01)
    p.add_argument("--fc3-cutoff", default="-3", help="-3 means third nearest neighbor in thirdorder_vasp.py.")

    p.add_argument("--stability-mesh", nargs=3, type=int, default=[12, 12, 12],
                   help="Phonopy mesh for FC2 imaginary-frequency report. Default: 12 12 12.")
    p.add_argument("--imag-tol", type=float, default=0.05,
                   help="Imaginary-frequency tolerance in THz; frequencies above -imag_tol are treated as tiny. Default: 0.05.")
    p.add_argument("--soft-mode-threshold", type=float, default=0.20,
                   help="Strong soft-mode threshold in THz. Default: 0.20.")
    p.add_argument("--optical-ratio-threshold", type=float, default=0.50,
                   help="Gamma translation_ratio threshold above which an imaginary mode is treated as optical-like. Default: 0.50.")
    p.add_argument("--soft-branch-qfraction", type=float, default=0.10,
                   help="Flag a soft branch if one band is significantly imaginary over this weighted q fraction. Default: 0.10.")
    p.add_argument("--many-imag-fraction", type=float, default=0.01,
                   help="Flag many imaginary modes if significant weighted negative modes exceed this fraction of all modes. Default: 0.01.")

    p.add_argument("--qgrid", nargs=3, type=int, default=[12, 12, 12])
    p.add_argument("--tmin", type=float, default=300.0)
    p.add_argument("--tmax", type=float, default=900.0)
    p.add_argument("--tstep", type=float, default=100.0)
    p.add_argument("--scalebroad", type=float, default=0.1)
    p.add_argument("--rta", action="store_true", help="Write convergence=.FALSE.; default is iterative.")
    p.add_argument("--no-diis", action="store_true", help="Write DIIS=.FALSE.")
    p.add_argument("--no-itergpu", action="store_true", help="Write itergpu=.FALSE.")
    p.add_argument("--nanowires", action="store_true")

    p.add_argument("--nac", action="store_true", help="Read OUTCAR next to POSCAR and write NAC into CONTROL.")
    p.add_argument("--outcar-name", default="OUTCAR")

    p.add_argument("--phonopy", default="phonopy")
    p.add_argument("--thirdorder", default="thirdorder_vasp.py")

    p.add_argument("--overwrite", action="store_true", help="Remove existing output directory first.")
    p.add_argument("--clean-intermediates", action=argparse.BooleanOptionalAction, default=True,
                   help="After success, delete intermediates and keep only CONTROL/2ND/3RD plus NNP_BTE_REPORT.")
    p.add_argument("--quiet", action="store_true")

    return p.parse_args()



def _parse_gpu_devices_for_workers(gpu_devices: str | None, parallel_jobs: int) -> list[str | None]:
    """Return one GPU visibility string per worker.

    Examples:
      --gpu-devices 0,1 with --parallel-jobs 2 -> ["0", "1"]
      CUDA_VISIBLE_DEVICES=2,3 and no --gpu-devices -> ["2", "3"]
      no environment and --parallel-jobs 2 -> ["0", "1"]
    """
    n = int(parallel_jobs)
    if n <= 0:
        raise ValueError(f"parallel_jobs must be positive, got {parallel_jobs}")

    if gpu_devices is not None and str(gpu_devices).strip():
        devices = [x.strip() for x in str(gpu_devices).split(",") if x.strip()]
    else:
        visible = os.environ.get("CUDA_VISIBLE_DEVICES", "").strip()
        if visible:
            devices = [x.strip() for x in visible.split(",") if x.strip()]
        else:
            devices = [str(i) for i in range(n)]

    if not devices:
        devices = [None]

    if len(devices) < n:
        print(
            f"WARNING: parallel_jobs={n} but only {len(devices)} GPU device(s) were provided/visible. "
            "Devices will be reused round-robin.",
            file=sys.stderr,
            flush=True,
        )

    return [devices[i % len(devices)] for i in range(n)]


def _split_cases_round_robin(cases: list[Path], nqueues: int) -> list[list[Path]]:
    queues: list[list[Path]] = [[] for _ in range(int(nqueues))]
    for i, case in enumerate(cases):
        queues[i % int(nqueues)].append(case)
    return queues


def _load_calculator_for_worker(args, verbose: bool):
    """Load one calculator inside a worker process, after GPU env is set."""
    global NNP_MODEL_PATH, NNP_DEVICE
    if args.model_path is not None:
        NNP_MODEL_PATH = Path(args.model_path).expanduser().resolve()
    if args.device is not None:
        NNP_DEVICE = str(args.device)

    if verbose:
        print("=== Loading ASE calculator in worker ===", flush=True)
        if args.calc_factory:
            print(f"Calculator source: external {args.calc_factory}", flush=True)
        else:
            print("Calculator source: internal make_calc()", flush=True)
            print(f"CALCULATOR_KIND : {CALCULATOR_KIND}", flush=True)
            if str(CALCULATOR_KIND).lower().strip() == "nequip_compiled":
                print(f"NNP_MODEL_PATH  : {NNP_MODEL_PATH}", flush=True)
                print(f"NNP_DEVICE      : {NNP_DEVICE}", flush=True)

    if args.calc_factory:
        factory = load_calc_factory(args.calc_factory)
        calc = factory()
    else:
        calc = make_calc()

    if verbose:
        print("Calculator loaded in worker.", flush=True)
    return calc


def _batch_worker_process(
    worker_id: int,
    gpu_device: str | None,
    cases: list[Path],
    root: Path,
    args,
    verbose: bool,
) -> dict:
    """Run one serial queue of cases inside one process bound to one GPU."""
    root = Path(root).resolve()
    safe_gpu = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(gpu_device) if gpu_device is not None else "none")
    worker_log = root / f"batch_worker_{worker_id:02d}_gpu_{safe_gpu}.log"

    successes: list[str] = []
    failures: list[tuple[str, str]] = []

    # Redirect worker chatter to a per-worker log. The parent prints compact summaries.
    with worker_log.open("w", encoding="utf-8") as log_fh, contextlib.redirect_stdout(log_fh), contextlib.redirect_stderr(log_fh):
        print(f"Worker {worker_id} starting at {datetime.datetime.now()}", flush=True)
        print(f"Assigned GPU device: {gpu_device}", flush=True)
        print(f"Cases in this queue: {len(cases)}", flush=True)

        if gpu_device is not None:
            # Must happen before the calculator imports torch/jax/cuda internals.
            os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_device)
            print(f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES')}", flush=True)

        calc = _load_calculator_for_worker(args, verbose=verbose)

        for idx, case in enumerate(cases, start=1):
            print(f"\n===== [worker {worker_id} | {idx}/{len(cases)}] {case} =====", flush=True)
            try:
                poscar = choose_structure_file(case, args.poscar_name, allow_contcar=not args.no_contcar)
                workdir = case_output_dir(case, args.output_style, args.outdir_name)
                generate_single_case(
                    poscar=poscar,
                    workdir=workdir,
                    calc=calc,
                    args=args,
                    verbose=verbose,
                )
                successes.append(str(workdir))
                print(f"=> SUCCESS: {workdir}", flush=True)
            except Exception as exc:
                failures.append((str(case), repr(exc)))
                print(f"=> FAILED: {case}: {exc}", flush=True)
                if not args.continue_on_error:
                    raise

        print(f"\nWorker {worker_id} finished at {datetime.datetime.now()}", flush=True)
        print(f"Success: {len(successes)}", flush=True)
        print(f"Failed : {len(failures)}", flush=True)

    return {
        "worker_id": int(worker_id),
        "gpu_device": str(gpu_device) if gpu_device is not None else None,
        "log": str(worker_log),
        "successes": successes,
        "failures": failures,
    }


def run_parallel_batch(args, verbose: bool) -> None:
    """Run batch mode with multiple independent GPU queues."""
    root = Path(args.batch_root).resolve()
    cases = find_case_directories(
        root=root,
        recursive=args.recursive,
        poscar_name=args.poscar_name,
        allow_contcar=not args.no_contcar,
    )
    if not cases:
        raise RuntimeError(f"No case directories found under {root}")

    # Preflight warnings are printed to the main terminal, so CONTCAR fallback is never silent.
    fallback_cases = [
        c for c in cases
        if not (c / args.poscar_name).is_file() and (c / "CONTCAR").is_file() and not args.no_contcar
    ]
    for c in fallback_cases:
        print(
            f"WARNING: {c} has no {args.poscar_name}; this run will read CONTCAR: {c / 'CONTCAR'}",
            file=sys.stderr,
            flush=True,
        )

    nworkers = max(1, int(args.parallel_jobs))
    gpu_assignments = _parse_gpu_devices_for_workers(args.gpu_devices, nworkers)
    queues = _split_cases_round_robin(cases, nworkers)

    print("=== Parallel batch mode ===")
    print(f"Root                : {root}")
    print(f"Cases found         : {len(cases)}")
    print(f"Workers/queues      : {nworkers}")
    print(f"GPU assignments     : {gpu_assignments}")
    print(f"Start method        : {args.mp_start_method}")
    print(f"Output style        : {args.output_style}")
    print(f"Clean intermediates : {args.clean_intermediates}")
    print(f"CONTCAR fallbacks   : {len(fallback_cases)}")
    print(datetime.datetime.now(), flush=True)

    ctx = mp.get_context(args.mp_start_method)
    all_successes: list[str] = []
    all_failures: list[tuple[str, str]] = []
    worker_summaries: list[dict] = []

    with cf.ProcessPoolExecutor(max_workers=nworkers, mp_context=ctx) as executor:
        futures = []
        for wid, queue in enumerate(queues):
            if not queue:
                continue
            futures.append(
                executor.submit(
                    _batch_worker_process,
                    wid,
                    gpu_assignments[wid],
                    queue,
                    root,
                    args,
                    verbose,
                )
            )

        for fut in cf.as_completed(futures):
            summary = fut.result()
            worker_summaries.append(summary)
            all_successes.extend(summary.get("successes", []))
            all_failures.extend(summary.get("failures", []))
            print(
                f"[worker {summary['worker_id']}] GPU={summary['gpu_device']} "
                f"OK={len(summary.get('successes', []))} "
                f"ERR={len(summary.get('failures', []))} "
                f"log={summary['log']}",
                flush=True,
            )

    print("\n=== Parallel batch summary ===")
    print(f"Success: {len(all_successes)}")
    for p in all_successes:
        print(f"  OK  {p}")
    print(f"Failed : {len(all_failures)}")
    for p, exc in all_failures:
        print(f"  ERR {p}: {exc}")
    print("Worker logs:")
    for s in sorted(worker_summaries, key=lambda x: x["worker_id"]):
        print(f"  worker {s['worker_id']}: {s['log']}")

def main() -> None:
    args = parse_args()
    verbose = not args.quiet

    if shutil.which(args.phonopy) is None:
        raise RuntimeError(f"phonopy executable not found: {args.phonopy}")
    if shutil.which(args.thirdorder) is None:
        raise RuntimeError(f"thirdorder_vasp.py executable not found: {args.thirdorder}")

    # Command-line overrides for the internal make_calc() configuration.
    global NNP_MODEL_PATH, NNP_DEVICE
    if args.model_path is not None:
        NNP_MODEL_PATH = Path(args.model_path).expanduser().resolve()
    if args.device is not None:
        NNP_DEVICE = str(args.device)

    if args.batch_root is not None and int(args.parallel_jobs) > 1:
        # In parallel mode, each worker process must load its own calculator after
        # CUDA_VISIBLE_DEVICES is set. Do not construct a calculator in the parent.
        run_parallel_batch(args, verbose=verbose)
        print(datetime.datetime.now())
        return

    if verbose:
        print("=== Loading ASE calculator once ===")
        if args.calc_factory:
            print(f"Calculator source: external {args.calc_factory}")
        else:
            print("Calculator source: internal make_calc()")
            print(f"CALCULATOR_KIND : {CALCULATOR_KIND}")
            if str(CALCULATOR_KIND).lower().strip() == "nequip_compiled":
                print(f"NNP_MODEL_PATH  : {NNP_MODEL_PATH}")
                print(f"NNP_DEVICE      : {NNP_DEVICE}")

    if args.calc_factory:
        factory = load_calc_factory(args.calc_factory)
        calc = factory()
    else:
        calc = make_calc()

    if verbose:
        print("Calculator loaded.")

    successes: list[Path] = []
    failures: list[tuple[Path, Exception]] = []

    if args.batch_root is not None:
        root = Path(args.batch_root).resolve()
        cases = find_case_directories(
            root=root,
            recursive=args.recursive,
            poscar_name=args.poscar_name,
            allow_contcar=not args.no_contcar,
        )
        if not cases:
            raise RuntimeError(f"No case directories found under {root}")

        fallback_cases = [
            c for c in cases
            if not (c / args.poscar_name).is_file() and (c / "CONTCAR").is_file() and not args.no_contcar
        ]
        for c in fallback_cases:
            print(
                f"WARNING: {c} has no {args.poscar_name}; this run will read CONTCAR: {c / 'CONTCAR'}",
                file=sys.stderr,
                flush=True,
            )

        if verbose:
            print(f"=== Batch mode ===")
            print(f"Root             : {root}")
            print(f"Cases found      : {len(cases)}")
            print(f"CONTCAR fallbacks: {len(fallback_cases)}")
            print(f"Output style     : {args.output_style}")
            print(f"Clean intermediates: {args.clean_intermediates}")
            print(f"NAC              : {args.nac}")
        print(datetime.datetime.now())
        for idx, case in enumerate(cases, start=1):
            print(f"\n===== [{idx}/{len(cases)}] {case} =====")
            try:
                poscar = choose_structure_file(case, args.poscar_name, allow_contcar=not args.no_contcar)
                workdir = case_output_dir(case, args.output_style, args.outdir_name)
                generate_single_case(
                    poscar=poscar,
                    workdir=workdir,
                    calc=calc,
                    args=args,
                    verbose=verbose,
                )
                successes.append(workdir)
                print(f"=> SUCCESS: {workdir}")
            except Exception as exc:
                failures.append((case, exc))
                print(f"=> FAILED: {case}: {exc}")
                if not args.continue_on_error:
                    raise

        print("\n=== Batch summary ===")
        print(f"Success: {len(successes)}")
        for p in successes:
            print(f"  OK  {p}")
        print(f"Failed : {len(failures)}")
        for p, exc in failures:
            print(f"  ERR {p}: {exc}")

    else:
        poscar = Path(args.poscar).resolve()
        workdir = Path(args.workdir).resolve()
        generate_single_case(
            poscar=poscar,
            workdir=workdir,
            calc=calc,
            args=args,
            verbose=verbose,
        )
        print(f"\n=> SUCCESS: {workdir}")
    print(datetime.datetime.now())

if __name__ == "__main__":
    main()
