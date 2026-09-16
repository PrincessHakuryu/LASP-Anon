#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Collect BTE.kappa_tensor from all *-BTE/T**K folders and append FC2 stability info.

Output columns include:
  fc2_stability: strictly green / yellow / red from report severity; otherwise unknown
  fc2_min_frequency_THz: minimum phonon frequency from NNP_BTE_REPORT.yaml if available

The tensor format assumed is ShengBTE-like:
  iter kxx kxy kxz kyx kyy kyz kzx kzy kzz
so kxx, kyy, kzz are tensor components 1, 5, 9 after the iteration column.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Any


def _try_import_yaml():
    try:
        import yaml  # type: ignore
        return yaml
    except Exception:
        return None


def parse_temperature_from_dirname(name: str):
    m = re.search(r"T[^0-9+-]*([0-9]+(?:\.[0-9]+)?)\s*K$", name, re.IGNORECASE)
    if not m:
        m = re.search(r"([0-9]+(?:\.[0-9]+)?)", name)
    return float(m.group(1)) if m else ""


def read_last_numeric_line(path: Path, min_cols: int):
    last_raw = None
    last_vals = None
    with path.open("r", errors="ignore") as f:
        for raw in f:
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.replace("D", "E").replace("d", "E").split()
            vals = []
            ok = True
            for x in parts:
                try:
                    vals.append(float(x))
                except ValueError:
                    ok = False
                    break
            if ok and len(vals) >= min_cols:
                last_raw = s
                last_vals = vals
    if last_vals is None:
        raise ValueError(f"No numeric line with at least {min_cols} columns found in {path}")
    return last_raw, last_vals


def _dig(d: Any, keys: list[str], default=""):
    cur = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def normalize_fc2_color(severity: str = "", label: str = "") -> str:
    """
    Strictly use the report's classification.severity value.

    Only exact green / yellow / red severities are accepted.
    Old or non-standard values such as green_yellow are not remapped.
    """
    sev = str(severity or "").strip().lower().replace("-", "_")
    if sev in {"green", "yellow", "red"}:
        return sev
    return "unknown"

def read_fc2_report(bte_dir: Path):
    """Return (color, min_frequency_THz). Missing reports become (unknown, '')."""
    yml = bte_dir / "NNP_BTE_REPORT.yaml"
    txt = bte_dir / "NNP_BTE_REPORT.txt"

    if yml.is_file():
        yaml = _try_import_yaml()
        if yaml is not None:
            try:
                data = yaml.safe_load(yml.read_text(errors="ignore")) or {}
                cls = _dig(data, ["harmonic_stability", "classification"], {})
                if not isinstance(cls, dict):
                    cls = {}
                severity = cls.get("severity", "")
                label = cls.get("label", "")
                minf = _dig(data, ["harmonic_stability", "mesh_summary", "min_frequency_THz"], "")
                # Some external check scripts use a flatter key layout.
                if minf == "":
                    minf = _dig(data, ["mesh_raw", "min_frequency_THz"], "")
                return normalize_fc2_color(severity, label), minf
            except Exception:
                pass

        # Lightweight regex fallback if PyYAML is unavailable or parsing failed.
        text = yml.read_text(errors="ignore")
    elif txt.is_file():
        text = txt.read_text(errors="ignore")
    else:
        return "unknown", ""

    sev_match = re.search(r"^\s*severity\s*:\s*([^\n#]+)", text, re.MULTILINE)
    lab_match = re.search(r"^\s*label\s*:\s*([^\n#]+)", text, re.MULTILINE)
    min_match = re.search(r"^\s*min_frequency_THz\s*:\s*([-+0-9.eEdD]+)", text, re.MULTILINE)
    severity = sev_match.group(1).strip() if sev_match else ""
    label = lab_match.group(1).strip() if lab_match else ""
    minf = min_match.group(1).replace("D", "E").replace("d", "E") if min_match else ""
    return normalize_fc2_color(severity, label), minf


def main():
    ap = argparse.ArgumentParser(description="Collect xx/yy/zz kappa from *-BTE/T**K/BTE.kappa_tensor with FC2 stability columns.")
    ap.add_argument("--root", default=".", help="Root directory. Default: current directory.")
    ap.add_argument("--out", default="kappa_tensor_all_T_with_fc2.csv", help="Output CSV filename.")
    ap.add_argument("--pattern", default="*-BTE", help="BTE directory glob pattern. Default: *-BTE")
    ap.add_argument("--include-missing", action="store_true", help="Include missing/error rows in CSV.")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    bte_dirs = sorted([p for p in root.glob(args.pattern) if p.is_dir()], key=lambda p: p.name)
    rows = []
    missing = []

    for bte_dir in bte_dirs:
        fc2_color, fc2_minf = read_fc2_report(bte_dir)
        temp_dirs = sorted(
            [p for p in bte_dir.iterdir() if p.is_dir() and re.match(r"^T.*K$", p.name, re.IGNORECASE)],
            key=lambda p: (parse_temperature_from_dirname(p.name) if parse_temperature_from_dirname(p.name) != "" else 1e99, p.name),
        )
        if not temp_dirs:
            missing.append(f"{bte_dir}: no T**K directories")
            if args.include_missing:
                rows.append({
                    "case": bte_dir.name, "temperature_dir": "", "temperature_K": "", "iteration": "",
                    "kappa_xx": "", "kappa_yy": "", "kappa_zz": "", "kappa_avg_diag": "",
                    "fc2_stability": fc2_color, "fc2_min_frequency_THz": fc2_minf,
                    "file": "", "status": "NO_T_DIR",
                })
            continue

        for tdir in temp_dirs:
            f = tdir / "BTE.kappa_tensor"
            temp_K = parse_temperature_from_dirname(tdir.name)
            if not f.is_file():
                missing.append(f"{tdir}: missing BTE.kappa_tensor")
                if args.include_missing:
                    rows.append({
                        "case": bte_dir.name, "temperature_dir": tdir.name, "temperature_K": temp_K, "iteration": "",
                        "kappa_xx": "", "kappa_yy": "", "kappa_zz": "", "kappa_avg_diag": "",
                        "fc2_stability": fc2_color, "fc2_min_frequency_THz": fc2_minf,
                        "file": str(f), "status": "MISSING",
                    })
                continue
            try:
                raw, vals = read_last_numeric_line(f, min_cols=9)
                if len(vals) >= 10:
                    iteration = raw.split()[0]
                    tensor = vals[1:10]
                elif len(vals) == 9:
                    iteration = ""
                    tensor = vals[0:9]
                else:
                    raise ValueError(f"Too few numeric columns: {len(vals)}")
                kxx, kyy, kzz = tensor[0], tensor[4], tensor[8]
                rows.append({
                    "case": bte_dir.name,
                    "temperature_dir": tdir.name,
                    "temperature_K": temp_K,
                    "iteration": iteration,
                    "kappa_xx": f"{kxx:.12g}",
                    "kappa_yy": f"{kyy:.12g}",
                    "kappa_zz": f"{kzz:.12g}",
                    "kappa_avg_diag": f"{(kxx + kyy + kzz) / 3.0:.12g}",
                    "fc2_stability": fc2_color,
                    "fc2_min_frequency_THz": fc2_minf,
                    "file": str(f),
                    "status": "OK",
                })
            except Exception as exc:
                missing.append(f"{f}: {exc}")
                if args.include_missing:
                    rows.append({
                        "case": bte_dir.name, "temperature_dir": tdir.name, "temperature_K": temp_K, "iteration": "",
                        "kappa_xx": "", "kappa_yy": "", "kappa_zz": "", "kappa_avg_diag": "",
                        "fc2_stability": fc2_color, "fc2_min_frequency_THz": fc2_minf,
                        "file": str(f), "status": f"ERROR: {exc}",
                    })

    columns = [
        "case", "temperature_dir", "temperature_K", "iteration",
        "kappa_xx", "kappa_yy", "kappa_zz", "kappa_avg_diag",
        "fc2_stability", "fc2_min_frequency_THz", "file", "status",
    ]
    out = Path(args.out)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns)
        w.writeheader()
        w.writerows(rows)

    miss_file = out.with_suffix(".missing.txt")
    miss_file.write_text("\n".join(missing) + ("\n" if missing else ""))

    print(f"Found BTE dirs : {len(bte_dirs)}")
    print(f"OK rows        : {sum(1 for r in rows if r.get('status') == 'OK')}")
    print(f"CSV written    : {out}")
    print(f"Missing/errors : {miss_file}")


if __name__ == "__main__":
    main()
