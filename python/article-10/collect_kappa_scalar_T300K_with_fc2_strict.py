#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Collect final scalar kappa from each folder/T300K/BTE.kappa_scalar and append FC2 stability info.

FC2 stability is read strictly from report classification.severity: green/yellow/red only.
Non-standard old values such as green_yellow become unknown.
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


def read_last_numeric_line(path: Path, min_cols: int = 1):
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
        raise ValueError(f"No numeric line found in {path}")
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

def read_fc2_report(folder: Path):
    yml = folder / "NNP_BTE_REPORT.yaml"
    txt = folder / "NNP_BTE_REPORT.txt"

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
                if minf == "":
                    minf = _dig(data, ["mesh_raw", "min_frequency_THz"], "")
                return normalize_fc2_color(severity, label), minf
            except Exception:
                pass
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
    ap = argparse.ArgumentParser(description="Collect final scalar kappa from each folder/T300K/BTE.kappa_scalar with FC2 stability columns.")
    ap.add_argument("--root", default=".", help="Root directory. Default: current directory.")
    ap.add_argument("--temp", default="T300K", help="Temperature directory. Default: T300K")
    ap.add_argument("--pattern", default="*", help="Folder glob pattern. Default: * ; use '*-BTE' for only BTE folders.")
    ap.add_argument("--out", default="kappa_scalar_T300K_with_fc2.csv", help="Output CSV filename.")
    ap.add_argument("--include-missing", action="store_true", help="Include missing/error rows in CSV.")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    dirs = sorted([p for p in root.glob(args.pattern) if p.is_dir()], key=lambda p: p.name)
    rows = []
    missing = []

    for d in dirs:
        fc2_color, fc2_minf = read_fc2_report(d)
        f = d / args.temp / "BTE.kappa_scalar"
        if not f.is_file():
            missing.append(f"{d}: missing {args.temp}/BTE.kappa_scalar")
            if args.include_missing:
                rows.append({
                    "case": d.name, "temperature_dir": args.temp, "iteration": "", "kappa_scalar": "",
                    "fc2_stability": fc2_color, "fc2_min_frequency_THz": fc2_minf,
                    "file": str(f), "status": "MISSING",
                })
            continue
        try:
            raw, vals = read_last_numeric_line(f, min_cols=1)
            if len(vals) >= 2:
                iteration = raw.split()[0]
                scalar = vals[1]
            else:
                iteration = ""
                scalar = vals[0]
            rows.append({
                "case": d.name,
                "temperature_dir": args.temp,
                "iteration": iteration,
                "kappa_scalar": f"{scalar:.12g}",
                "fc2_stability": fc2_color,
                "fc2_min_frequency_THz": fc2_minf,
                "file": str(f),
                "status": "OK",
            })
        except Exception as exc:
            missing.append(f"{f}: {exc}")
            if args.include_missing:
                rows.append({
                    "case": d.name, "temperature_dir": args.temp, "iteration": "", "kappa_scalar": "",
                    "fc2_stability": fc2_color, "fc2_min_frequency_THz": fc2_minf,
                    "file": str(f), "status": f"ERROR: {exc}",
                })

    columns = ["case", "temperature_dir", "iteration", "kappa_scalar", "fc2_stability", "fc2_min_frequency_THz", "file", "status"]
    out = Path(args.out)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns)
        w.writeheader()
        w.writerows(rows)

    miss_file = out.with_suffix(".missing.txt")
    miss_file.write_text("\n".join(missing) + ("\n" if missing else ""))

    print(f"Scanned folders : {len(dirs)}")
    print(f"OK rows         : {sum(1 for r in rows if r.get('status') == 'OK')}")
    print(f"CSV written     : {out}")
    print(f"Missing/errors  : {miss_file}")


if __name__ == "__main__":
    main()
