#!/usr/bin/env python3
"""
Usage:
  python master_script.py

Notes (NYU Greene / SLURM):
- This script is meant to be run on a login node (it only submits/monitors jobs and edits small files).
- ORCA work still runs inside your SLURM job script: orca.sbatch
"""
from __future__ import annotations

import os
import sys
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

NATOM = 37
MAX_STEPS = 40

HERE = Path.cwd()

def die(msg: str, code: int = 1) -> None:
    print(msg, file=sys.stderr)
    raise SystemExit(code)

def append_line(path: Path, line: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(line.rstrip("\n") + "\n")

def write_line(path: Path, line: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write(line.rstrip("\n") + "\n")

def file_contains(path: Path, needle: str) -> bool:
    if not path.exists():
        return False
    try:
        text = path.read_text(errors="ignore")
    except Exception:
        return False
    return needle in text

def run_py_sub_script(job_name: str, natom: int) -> None:
    # Calls the python version of sub_script (must be in same dir or on PATH)
    cmd = [sys.executable, str(HERE / "sub_script.py"), job_name, str(natom)]
    res = subprocess.run(cmd)
    if res.returncode != 0:
        die(f"sub_script.py failed for {job_name} (exit {res.returncode})", res.returncode)

def concat_files(out_path: Path, parts: list[Path]) -> None:
    with out_path.open("wb") as w:
        for p in parts:
            if not p.exists():
                die(f"Missing required file for concatenation: {p}")
            w.write(p.read_bytes())

def main() -> None:
    # Checks
    if not (HERE / "ProgFile").exists():
        die("ProgFile missing")
    else:
        print("ProgFile Exists - OK", flush=True)

    if not (HERE / "Job0_A.inp").exists():
        die("First Input Missing")
    else:
        print("First Input OK", flush=True)

    #(HERE / "JOBS").mkdir(exist_ok=True)

    report = HERE / "ReportFile"
    if report.exists():
        # mimic: cat ReportFile >> reportfile_old
        old = HERE / "reportfile_old"
        with old.open("ab") as w:
            w.write(report.read_bytes())

    # mimic: date > ReportFile
    write_line(report, datetime.now().ctime())

    num = 0
    while num < MAX_STEPS:
        job_name = f"Job{num}"

        # Run the "sub" step (submits ORCA jobs, extracts ab_initio, runs MECP.py)
        run_py_sub_script(job_name, NATOM)

        # Copy gbw -> start files for next iteration
        a_gbw = HERE / f"{job_name}_A.gbw"
        b_gbw = HERE / f"{job_name}_B.gbw"
        if not a_gbw.exists() or not b_gbw.exists():
            die(f"Missing GBW files after {job_name}: {a_gbw.name} / {b_gbw.name}")

        shutil.copyfile(a_gbw, HERE / "start_A.gbw")
        shutil.copyfile(b_gbw, HERE / "start_B.gbw")

        # Check convergence / error from ReportFile
        if file_contains(report, "CONVERGED"):
            print(f"MECP optimization has converged at Step {num}", flush=True)
            append_line(report, datetime.now().ctime())
            return

        print(f"Step Number {num} -- MECP not yet converged", flush=True)

        if file_contains(report, "ERROR"):
            die("An error has occurred, possibly in the ORCA Job")

        # Prepare next inputs
        num += 1
        concat_files(HERE / f"Job{num}_A.inp", [HERE / "Input_Header_A", HERE / "geom", HERE / "Close_Header"])
        concat_files(HERE / f"Job{num}_B.inp", [HERE / "Input_Header_B", HERE / "geom", HERE / "Close_Header"])

    die(f"Reached MAX_STEPS={MAX_STEPS} without convergence.", 2)

if __name__ == "__main__":
    main()
