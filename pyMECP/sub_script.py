#!/usr/bin/env python3
"""
Usage:
  python sub_script.py <JobName> <natom>

What it does (same as csh):
- submits A and B ORCA jobs via sbatch (orca.sbatch)
- waits for completion
- checks for "****ORCA TERMINATED NORMALLY****" in outputs
- builds ab_initio file (extract_energy, extract_symbol, extract_gradient)
- runs MECP.py
- appends AddtoReportFile -> ReportFile and removes AddtoReportFile

Notes:
- This script assumes:
  * orca.sbatch is in the current working directory
  * ORCA output files are <JobName>_A.out and <JobName>_B.out
"""
from __future__ import annotations

import os
import re
import sys
import time
import subprocess
from pathlib import Path
from time import sleep

HERE = Path.cwd()

ORCA_OK_LINE = "****ORCA TERMINATED NORMALLY****"

def die(msg: str, code: int = 1) -> None:
    print(msg, file=sys.stderr)
    raise SystemExit(code)

def run(cmd: str, *, check: bool = False, capture: bool = False) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        shell=True,
        check=check,
        text=True,
        capture_output=capture,
    )

def submit_job(inp: str, gbw: str) -> str:
    """
    Submit sbatch job and return jobid (string). Raises if submission fails or jobid can't be parsed.
    """
    cp = run(f"sbatch orca.sbatch {inp} {gbw}", capture=True)
    if cp.returncode != 0:
        die(f"sbatch failed: {cp.stderr.strip()}")
    m = re.search(r"Submitted batch job\s+(\d+)", cp.stdout)
    if not m:
        die(f"Could not parse jobid from sbatch output: {cp.stdout.strip()}")
    return m.group(1)

def wait_for_job(jobid: str, poll_s: int = 10) -> None:
    """
    Wait until squeue no longer shows the job.
    """
    while True:
        cp = run(f"squeue -j {jobid} -h", capture=True)
        # If job is not in queue, stdout should be empty.
        if cp.returncode == 0 and cp.stdout.strip() == "":
            return
        time.sleep(poll_s)

def wait_for_file(path: Path, timeout_s: int = 900, poll_s: int = 30) -> None:
    t0 = time.time()
    while time.time() - t0 < timeout_s:
        if path.exists() and path.stat().st_size > 0:
            return
        time.sleep(poll_s)
    die(f"Timed out waiting for file: {path}")

def wait_for_orca_ok(out_path: Path, timeout_s: int = 7200, poll_s: int = 10) -> None:
    """Wait until the ORCA normal termination line appears in the output file."""
    t0 = time.time()
    while time.time() - t0 < timeout_s:
        if out_path.exists():
            try:
                txt = out_path.read_text(errors="ignore")
            except Exception:
                txt = ""
            if ORCA_OK_LINE in txt:
                return
        time.sleep(poll_s)
    # helpful debug before dying
    print(f"[ERROR] Timed out waiting for ORCA normal termination in {out_path}", file=sys.stderr)
    if out_path.exists():
        try:
            tail = run(f"tail -n 80 {out_path.name}", capture=True)
            print(tail.stdout, file=sys.stderr)
        except Exception:
            pass
    die(f"ORCA did not report normal termination in time: {out_path}")


def file_contains(path: Path, needle: str) -> bool:
    if not path.exists():
        return False
    try:
        txt = path.read_text(errors="ignore")
    except Exception:
        return False
    return needle in txt

def append_report_error() -> None:
    with (HERE / "ReportFile").open("a", encoding="utf-8") as f:
        f.write("ERROR\n")

# --- Extraction Logic ---

def get_energy(engrad_path: Path) -> str:
    """Reads the energy from the .engrad file (usually the line after # THE CURRENT TOTAL ENERGY)."""
    lines = engrad_path.read_text().splitlines()
    for i, line in enumerate(lines):
        if "# The current total energy in Eh" in line:
            return lines[i+2].strip()
    die(f"Could not find energy in {engrad_path}")

def get_gradients(engrad_path: Path, natom: int) -> list[list[str]]:
    """Reads gradients from .engrad: 3 lines (X, Y, Z) per atom."""
    lines = engrad_path.read_text().splitlines()
    start_idx = -1
    for i, line in enumerate(lines):
        if "# The current gradient in Eh/bohr" in line:
            start_idx = i + 2 # Skip header and the 'number of gradients' line
            break
    
    if start_idx == -1:
        die(f"Could not find gradients in {engrad_path}")
    
    grads = []
    for a in range(natom):
        # Extract 3 consecutive lines for X, Y, Z
        idx = start_idx + (a * 3)
        grads.append([lines[idx].strip(), lines[idx+1].strip(), lines[idx+2].strip()])
    return grads

def get_symbols(out_path: Path, natom: int) -> list[str]:
    """Extracts atom symbols from the .out file."""
    lines = out_path.read_text(errors="ignore").splitlines()
    symbols = []
    for i, line in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
            # Skip the dashed line header
            for j in range(1, natom + 1):
                parts = lines[i + 1 + j].split()
                if parts:
                    symbols.append(parts[0])
            return symbols
    die(f"Could not find Cartesian Coordinates in {out_path}")

def build_ab_initio(name: str, natom: int) -> None:
    """Consolidates data into the 'ab_initio' file."""
    ab_initio_path = HERE / "ab_initio"
    
    with ab_initio_path.open("w", encoding="utf-8") as f:
        for suffix, label in [("A", "First"), ("B", "Second")]:
            out_file = HERE / f"{name}_{suffix}.out"
            engrad_file = HERE / f"{name}_{suffix}.engrad"
            
            if not engrad_file.exists():
                die(f"Missing engrad file: {engrad_file}")

            # 1. Energy
            f.write(f"Energy of the {label} State\n")
            f.write(f"{get_energy(engrad_file)}\n")
            
            # 2. Gradients with Symbols
            f.write(f"Gradient of the {label} State\n")
            symbols = get_symbols(out_file, natom)
            grads = get_gradients(engrad_file, natom)
            
            for sym, g_xyz in zip(symbols, grads):
                # Formats: Symbol  GradX  GradY  GradZ
                f.write(f"{sym}   {'   '.join(g_xyz)}\n")

def main() -> None:
    if len(sys.argv) != 3:
        die("Usage: python sub_script.py <JobName> <natom>")

    name = sys.argv[1]
    natom = int(sys.argv[2])

    a_out = HERE / f"{name}_A.out"
    b_out = HERE / f"{name}_B.out"

    jobids: list[str] = []

    # Mimic original conditional submission for Job0
    if name == "Job0":
        if not a_out.exists():
            jobids.append(submit_job(f"{name}_A.inp", "start_A.gbw"))
        if not b_out.exists():
            jobids.append(submit_job(f"{name}_B.inp", "start_B.gbw"))
    else:
        jobids.append(submit_job(f"{name}_A.inp", "start_A.gbw"))
        jobids.append(submit_job(f"{name}_B.inp", "start_B.gbw"))

    # Wait for submitted jobs
    #for jid in jobids:
        #wait_for_job(jid)

    # Safety: ensure output files exist (in case job finished but FS latency)
    #wait_for_orca_ok(a_out, timeout_s=int(os.environ.get("ORCA_OK_TIMEOUT_S", "7200")))
    #wait_for_orca_ok(b_out, timeout_s=int(os.environ.get("ORCA_OK_TIMEOUT_S", "7200")))

    sleep(1200)

    # Optional: print the termination lines like the csh script did via grep
    if a_out.exists():
        run(f"grep \"{ORCA_OK_LINE}\" {a_out.name}")
    if b_out.exists():
        run(f"grep \"{ORCA_OK_LINE}\" {b_out.name}")

    # Validate success
    if not file_contains(a_out, ORCA_OK_LINE):
        print("There has been a problem in ORCA", file=sys.stderr)
        append_report_error()
        raise SystemExit(1)

    if not file_contains(b_out, ORCA_OK_LINE):
        print("There has been a problem in ORCA", file=sys.stderr)
        append_report_error()
        raise SystemExit(1)

    # Build ab_initio using your existing awk scripts
    build_ab_initio(name, natom)

    # Run MECP.py with the same python as this script
    mecp = subprocess.run([sys.executable, "MECP.py"])
    if mecp.returncode != 0:
        print("Problem with pyMECP program", file=sys.stderr)
        append_report_error()
        raise SystemExit(mecp.returncode)

    # Append AddtoReportFile into ReportFile, then delete it
    add = HERE / "AddtoReportFile"
    report = HERE / "ReportFile"
    if add.exists():
        with report.open("ab") as w:
            w.write(add.read_bytes())
        add.unlink(missing_ok=True)

if __name__ == "__main__":
    main()
