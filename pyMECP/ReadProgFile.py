
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple
import numpy as np
import re


_FLOAT_RE = re.compile(
    r"""[-+]? (?:\d+\.\d*|\.\d+|\d+) (?:[eEdD][-+]?\d+)?""",
    re.VERBOSE,
)

def _floats_in(s: str) -> List[float]:
    return [float(x.replace("D", "E").replace("d", "e")) for x in _FLOAT_RE.findall(s)]

def _first_int_in(s: str) -> Optional[int]:
    m = re.search(r"[-+]?\d+", s)
    return int(m.group(0)) if m else None

def _find_line(lines: List[str], key: str) -> Optional[int]:
    key = key.lower()
    for i, ln in enumerate(lines):
        if key in ln.lower():
            return i
    return None

def _read_int_same_or_next(lines: List[str], idx: int, max_look=3) -> Optional[int]:
    """Read an integer on the same line or within the next few lines."""
    for j in range(idx, min(idx + 1 + max_look, len(lines))):
        v = _first_int_in(lines[j])
        if v is not None:
            return v
    return None

def _collect_floats_after(lines: List[str], start_idx: int, need: int) -> Tuple[np.ndarray, int]:
    """Collect `need` floats starting from the line AFTER start_idx."""
    out: List[float] = []
    i = start_idx + 1
    while i < len(lines) and len(out) < need:
        out.extend(_floats_in(lines[i]))
        i += 1
    if len(out) < need:
        raise ValueError(f"Expected {need} numbers after line {start_idx}, got {len(out)}.")
    return np.asarray(out[:need], float), i

def _collect_geom_labeled(lines: List[str], start_idx: int, natom: int) -> Tuple[List[str], np.ndarray, int]:
    """Read natom lines of: label + 3 floats (e.g., '(A2,3F...)')."""
    labels: List[str] = []
    coords: List[float] = []
    i = start_idx + 1
    for _ in range(natom):
        if i >= len(lines):
            raise ValueError("Unexpected EOF while reading labeled geometry.")
        ln = lines[i].rstrip("\n")
        lab = (ln[:2].strip() or "X")[:2]
        vals = _floats_in(ln[2:]) or _floats_in(ln)  # be forgiving
        if len(vals) < 3:
            raise ValueError(f"Not enough numbers on geometry line {i}: '{ln}'")
        labels.append(lab)
        coords.extend(vals[:3])
        i += 1
    return labels, np.asarray(coords, float), i

# --- data structure -----------------------------------------------------------

@dataclass
class ProgState:
    natom: int
    nx: int
    nstep: Optional[int] = None
    full: Optional[bool] = None
    AtNum: Optional[List[str]] = None
    X_next: Optional[np.ndarray] = None
    X_prev: Optional[np.ndarray] = None
    Ea: Optional[float] = None
    Eb: Optional[float] = None
    Ga: Optional[np.ndarray] = None
    Gb: Optional[np.ndarray] = None
    Geff: Optional[np.ndarray] = None
    HI: Optional[np.ndarray] = None

# --- main reader --------------------------------------------------------------

def read_prog_file(path: str = "ProgFile") -> ProgState:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Cannot open '{path}'")

    lines = p.read_text(errors="ignore").splitlines()

    # Number of Atoms (allow value on next line)
    i_nat = _find_line(lines, "Number of Atoms")
    if i_nat is None:
        raise ValueError("Missing 'Number of Atoms' in ProgFile.")
    natom = _read_int_same_or_next(lines, i_nat)
    if not natom or natom <= 0:
        raise ValueError("Failed to parse a positive integer for Number of Atoms.")
    nx = 3 * natom
    st = ProgState(natom=natom, nx=nx)

    # Steps already run (value may be on next line)
    i_step = _find_line(lines, "Number of Steps already Run")
    if i_step is not None:
        st.nstep = _read_int_same_or_next(lines, i_step)

        # Next Geometry to Compute (prefer labeled format)
    i_next = _find_line(lines, "Next Geometry to Compute")
    # Accept older/newer variants:
    if i_next is None:
        i_next = _find_line(lines, "Current Geometry")
    if i_next is None:
        i_next = _find_line(lines, "Geometry:")  # single generic block
    
    if i_next is not None:
        try:
            labs, xyz, _ = _collect_geom_labeled(lines, i_next, natom)
            st.AtNum = labs
            st.X_next = xyz
        except Exception:
            arr, _ = _collect_floats_after(lines, i_next, nx)
            st.X_next = arr

    # If this file only had one 'Geometry' block and there is no explicit
    # 'Previous Geometry', seed previous as zeros for step 0.
    if st.X_prev is None and (st.nstep in (None, 0)):
        st.X_prev = np.zeros(nx, float)


    # Next Geometry to Compute (prefer labeled format)
    i_next = _find_line(lines, "Next Geometry to Compute")
    if i_next is not None:
        try:
            labs, xyz, _ = _collect_geom_labeled(lines, i_next, natom)
            st.AtNum = labs
            st.X_next = xyz
        except Exception:
            arr, _ = _collect_floats_after(lines, i_next, nx)
            st.X_next = arr

    # Previous Geometry (numbers only in your file)
    i_prev = _find_line(lines, "Previous Geometry")
    if i_prev is not None:
        try:
            arr, _ = _collect_floats_after(lines, i_prev, nx)
            st.X_prev = arr
        except Exception:
            # fallback: maybe labeled
            _, xyz, _ = _collect_geom_labeled(lines, i_prev, natom)
            st.X_prev = xyz

    # Energies
    i_E = _find_line(lines, "Energies of First, Second State")
    if i_E is not None:
        vals, _ = _collect_floats_after(lines, i_E, 2)
        st.Ea, st.Eb = float(vals[0]), float(vals[1])

    # Gradients
    i_Ga = _find_line(lines, "Gradient of First State")
    if i_Ga is None:
        i_Ga = _find_line(lines, "Gradient of the First State")
    if i_Ga is not None:
        st.Ga, _ = _collect_floats_after(lines, i_Ga, nx)

    i_Gb = _find_line(lines, "Gradient of Second State")
    if i_Gb is None:
        i_Gb = _find_line(lines, "Gradient of the Second State")
    if i_Gb is not None:
        st.Gb, _ = _collect_floats_after(lines, i_Gb, nx)

    # Effective Gradient
    i_Geff = _find_line(lines, "Effective Gradient")
    if i_Geff is not None:
        st.Geff, _ = _collect_floats_after(lines, i_Geff, nx)

    # Approximate Inverse Hessian
    i_HI = _find_line(lines, "Approximate Inverse Hessian")
    if i_HI is not None:
        flat, _ = _collect_floats_after(lines, i_HI, nx * nx)
        st.HI = flat.reshape(nx, nx)

    return st

if __name__ == "__main__":
    try:
        st = read_prog_file("")
    except Exception as e:
        print(f"[ReadProgFile] Error: {e}")
    else:
        print(f"natom={st.natom}, nx={st.nx}, nstep={st.nstep}, full={st.full}")
        if st.X_next is not None:
            print("X_next len:", st.X_next.size, "max|.|:", float(np.max(np.abs(st.X_next))))
        if st.X_prev is not None:
            print("X_prev len:", st.X_prev.size, "max|.|:", float(np.max(np.abs(st.X_prev))))
        if (st.Ea is not None) and (st.Eb is not None):
            print(f"Ea={st.Ea:.8f}  Eb={st.Eb:.8f}  dE={st.Ea-st.Eb:.6e}")
        if st.Geff is not None:
            print("Geff len:", st.Geff.size, "max|.|:", float(np.max(np.abs(st.Geff))))
        if st.HI is not None:
            d = np.diag(st.HI)
            print("HI shape:", st.HI.shape, "diag min/max:", float(np.min(d)), float(np.max(d)))

