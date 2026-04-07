
from __future__ import annotations
from pathlib import Path
from typing import Tuple, List
import numpy as np
import re
from Effective_gradient import effective_gradient

BOHR_ANG = 0.529177  # Å per Bohr

_FLOAT_RE = re.compile(
    r"""[-+]? (?:\d+\.\d*|\.\d+|\d+) (?:[eEdD][-+]?\d+)?""",
    re.VERBOSE,
)

def _floats_in(s: str) -> List[float]:
    out = []
    for t in _FLOAT_RE.findall(s):
        out.append(float(t.replace("D","E").replace("d","e")))
    return out

def read_input(natom: int, nx: int, path: str = "ab_initio",
               assume_ga_is_force: bool = False,
               assume_gb_is_force: bool = False) -> Tuple[float, np.ndarray, float, np.ndarray]:
    """
    Parse an 'ab_initio' file with sections:
      Energy of the First State
      <Ea>
      Gradient of the First State
      <natom lines: label + 3 numbers>
      Energy of the Second State
      <Eb>
      Gradient of the Second State
      <natom lines: label + 3 numbers>
    Returns gradients in Hartree/Å.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Cannot open '{path}'")
    lines = p.read_text(errors="ignore").splitlines()

    # Helper to find a line index containing a key (case-insensitive)
    def find(key: str, start: int = 0) -> int:
        key_l = key.lower()
        for i in range(start, len(lines)):
            if key_l in lines[i].lower():
                return i
        raise ValueError(f"Missing section header: '{key}'")

    # 1) Ea
    i_energy_a = find("Energy of the First State", 0)
    # first numeric line after the header
    Ea = None
    i = i_energy_a + 1
    while i < len(lines) and Ea is None:
        nums = _floats_in(lines[i])
        if nums:
            Ea = float(nums[0])
        i += 1
    if Ea is None:
        raise ValueError("Failed to read Ea")

    # 2) Ga
    i_grad_a = find("Gradient of the First State", i_energy_a)
    Ga_vals: List[float] = []
    j = i_grad_a + 1
    for k in range(natom):
        if j >= len(lines):
            raise ValueError(f"Unexpected EOF while reading Ga at atom {k+1}/{natom}")
        nums = _floats_in(lines[j])
        if len(nums) < 3:
            raise ValueError(f"Not enough numbers on Ga line {j+1}: '{lines[j]}'")
        Ga_vals.extend(nums[:3])
        j += 1
    if len(Ga_vals) != nx:
        raise ValueError(f"Ga length mismatch: expected {nx}, got {len(Ga_vals)}")
    Ga = np.asarray(Ga_vals, float)

    # 3) Eb
    i_energy_b = find("Energy of the Second State", i_grad_a)
    Eb = None
    i = i_energy_b + 1
    while i < len(lines) and Eb is None:
        nums = _floats_in(lines[i])
        if nums:
            Eb = float(nums[0])
        i += 1
    if Eb is None:
        raise ValueError("Failed to read Eb")

    # 4) Gb
    i_grad_b = find("Gradient of the Second State", i_energy_b)
    Gb_vals: List[float] = []
    j = i_grad_b + 1
    for k in range(natom):
        if j >= len(lines):
            raise ValueError(f"Unexpected EOF while reading Gb at atom {k+1}/{natom}")
        nums = _floats_in(lines[j])
        if len(nums) < 3:
            raise ValueError(f"Not enough numbers on Gb line {j+1}: '{lines[j]}'")
        Gb_vals.extend(nums[:3])
        j += 1
    if len(Gb_vals) != nx:
        raise ValueError(f"Gb length mismatch: expected {nx}, got {len(Gb_vals)}")
    Gb = np.asarray(Gb_vals, float)

    # Unit conversion: Hartree/Bohr -> Hartree/Å
    Ga = Ga / BOHR_ANG
    Gb = Gb / BOHR_ANG

    # If the file provided forces, flip sign to get gradients
    if assume_ga_is_force:
        Ga = -Ga
    if assume_gb_is_force:
        Gb = -Gb

    return Ea, Ga, Eb, Gb

if __name__ == "__main__":
    # quick smoke: adjust path & natom as needed
    natom = 5
    nx = 3 * natom
    path = ""
    Ea, Ga, Eb, Gb = read_input(natom, nx, path)
    print("Ea:", Ea)
    print("Eb:", Eb)
    print("Ga shape:", Ga.shape, "Gb shape:", Gb.shape)
    print("Ga[:6]:", Ga[:6])
    print("Gb[:6]:", Gb[:6])
    G2, ParG, PerpG = effective_gradient(Ea, Eb, Ga, Gb)
    print("G2:", G2)
    print("ParG:", ParG)
    print("PeroG:",PerpG)
