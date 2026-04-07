"""
    SUBROUTINE WriteGeomFile(Natom, Nx, AtNum, X)

Behavior:
- Writes a file named "geom" containing Natom lines.
- Each line has a 2-character atom label followed by three coordinates from X.
- Coordinates come from X in blocks of 3 (x, y, z) for each atom, i.e. X[3*i : 3*i+3].
- Adds a blank line at the end of the file.
"""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Sequence, Union
import numpy as np

Number = Union[int, float, np.number]

def write_geom_file(
    Natom: int,
    Nx: int,
    AtNum: Sequence[str],
    X: Iterable[Number],
    path: str = "",
    strict_lengths: bool = True,
) -> None:
    """
    Parameters
    ----------
    Natom : int
        Number of atoms.
    Nx : int
        Number of coordinates (expected 3*Natom).
    AtNum : Sequence[str]
        Atom labels; only first 2 characters of each are written.
    X : Iterable[Number]
        Flat iterable of length Nx containing coordinates in Angstroms (x1,y1,z1,x2,y2,z2,...).
    path : str
        Output filename (default: "geom").
    strict_lengths : bool
        If True (default), require len(X) == Nx == 3*Natom. If False, will write min(Natom, len(X)//3) atoms.
    """
    X = np.asarray(list(X), dtype=float).ravel()
    if strict_lengths:
        if Nx != 3 * Natom:
            raise ValueError(f"Nx must be 3*Natom (got Nx={Nx}, Natom={Natom}).")
        if X.size != Nx:
            raise ValueError(f"X length mismatch: expected {Nx}, got {X.size}.")
        count_atoms = Natom
    else:
        # be forgiving
        count_atoms = min(Natom, X.size // 3)

    p = Path(path)
    with p.open("w", encoding="utf-8") as fh:
        for i in range(count_atoms):
            k = 3 * i
            label = (AtNum[i] if i < len(AtNum) else "X ")[:2]
            x, y, z = X[k:k+3]
            fh.write(f"{label:>2s}{x:14.8f}{y:14.8f}{z:14.8f}\n")
        fh.write("\n")

if __name__ == "__main__":
    # Minimal smoke test
    Natom = 2
    Nx = 3 * Natom
    AtNum = ["C ", "H "]
    X = [0.0, 1.23456789, -2.0,  1.0, -0.5, 3.14159265]
    write_geom_file(Natom, Nx, AtNum, X, path="")
    print("Wrote 'geom' with", Natom, "atoms.")
