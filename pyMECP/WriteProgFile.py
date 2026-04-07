from pathlib import Path
from typing import Sequence, Iterable, Union
import numpy as np

Number = Union[int, float, np.number]

def _as_vec(x: Iterable[Number], n: int, name: str) -> np.ndarray:
    arr = np.asarray(list(x), dtype=float).ravel()
    if arr.size != n:
        raise ValueError(f"{name} must have length {n}, got {arr.size}")
    return arr

def _as_mat(m: Iterable[Iterable[Number]], n: int, name: str) -> np.ndarray:
    M = np.asarray(m, dtype=float)
    if M.shape != (n, n):
        raise ValueError(f"{name} must have shape ({n},{n}), got {M.shape}")
    return M

def write_prog_file(
    Natom: int,
    Nx: int,
    AtNum: Sequence[str],
    Nstep: int,
    X_prev: Iterable[Number],       
    X_next: Iterable[Number],   
    HI: Iterable[Iterable[Number]],
    Ea: Number,
    Eb: Number,
    Ga: Iterable[Number],
    Gb: Iterable[Number],
    G_eff: Iterable[Number],
    path: str = "ProgFile",
    full_flag: int = 1,
) -> None:
    """Write a ProgFile."""
    if Nx != 3 * Natom:
        raise ValueError(f"Nx must equal 3*Natom (got Nx={Nx}, Natom={Natom})")

    # Normalize arrays
    X_prev = _as_vec(X_prev, Nx, "X_prev")
    X_next = _as_vec(X_next, Nx, "X_next")
    Ga     = _as_vec(Ga,     Nx, "Ga")
    Gb     = _as_vec(Gb,     Nx, "Gb")
    G_eff  = _as_vec(G_eff,  Nx, "G_eff")
    HI     = _as_mat(HI,     Nx, "HI")

    F20_12 = "{v:20.12f}"
    def f20(v): return F20_12.format(v=float(v))

    p = Path(path)
    with p.open("w", encoding="utf-8") as fh:
        # Headers & counts
        fh.write("Progress File for MECP Optimization\n")
        fh.write("Number of Atoms:\n")
        fh.write(f"{Natom}\n")
        fh.write("Number of Steps already Run\n")
        fh.write(f"{Nstep}\n")
        fh.write("Is this a full ProgFile ?\n")
        fh.write(f"{int(full_flag)}\n")

        # Next Geometry to Compute
        fh.write("Next Geometry to Compute:\n")
        for i in range(Natom):
            k = 3 * i
            lab = (AtNum[i] if i < len(AtNum) else "X ")[:2]
            x, y, z = X_next[k:k+3]
            fh.write(f"{lab:>2s}{f20(x)}{f20(y)}{f20(z)}\n")

        # Previous Geometry
        fh.write("Previous Geometry:\n")
        for i in range(Natom):
            k = 3 * i
            x, y, z = X_prev[k:k+3]
            fh.write(f"{f20(x)}{f20(y)}{f20(z)}\n")

        # Energies
        fh.write("Energies of First, Second State at that Geometry:\n")
        fh.write(f"{f20(Ea)}\n")
        fh.write(f"{f20(Eb)}\n")

        # Gradients
        fh.write("Gradient of First State at that Geometry:\n")
        for v in Ga: fh.write(f"{f20(v)}\n")
        fh.write("Gradient of Second State at that Geometry:\n")
        for v in Gb: fh.write(f"{f20(v)}\n")

        # Effective Gradient
        fh.write("Effective Gradient at that Geometry:\n")
        for v in G_eff: fh.write(f"{f20(v)}\n")

        # Approximate Inverse Hessian
        fh.write("Approximate Inverse Hessian at that Geometry:\n")
        for i in range(Nx):
            for j in range(Nx):
                fh.write(f"{f20(HI[i, j])}\n")
