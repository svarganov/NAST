
from __future__ import annotations
import numpy as np
from pathlib import Path
from typing import Iterable, Tuple, List

def _rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).ravel()
    if x.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(x * x)))

def _maxabs(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).ravel()
    return float(np.max(np.abs(x))) if x.size else 0.0

def test_convergence(
    N: int,
    Natom: int,
    Nstep: int,
    AtNum: Iterable[str],
    Ea: float,
    Eb: float,
    X_2: Iterable[float],
    X_3: Iterable[float],
    ParG: Iterable[float],
    PerpG: Iterable[float],
    G: Iterable[float],
    # Thresholds
    avg_dx_thresh: float = 2.0e-3,
    max_dx_thresh: float = 4.0e-3,
    avg_grad_thresh: float = 1.0e-4,
    max_grad_thresh: float = 3.0e-4,
    delta_e_thresh: float = 1.0e-4,
    report_path: str = "ReportFile",
    write_geometry_on_every_step: bool = True,
) -> Tuple[int, dict]:
    """
    Computes five criteria (Av.DeltaX, Max.DeltaX, Av.Grad., Max.Grad., DeltaE)

    Returns
    -------
    Conv : int
        1 if all criteria satisfied, else 0.
    metrics : dict
        Dictionary of computed scalar metrics and boolean flags.
    """
    AtNum = list(AtNum)
    X_2 = np.asarray(X_2, dtype=float).reshape(N)
    X_3 = np.asarray(X_3, dtype=float).reshape(N)
    ParG = np.asarray(ParG, dtype=float).reshape(N)
    PerpG = np.asarray(PerpG, dtype=float).reshape(N)
    G = np.asarray(G, dtype=float).reshape(N)

    # Displacement
    dX = X_3 - X_2
    DXRMS = _rms(dX)
    DXMax = _maxabs(dX)

    # Gradients: use the *total* gradient for Av/Max, but also compute parallel/perpendicular RMS for reference
    GRMS = _rms(G)
    GMax = _maxabs(G)
    PGRMS = _rms(PerpG)
    PpGRMS = _rms(ParG)

    # Energy gap between states (MECP criterion): target is to make |Ea - Eb| small
    TDE = abs(Ea - Eb)

    # Pass/fail flags
    PConv = [
        DXRMS <= avg_dx_thresh,
        DXMax <= max_dx_thresh,
        GRMS  <= avg_grad_thresh,
        GMax  <= max_grad_thresh,
        TDE   <= delta_e_thresh,
    ]
    flags = ["YES" if ok else "NO " for ok in PConv]
    Conv = 1 if all(PConv) else 0

    # Write to report
    rp = Path(report_path)
    with rp.open("a", encoding="utf-8") as fh:
        fh.write("\n")
        fh.write(f"Step {Nstep:3d}: Convergence Check\n")
        fh.write("  Criteria (thresholds in parentheses)\n")
        fh.write(f"    Av.DeltaX  = {DXRMS:12.6e}  ({avg_dx_thresh:12.6e})  [{flags[0]}]\n")
        fh.write(f"    Max.DeltaX = {DXMax:12.6e}  ({max_dx_thresh:12.6e})  [{flags[1]}]\n")
        fh.write(f"    Av.Grad    = {GRMS:12.6e}  ({avg_grad_thresh:12.6e})  [{flags[2]}]\n")
        fh.write(f"    Max.Grad   = {GMax:12.6e}  ({max_grad_thresh:12.6e})  [{flags[3]}]\n")
        fh.write(f"    DeltaE     = {TDE:12.6e}  ({delta_e_thresh:12.6e})  [{flags[4]}]\n")
        fh.write("  (Extra) PerpGrad.RMS = %12.6e,  ParGrad.RMS = %12.6e\n" % (PGRMS, PpGRMS))

        if Conv == 1:
            fh.write("  >>> CONVERGED: all criteria satisfied. <<<\n\n")
        else:
            fh.write("  Not converged.\n\n")

        # Geometry printout
        if write_geometry_on_every_step or (Conv == 1):
            try:
                # Expect length(X_3) == 3*Natom
                if len(X_3) != 3 * Natom:
                    raise ValueError("X_3 length is not 3*Natom; cannot format geometry by atoms.")
                fh.write(f"Geometry at Step {Nstep}\n")
                for i in range(Natom):
                    k = 3 * i
                    label = (AtNum[i] if i < len(AtNum) else "X ")[:2]
                    x, y, z = X_3[k:k+3]
                    fh.write(f"{label:>2s}{x:15.7f}{y:15.7f}{z:15.7f}\n")
                fh.write("\n")
            except Exception as e:
                fh.write(f"[Note] Geometry print skipped: {e}\n\n")

    metrics = dict(
        DXRMS=DXRMS, DXMax=DXMax, GRMS=GRMS, GMax=GMax, TDE=TDE,
        PGRMS=PGRMS, PpGRMS=PpGRMS, flags=flags, Conv=Conv
    )
    return Conv, metrics

if __name__ == "__main__":
    # Minimal smoke test
    Natom = 2
    N = 3 * Natom
    AtNum = ["C ", "H "]
    Ea, Eb = -100.0, -100.00002
    X_2 = np.zeros(N)
    X_3 = np.array([1e-4, -2e-4, 1e-4, -1e-4, 1e-4, -1e-4])
    G   = np.array([2e-4, -1e-4, 1e-4,  1e-4, -2e-4,  3e-4])
    ParG = 0.5 * G
    PerpG = 0.5 * G
    Conv, metrics = test_convergence(N, Natom, 1, AtNum, Ea, Eb, X_2, X_3, ParG, PerpG, G)
    print("Conv:", Conv)
    for k, v in metrics.items():
        if k != "flags":
            print(f"{k}:", v)
    print("flags:", metrics["flags"])
    print("Report written to 'ReportFile'")
