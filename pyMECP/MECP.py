# MECP.py
"""
Workflow (one optimization step per call):
  1) Read previous state from ProgFile (geometry, inverse Hessian, etc.)
  2) Read Ea/Eb and Ga/Gb from 'ab_initio' (created from ORCA outputs)
  3) Build effective MECP gradient
  4) Update geometry (BFGS-like; step capped inside updatex.py)
  5) Check convergence
  6) Write: 'geom' (next coords), 'ProgFile' (updated state),
            'AddtoReportFile' (short text that master_script greps)
"""

from __future__ import annotations
from pathlib import Path
import numpy as np

from ReadProgFile import read_prog_file, ProgState
from ReadInput import read_input
from Effective_gradient import effective_gradient
from updatex import update_x
from TestConvergence import test_convergence
from WriteGeomFile import write_geom_file
from WriteProgFile import write_prog_file

HERE = Path.cwd()

def main():
    # --- 1) Restore previous optimizer state from ProgFile --------------------
    if not Path("ProgFile").exists():
        raise RuntimeError(
            "ProgFile not found. Run your initialization step first "
            "(create initial 'geom' and 'ProgFile')."
        )

    st: ProgState = read_prog_file(HERE/"ProgFile")
    Natom = st.natom
    Nx = 3 * Natom

    # tolerate partial ProgFile content
    nstep = int(st.nstep or 0)
    AtNum = st.AtNum or ["X"] * Natom
    X1 = st.X_prev if st.X_prev is not None else np.zeros(Nx)
    X2 = st.X_next if st.X_next is not None else np.zeros(Nx)
    HI1 = st.HI if st.HI is not None else np.eye(Nx)
    G1 = st.Geff if st.Geff is not None else np.zeros(Nx)
    FFile = int(bool(st.full)) if (st.full is not None) else 0

    # --- 2) Read energies/gradients aggregated into 'ab_initio' --------------
    # (Your sub_script creates 'ab_initio' from ORCA outputs via awk.)
    Ea, Ga, Eb, Gb = read_input(natom=Natom, nx=Nx, path=HERE/"ab_initio")

    # --- 3) Build effective MECP gradient ------------------------------------
    G2, ParG, PerpG = effective_gradient(Ea, Eb, Ga, Gb)

    # --- 4) Take an optimization step (BFGS-like with capping) ----------------
    X3, HI2, step = update_x(
        N=Nx, Nstep=nstep + 1, FFile=FFile,
        X_1=X1, X_2=X2, G_1=G1, G_2=G2, HI1=HI1
    )

    # --- 5) Convergence check --------------------------------------------------
    Conv = test_convergence(
        Nx, Natom, nstep + 1, np.array(AtNum),
        Ea, Eb, X2, X3, ParG, PerpG, G2
    )

    # --- 6) Write outputs expected by the shell workflow ----------------------
    # 6a) New geometry for the next ORCA job
    write_geom_file(Natom, Nx, AtNum, X3, path=HERE/"geom")

    # 6b) Updated progress/state
    write_prog_file(
        Natom=Natom, Nx=Nx, AtNum=AtNum, Nstep=nstep + 1,
        X_prev=X2, X_next=X3, HI=HI2,
        Ea=Ea, Eb=Eb, Ga=Ga, Gb=Gb, G_eff=G2,
        path=HERE/"ProgFile", full_flag=1
    )

    # 6c) Short summary; master_script greps for "CONVERGED"
    with open("AddtoReportFile", "w", encoding="utf-8") as f:
        f.write(f"Step {nstep+1}\n")
        f.write(f"Energy of the First State:  {Ea: .8f}\n")
        f.write(f"Energy of the Second State: {Eb: .8f}\n")
        f.write(f"Difference in E:            {Ea - Eb: .8f}\n")
        f.write(f"Max |Grad|:                 {np.max(np.abs(G2)): .6e}\n")
        f.write(f"Max |Step|:                 {np.max(np.abs(step)): .6e}\n")
        if Conv == 0:
            f.write("The MECP Optimization has CONVERGED at that geometry !!!\n")
            f.write("CONVERGED\n")

#    print(
#        f"[MECP.py] step {nstep}: dE={Ea-Eb:.6e}, "
#        f"max|G|={np.max(np.abs(G2)):.3e}, max|step|={np.max(np.abs(step)):.3e}, "
#        f"{'CONVERGED' if Conv else 'continuing'}"
#    )

    print(
        f"[MECP.py] step {nstep}: dE={Ea-Eb:.6e}, "
        f"max|G|={np.max(np.abs(G2)):.3e}, max|step|={np.max(np.abs(step)):.3e} ")

if __name__ == "__main__":
    main()
