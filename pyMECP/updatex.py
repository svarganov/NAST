import numpy as np

def update_x(N, Nstep, FFile, X_1, X_2, G_1, G_2, HI1):
    """
    Parameters
    ----------
    N : int
        Dimension of the coordinate/gradient vectors.
    Nstep : int
        Current step number.
    FFile : int
        Flag (Nstep==0 and FFile==0 triggers identity HI_2).
    X_1, X_2 : (N,) array_like
        Consecutive coordinates (previous, current).
    G_1, G_2 : (N,) array_like
        Gradients at X_1 and X_2 respectively.
    HI1 : (N,N) array_like
        Current inverse Hessian approximation H^{-1}.

    Returns
    -------
    X_3 : (N,) ndarray
        Next coordinates after applying the BFGS step with componentwise capping.
    HI_2 : (N,N) ndarray
        Updated inverse Hessian approximation.
    ChgeX : (N,) ndarray
        The raw step that was applied (may be scaled if too large).
    """
    X_1 = np.asarray(X_1, dtype=float).reshape(N)
    X_2 = np.asarray(X_2, dtype=float).reshape(N)
    G_1 = np.asarray(G_1, dtype=float).reshape(N)
    G_2 = np.asarray(G_2, dtype=float).reshape(N)
    HI1 = np.asarray(HI1, dtype=float).reshape(N, N)

    STPMX = 0.1  # per-component max step

    # Initialize HI_2 ~ Identity on the very first call (mirrors IF (Nstep==0 .and. FFile==0))
    if (Nstep == 0) and (FFile == 0):
        HI2 = np.eye(N, dtype=float)
    else:
        # BFGS inverse-Hessian update
        s = (X_2 - X_1)          # Δx
        y = (G_2 - G_1)          # Δg

        # Safeguards
        ys = float(np.dot(y, s))
        if ys <= 1e-14:
            # Skip update if not positive-definite
            HI2 = HI1.copy()
        else:
            rho = 1.0 / ys
            Hy = HI1 @ y
            yHy = float(y @ Hy)
            # BFGS formula for inverse Hessian:
            term1 = (1.0 + yHy * rho) * (rho * np.outer(s, s))
            term2 = rho * (np.outer(s, Hy) + np.outer(Hy, s))
            HI2 = HI1 + term1 - term2

    # Compute the quasi-Newton step: ChgeX = - H^{-1} * g_current
    ChgeX = - (HI2 @ G_2)

    # Componentwise cap (find largest |ChgeX(i)|; if > STPMX, scale all components)
    lgstst = float(np.max(np.abs(ChgeX))) if ChgeX.size else 0.0
    if lgstst > STPMX:
        ChgeX = (ChgeX / lgstst) * STPMX

    # Advance
    X_3 = X_2 + ChgeX

    return X_3, HI2, ChgeX

if __name__ == "__main__":
    # Tiny self-test with random numbers
    N = 5
    X_1 = np.zeros(N)
    X_2 = np.random.randn(N) * 0.01
    G_1 = np.random.randn(N)
    G_2 = np.random.randn(N)
    HI1 = np.eye(N)
    X_3, HI2, step = update_x(N, Nstep=1, FFile=1, X_1=X_1, X_2=X_2, G_1=G_1, G_2=G_2, HI1=HI1)
    print("X_3:", X_3)
    print("||step||_inf:", np.max(np.abs(step)))
