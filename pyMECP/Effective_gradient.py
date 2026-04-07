from __future__ import annotations
from typing import Tuple
import numpy as np

def effective_gradient(Ea: float,
                       Eb: float,
                       Ga: np.ndarray,
                       Gb: np.ndarray,
                       facPP: float = 140.0,
                       facP: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    Ga = np.asarray(Ga, dtype=float)
    Gb = np.asarray(Gb, dtype=float)
    if Ga.shape != Gb.shape:
        raise ValueError(f"Ga and Gb must have the same shape; got {Ga.shape} vs {Gb.shape}")
    PerpG = Ga - Gb
    norm_perp = np.linalg.norm(PerpG)
    if norm_perp == 0.0:
        ParG = Ga.copy()
        G = facP * ParG
        return G, ParG, PerpG
    pp = np.dot(Ga, PerpG) / norm_perp
    ParG = Ga - (PerpG / norm_perp) * pp
    G = (Ea - Eb) * facPP * PerpG + facP * ParG
    return G, ParG, PerpG

if __name__ == "__main__":
    Ea, Eb = -100.123, -100.100
    Ga = np.array([0.01, 0.02, -0.03])
    Gb = np.array([0.005, -0.01, -0.02])
    G, ParG, PerpG = effective_gradient(Ea, Eb, Ga, Gb)
    print("PerpG:", PerpG)
    print("ParG:", ParG)
    print("Effective G:", G)

