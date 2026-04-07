"""
Initializes:
- Ea, Eb = 0.0
- Ga, Gb = zero vectors (length n)
- H = n×n matrix with diagonal = 0.7, off-diagonals = 0.0
"""
from __future__ import annotations
from typing import Tuple
import numpy as np

def initialize(n: int) -> Tuple[np.ndarray, float, float, np.ndarray, np.ndarray]:
    H = np.zeros((n, n), dtype=float)
    Ea = 0.0
    Eb = 0.0
    Ga = np.zeros(n, dtype=float)
    Gb = np.zeros(n, dtype=float)
    for i in range(n):
        H[i, i] = 0.7
        for j in range(i + 1, n):
            H[i, j] = 0.0
            H[j, i] = 0.0
    return H, Ea, Eb, Ga, Gb

if __name__ == "__main__":
    n = 4
    H, Ea, Eb, Ga, Gb = initialize(n)
    print("H:\n", H)
    print("Ea:", Ea, "Eb:", Eb)
    print("Ga:", Ga, "Gb:", Gb)
