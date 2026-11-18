"""模型のハミルトニアンと観測量（行列表示）。"""
import numpy as np

from .operators import X, Z, site_op


def tfim_hamiltonian(n, J, h, pbc=False):
    """1次元横磁場イジングモデル H = -J Σ Z_i Z_{i+1} - h Σ X_i。

    pbc=True: 周期境界条件（Z_{n-1} Z_0 のボンドも含む）
    pbc=False: 開放境界条件
    """
    dim = 2**n
    H = np.zeros((dim, dim), dtype=complex)
    n_bonds = n if pbc else n - 1
    for i in range(n_bonds):
        H -= J * site_op(Z, i, n) @ site_op(Z, (i + 1) % n, n)
    for i in range(n):
        H -= h * site_op(X, i, n)
    return H


def magnetization_x(n):
    """M_x = (1/n) Σ X_j"""
    return sum(site_op(X, j, n) for j in range(n)) / n
