"""Jordan-Wigner 変換による自由フェルミオン解（TFIM、開放境界条件）。

規約: |+> = フェルミオンなし、|-> = フェルミオンあり（X_j = 1 - 2 c_j† c_j）。
"""
import numpy as np


def tfim_bdg_matrix(n, J, h):
    """BdG 行列 H_BdG = [[A, B], [-B, -A]]（H = (1/2) Ψ† H_BdG Ψ）。"""
    A = np.zeros((n, n))
    B = np.zeros((n, n))
    for j in range(n):
        A[j, j] = 2 * h
    for j in range(n - 1):
        A[j, j + 1] = A[j + 1, j] = -J
        B[j, j + 1] = -J
        B[j + 1, j] = J
    return np.block([[A, B], [-B, -A]])


def quasiparticle_energies(n, J, h):
    """1粒子エネルギー ε_k ≥ 0（H_BdG の正の固有値）。"""
    e = np.linalg.eigvalsh(tfim_bdg_matrix(n, J, h))
    return e[n:]


def ground_energy(n, J, h):
    """E0 = -(1/2) Σ ε_k"""
    return -0.5 * np.sum(quasiparticle_energies(n, J, h))


def mx_quench(n, J, h, ts):
    """|+...+> からのクエンチ後の m_x(t) = 1 - (2/n) Σ |U12(t)|^2。"""
    w, W = np.linalg.eigh(tfim_bdg_matrix(n, J, h))
    out = []
    for t in ts:
        U = (W * np.exp(-1j * w * t)) @ W.conj().T
        U12 = U[:n, n:]
        out.append(1 - 2 / n * np.sum(np.abs(U12) ** 2))
    return np.array(out)
