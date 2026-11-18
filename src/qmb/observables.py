"""模型によらない観測量：エネルギーギャップ、相関関数、もつれエントロピー。

状態ベクトルは np.kron の並び順（サイト0が最上位ビット）を前提とする。
"""
import numpy as np

from .operators import site_op


def energy_gap(H, k=1):
    """基底状態と k 番目の励起状態のエネルギー差 E_k - E_0（厳密対角化）。"""
    w = np.linalg.eigvalsh(H)
    return w[k] - w[0]


def two_point_correlation(psi, op_i, i, op_j, j, n):
    """2点相関関数 <ψ| A_i B_j |ψ>。例: Z_i Z_j の相関なら op_i = op_j = Z。"""
    O = site_op(op_i, i, n) @ site_op(op_j, j, n)
    return np.real(np.vdot(psi, O @ psi))


def entanglement_entropy(psi, n, cut):
    """サイト 0..cut-1 と cut..n-1 に分けたときのもつれエントロピー S = -Σ λ^2 log λ^2。

    状態ベクトルを (2^cut) x (2^(n-cut)) の行列に並べ替え、特異値分解（シュミット分解）する。
    """
    if not 0 < cut < n:
        raise ValueError("cut must satisfy 0 < cut < n")
    s = np.linalg.svd(psi.reshape(2**cut, 2**(n - cut)), compute_uv=False)
    p = s**2
    p = p[p > 1e-15]
    return float(-np.sum(p * np.log(p)))
