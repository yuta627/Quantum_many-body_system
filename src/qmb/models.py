"""模型ごとの定義：ハミルトニアンの行列表示、観測量、自由フェルミオン表示の行列 A, B など。"""
import numpy as np

from . import free_fermion as ff
from .operators import X, Z, site_op

# ---------------------------------------------------------------------
# 横磁場イジングモデル（TFIM）: H = -J Σ Z_i Z_{i+1} - h Σ X_i
# ---------------------------------------------------------------------


def tfim_hamiltonian(n, J, h, pbc=False):
    """TFIM のハミルトニアンの行列（2^n x 2^n）。

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


def tfim_free_fermion(n, J, h):
    """TFIM（開放境界）の Jordan-Wigner 変換後の行列 A, B と定数項を返す。

    規約: |+> = 空、|-> = 占有（X_j = 1 - 2 c_j† c_j）。
    H = Σ c† A c + (1/2)(c† B c† + h.c.) + const で
      A_jj = 2h, A_{j,j+1} = A_{j+1,j} = -J, B_{j,j+1} = -J, B_{j+1,j} = +J, const = -hn。
    """
    A = np.zeros((n, n))
    B = np.zeros((n, n))
    for j in range(n):
        A[j, j] = 2 * h
    for j in range(n - 1):
        A[j, j + 1] = A[j + 1, j] = -J
        B[j, j + 1] = -J
        B[j + 1, j] = J
    return A, B, -h * n


def tfim_ground_energy_ff(n, J, h):
    """TFIM（開放境界）の基底エネルギーを自由フェルミオン表示で求める。"""
    return ff.ground_energy(*tfim_free_fermion(n, J, h))


def tfim_mx_quench_ff(n, J, h, ts):
    """|+...+>（フェルミオンの真空）からのクエンチ後の m_x(t) = 1 - (2/n) Σ <n_j(t)>。"""
    A, B, _ = tfim_free_fermion(n, J, h)
    occ = ff.vacuum_occupations(ff.bdg_matrix(A, B), ts)
    return 1 - 2 / n * occ.sum(axis=1)
