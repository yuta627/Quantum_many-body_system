"""自由フェルミオン（フェルミオンの2次形式）の一般的な扱い。模型によらない。

ハミルトニアンを
    H = Σ_ij c_i† A_ij c_j + (1/2) Σ_ij (c_i† B_ij c_j† + h.c.) + const
（A はエルミート行列、B は反対称行列）と書き、Nambu 表示 Ψ = (c_1..c_n, c_1†..c_n†)^T で
    H = (1/2) Ψ† H_BdG Ψ + (1/2) tr A + const,   H_BdG = [[A, B], [-B*, -A*]]
とする。模型ごとの A, B は qmb.models などで作る。
"""
import numpy as np


def bdg_matrix(A, B):
    """BdG 行列 H_BdG = [[A, B], [-B*, -A*]] を返す。"""
    A = np.asarray(A)
    B = np.asarray(B)
    return np.block([[A, B], [-B.conj(), -A.conj()]])


def quasiparticle_energies(H_bdg):
    """1粒子エネルギー ε_k ≥ 0（H_BdG の固有値は ±ε_k の対なので、大きい方の半分）。"""
    n = H_bdg.shape[0] // 2
    return np.linalg.eigvalsh(H_bdg)[n:]


def ground_energy(A, B, const=0.0):
    """基底エネルギー E0 = -(1/2) Σ ε_k + (1/2) tr A + const。"""
    eps = quasiparticle_energies(bdg_matrix(A, B))
    return -0.5 * np.sum(eps) + 0.5 * np.real(np.trace(A)) + const


def nambu_propagator(H_bdg, ts):
    """Heisenberg 描像の時間発展 Ψ(t) = U(t) Ψ、U(t) = e^{-i H_BdG t} を各時刻で返す。"""
    w, W = np.linalg.eigh(H_bdg)
    return np.array([(W * np.exp(-1j * w * t)) @ W.conj().T for t in ts])


def vacuum_occupations(H_bdg, ts):
    """フェルミオンの真空 |0> から時間発展したときの各サイトの粒子数 <n_j(t)>。

    c_j(t) = Σ_k [(U11)_jk c_k + (U12)_jk c_k†] より <0|n_j(t)|0> = Σ_k |(U12)_jk|^2。
    戻り値の形は (len(ts), n)。
    """
    n = H_bdg.shape[0] // 2
    U = nambu_propagator(H_bdg, ts)
    U12 = U[:, :n, n:]
    return np.sum(np.abs(U12) ** 2, axis=2)


def number_conserving_occupations(A, n0, ts):
    """粒子数を保存する2次形式（B = 0）で、Fock 状態から時間発展したときの <n_j(t)>。

    H = Σ c† A c（+ 定数）なら c_j(t) = Σ_k U_jk c_k、U(t) = e^{-iAt} なので、初期状態が
    各サイトの粒子数が決まった状態（n0[k] = 0 または 1）なら <n_j(t)> = Σ_k |U_jk(t)|^2 n0[k]。
    戻り値の形は (len(ts), n)。
    """
    w, W = np.linalg.eigh(np.asarray(A))
    n0 = np.asarray(n0, dtype=float)
    out = []
    for t in ts:
        U = (W * np.exp(-1j * w * t)) @ W.conj().T
        out.append(np.abs(U) ** 2 @ n0)
    return np.array(out)
