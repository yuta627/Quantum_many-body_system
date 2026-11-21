"""厳密対角化。"""
import numpy as np


def ground_state(H):
    """基底エネルギー E0 と基底状態ベクトルを返す。"""
    w, v = np.linalg.eigh(H)
    return w[0], v[:, 0]


def time_evolve_expectation(H, psi0, A, ts):
    """<ψ(t)|A|ψ(t)> を各時刻で計算する（ψ(t) = e^{-iHt} ψ0）。

    1. H = U D U† と対角化（1回だけ）
    2. c = U† ψ0（固有基底で展開）
    3. ψ(t) = U e^{-iDt} c（位相を付けて戻す）
    4. 期待値 <ψ(t)|A|ψ(t)>
    """
    w, v = np.linalg.eigh(H)
    c = v.conj().T @ psi0
    out = []
    for t in ts:
        psi_t = v @ (np.exp(-1j * w * t) * c)
        out.append(np.real(psi_t.conj() @ A @ psi_t))
    return np.array(out)


# ---------------------------------------------------------------------
# 保存量によるブロック対角化（methods/01_exact_diag の4節）
# ---------------------------------------------------------------------


def sector_basis(n, m):
    """|1> の数がちょうど m 個の計算基底の状態の添字（昇順）を返す。

    添字は np.kron の並び順（サイト0が最上位ビット）での整数。たとえば n=2, m=1 なら
    |01> = 1, |10> = 2 なので [1, 2]。Z の総和 Σ_j Z_j = n - 2m が保存する模型では、
    H はこの添字の集合ごとのブロックに分かれる。
    """
    return np.array([s for s in range(2**n) if bin(s).count("1") == m], dtype=np.int64)


def project_to_sector(psi, basis):
    """全空間の状態ベクトル psi を、sector_basis で選んだ部分空間の成分だけにする。"""
    return np.asarray(psi)[basis]
