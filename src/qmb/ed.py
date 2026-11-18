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
