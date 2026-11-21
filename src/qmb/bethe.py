"""Bethe 仮説の方程式（対数形）を数値的に解く。模型によらない部分。

粒子（マグノンなど）M 個の Bethe 方程式を、ラピディティ λ_j を使った対数形
    n θ1(λ_j) = 2π I_j + Σ_{l≠j} θ2(λ_j - λ_l)      (j = 1, ..., M)
で書く（methods/06_bethe_ansatz を参照）。θ1 は1粒子の運動量に対応する関数、θ2 は
2粒子の散乱位相、I_j は整数または半整数の量子数である。模型ごとの θ1, θ2 と1粒子
エネルギーは qmb.models などで作る。
"""
import numpy as np
from scipy.optimize import root


def ground_state_quantum_numbers(m):
    """基底状態の量子数 I_j = -(m-1)/2, ..., (m-1)/2（間隔1で原点について対称）。"""
    return np.arange(m) - (m - 1) / 2


def bethe_residual(lam, n, I, theta1, theta2):
    """対数形の Bethe 方程式の残差 n θ1(λ_j) - 2π I_j - Σ_l θ2(λ_j - λ_l)。

    l = j の項は θ2(0) = 0 なので和に含めてよい。
    """
    lam = np.asarray(lam, dtype=float)
    return n * theta1(lam) - 2 * np.pi * np.asarray(I) - theta2(lam[:, None] - lam[None, :]).sum(axis=1)


def solve_bethe_equations(n, I, theta1, theta2, x0, theta1_inv=None, n_iter=0, mixing=0.5, tol=1e-12):
    """対数形の Bethe 方程式を解いてラピディティ λ_j を返す。

    x0          : 初期値
    theta1_inv  : θ1 の逆関数。与えると、まず反復
                  λ_j ← θ1^{-1}((2π I_j + Σ_l θ2(λ_j - λ_l)) / n)
                  を n_iter 回（mixing で混ぜながら）行ってから、ニュートン法系の解法で仕上げる。
    戻り値      : (λ, 残差の最大値)
    """
    lam = np.asarray(x0, dtype=float)
    I = np.asarray(I)
    if theta1_inv is not None:
        for _ in range(n_iter):
            new = theta1_inv((2 * np.pi * I + theta2(lam[:, None] - lam[None, :]).sum(axis=1)) / n)
            if np.max(np.abs(new - lam)) < tol:
                lam = new
                break
            lam = (1 - mixing) * lam + mixing * new
    sol = root(lambda x: bethe_residual(x, n, I, theta1, theta2), lam, method="hybr", options={"xtol": 1e-14})
    res = np.max(np.abs(bethe_residual(sol.x, n, I, theta1, theta2)))
    return sol.x, res
