"""自由フェルミオンの一般的な関数（src/qmb/free_fermion.py）のテスト。

ランダムな 2次形式（A: エルミート、B: 反対称）について、フォック空間で作ったハミルトニアンの
厳密対角化と、BdG 行列による結果が一致することを確かめる。
"""
from functools import reduce

import numpy as np
import pytest

from qmb import free_fermion as ff

N = 4
I2 = np.eye(2)
ZM = np.diag([1.0, -1.0])
A_LOWER = np.array([[0, 1], [0, 0]], dtype=complex)  # |0><1|（|0> = 空、|1> = 占有）


def fermion_ops(n):
    """Jordan-Wigner 変換で作った消滅演算子 c_0, ..., c_{n-1}（フォック空間の行列）。"""
    return [reduce(np.kron, [ZM] * j + [A_LOWER] + [I2] * (n - j - 1)) for j in range(n)]


def random_quadratic(n, seed=0):
    rng = np.random.default_rng(seed)
    M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    A = (M + M.conj().T) / 2
    K = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    B = (K - K.T) / 2
    return A, B


def fock_hamiltonian(A, B, const=0.0):
    c = fermion_ops(A.shape[0])
    cd = [x.conj().T for x in c]
    n = len(c)
    H = const * np.eye(2**n, dtype=complex)
    for i in range(n):
        for j in range(n):
            H += A[i, j] * cd[i] @ c[j]
            H += 0.5 * (B[i, j] * cd[i] @ cd[j] + np.conj(B[i, j]) * c[j] @ c[i])
    return H


def test_ground_energy_random_quadratic():
    A, B = random_quadratic(N)
    const = 0.3
    E_fock = np.linalg.eigvalsh(fock_hamiltonian(A, B, const))[0]
    assert ff.ground_energy(A, B, const) == pytest.approx(E_fock, abs=1e-10)


def test_vacuum_occupations_random_quadratic():
    A, B = random_quadratic(N, seed=1)
    H = fock_hamiltonian(A, B)
    c = fermion_ops(N)
    vac = np.zeros(2**N, dtype=complex)
    vac[0] = 1.0  # 全サイト空
    ts = np.linspace(0, 2, 11)
    occ = ff.vacuum_occupations(ff.bdg_matrix(A, B), ts)
    w, v = np.linalg.eigh(H)
    for k, t in enumerate(ts):
        psi = v @ (np.exp(-1j * w * t) * (v.conj().T @ vac))
        exact = [np.real(np.vdot(psi, c[j].conj().T @ c[j] @ psi)) for j in range(N)]
        np.testing.assert_allclose(occ[k], exact, atol=1e-10)
