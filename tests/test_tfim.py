"""TFIM: 3つの手法（厳密対角化・JW変換・量子回路）が互いに一致することを確かめる。"""
import numpy as np
import pytest

from qmb import free_fermion as ff
from qmb.ed import ground_state, time_evolve_expectation
from qmb.models import (
    magnetization_x,
    tfim_free_fermion,
    tfim_ground_energy_ff,
    tfim_hamiltonian,
    tfim_mx_quench_ff,
)
from qmb.operators import plus_state

TS = np.linspace(0, 5, 51)


def test_n2_hand_calculation():
    """n=2: E0 = -sqrt(J^2+4h^2)、m_x(t) = cos^2(sqrt(2) t)（J=1, h=0.5）"""
    J, h = 1.0, 0.5
    E0, _ = ground_state(tfim_hamiltonian(2, J, h))
    assert E0 == pytest.approx(-np.sqrt(J**2 + 4 * h**2))
    mx = time_evolve_expectation(tfim_hamiltonian(2, J, h), plus_state(2), magnetization_x(2), TS)
    np.testing.assert_allclose(mx, np.cos(np.sqrt(2) * TS) ** 2, atol=1e-12)


@pytest.mark.parametrize("n,h", [(2, 0.5), (6, 0.7), (8, 1.3)])
def test_ed_vs_jordan_wigner(n, h):
    J = 1.0
    H = tfim_hamiltonian(n, J, h)
    assert ground_state(H)[0] == pytest.approx(tfim_ground_energy_ff(n, J, h), abs=1e-10)
    mx_ed = time_evolve_expectation(H, plus_state(n), magnetization_x(n), TS)
    np.testing.assert_allclose(tfim_mx_quench_ff(n, J, h, TS), mx_ed, atol=1e-10)


def test_bdg_spectrum_n2():
    """n=2: ε± = sqrt(J^2+4h^2) ± J"""
    J, h = 1.0, 0.5
    A, B, _ = tfim_free_fermion(2, J, h)
    eps = ff.quasiparticle_energies(ff.bdg_matrix(A, B))
    np.testing.assert_allclose(eps, [np.sqrt(2) - 1, np.sqrt(2) + 1])


qiskit = pytest.importorskip("qiskit")


def test_trotter_converges_and_orders_agree():
    from qmb.circuits import mx_trotter

    n, J, h, t_max = 6, 1.0, 0.7, 5.0
    ts, m1 = mx_trotter(n, J, h, t_max, 200, order=1)
    _, m2 = mx_trotter(n, J, h, t_max, 200, order=2)
    exact = time_evolve_expectation(tfim_hamiltonian(n, J, h), plus_state(n), magnetization_x(n), ts)
    # 1次と2次で m_x は完全に一致する
    np.testing.assert_allclose(m1, m2, atol=1e-12)
    # 誤差は O(δt^2)
    assert np.max(np.abs(m1 - exact)) < 2e-3


def test_vqe_n2_reaches_exact_energy():
    from qmb.circuits import vqe

    J, h = 1.0, 0.5
    energies, _ = vqe(2, J, h, p_max=1)
    assert energies[-1] == pytest.approx(-np.sqrt(J**2 + 4 * h**2), abs=1e-6)
