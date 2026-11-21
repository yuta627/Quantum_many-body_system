"""XXZ 模型: 厳密対角化・Jordan-Wigner 変換（Δ=0）・Bethe 仮説・量子回路が互いに一致することを確かめる。"""
import numpy as np
import pytest
from scipy.special import j0

from qmb.ed import ground_state, project_to_sector, sector_basis, time_evolve_expectation
from qmb.models import (
    staggered_magnetization,
    staggered_magnetization_sector,
    xxz_bethe_ground_state,
    xxz_e0_thermodynamic,
    xxz_free_fermion,
    xxz_hamiltonian,
    xxz_hamiltonian_sector,
    xxz_ms_quench_ff,
)
from qmb.free_fermion import ground_energy
from qmb.operators import neel_state

TS = np.linspace(0, 3, 31)


def test_n2_hand_calculation():
    """n=2: E0 = -J(2+Δ)、Néel 状態からの m_s(t) = cos(4Jt)（Δ によらない）"""
    J = 1.0
    for Delta in [-0.5, 0.0, 0.7, 2.0]:
        H = xxz_hamiltonian(2, J, Delta)
        assert ground_state(H)[0] == pytest.approx(-J * (2 + Delta))
        ms = time_evolve_expectation(H, neel_state(2), staggered_magnetization(2), TS)
        np.testing.assert_allclose(ms, np.cos(4 * J * TS), atol=1e-12)


@pytest.mark.parametrize("n,pbc", [(4, False), (5, True), (6, True)])
def test_sector_blocks_reproduce_spectrum(n, pbc):
    J, Delta = 1.0, 0.7
    w = np.linalg.eigvalsh(xxz_hamiltonian(n, J, Delta, pbc))
    ws = np.concatenate([np.linalg.eigvalsh(xxz_hamiltonian_sector(n, J, Delta, m, pbc).toarray()) for m in range(n + 1)])
    np.testing.assert_allclose(np.sort(ws), w, atol=1e-10)


def test_sector_quench_matches_full_space():
    n, J, Delta = 6, 1.0, 1.5
    m = n // 2
    full = time_evolve_expectation(xxz_hamiltonian(n, J, Delta), neel_state(n), staggered_magnetization(n), TS)
    psi0 = project_to_sector(neel_state(n), sector_basis(n, m))
    sec = time_evolve_expectation(
        xxz_hamiltonian_sector(n, J, Delta, m).toarray(), psi0, staggered_magnetization_sector(n, m).toarray(), TS)
    np.testing.assert_allclose(sec, full, atol=1e-12)


@pytest.mark.parametrize("n", [2, 6, 9])
def test_jordan_wigner_delta0(n):
    J = 1.0
    H = xxz_hamiltonian(n, J, 0.0)
    assert ground_state(H)[0] == pytest.approx(ground_energy(*xxz_free_fermion(n, J)), abs=1e-10)
    ms_ed = time_evolve_expectation(H, neel_state(n), staggered_magnetization(n), TS)
    np.testing.assert_allclose(xxz_ms_quench_ff(n, J, TS), ms_ed, atol=1e-10)


def test_jordan_wigner_delta0_bessel():
    """Δ=0 の長い鎖の中央付近では m_s(t) ≈ J_0(8Jt)（短い時間）"""
    ts = np.linspace(0, 2, 21)
    np.testing.assert_allclose(xxz_ms_quench_ff(400, 1.0, ts), j0(8 * ts), atol=1e-2)


@pytest.mark.parametrize("Delta", [-0.6, 0.0, 0.5, 1.0, 2.0])
@pytest.mark.parametrize("n", [2, 6, 10])
def test_bethe_vs_ed(n, Delta):
    J = 1.0
    E_ed = ground_state(xxz_hamiltonian(n, J, Delta, pbc=True))[0]
    E_bethe, _, res = xxz_bethe_ground_state(n, J, Delta)
    assert res < 1e-9
    assert E_bethe == pytest.approx(E_ed, abs=1e-9)


@pytest.mark.parametrize("Delta", [-0.5, 0.0, 0.5, 1.0, 2.0])
def test_bethe_thermodynamic_limit(Delta):
    E, _, _ = xxz_bethe_ground_state(400, 1.0, Delta)
    assert E / 400 == pytest.approx(xxz_e0_thermodynamic(1.0, Delta), abs=1e-4)


def test_thermodynamic_limit_known_values():
    assert xxz_e0_thermodynamic(1.0, 0.0) == pytest.approx(-4 / np.pi)
    assert xxz_e0_thermodynamic(1.0, 0.5) == pytest.approx(-1.5)
    assert xxz_e0_thermodynamic(1.0, 1.0) == pytest.approx(1 - 4 * np.log(2))


qiskit = pytest.importorskip("qiskit")


def test_trotter_n2_is_exact():
    """n=2 はボンドが1本だけなので、Trotter 分解しても誤差がない。"""
    from qmb.circuits import ms_trotter

    ts, ms = ms_trotter(2, 1.0, 0.7, 3.0, 7)
    np.testing.assert_allclose(ms, np.cos(4 * ts), atol=1e-12)


def test_trotter_converges():
    from qmb.circuits import ms_trotter

    n, J, Delta, t_max = 6, 1.0, 0.5, 2.0
    ts, m1 = ms_trotter(n, J, Delta, t_max, 200, order=1)
    _, m2 = ms_trotter(n, J, Delta, t_max, 200, order=2)
    exact = time_evolve_expectation(xxz_hamiltonian(n, J, Delta), neel_state(n), staggered_magnetization(n), ts)
    assert np.max(np.abs(m1 - exact)) < 5e-3
    assert np.max(np.abs(m2 - exact)) < np.max(np.abs(m1 - exact))


def test_circuit_operators_match_matrices():
    from qiskit.quantum_info import Statevector

    from qmb.circuits import neel_circuit, staggered_op, xxz_op

    n, J, Delta = 5, 1.0, 0.7
    E_op = np.linalg.eigvalsh(xxz_op(n, J, Delta).to_matrix())
    np.testing.assert_allclose(E_op, np.linalg.eigvalsh(xxz_hamiltonian(n, J, Delta)), atol=1e-10)
    assert Statevector(neel_circuit(n)).expectation_value(staggered_op(n)).real == pytest.approx(
        np.real(np.vdot(neel_state(n), staggered_magnetization(n) @ neel_state(n))))


def test_vqe_singlet_is_exact_for_n2():
    from qiskit.quantum_info import Statevector

    from qmb.circuits import singlet_circuit, xxz_op

    J, Delta = 1.0, 0.7
    E = Statevector(singlet_circuit(2)).expectation_value(xxz_op(2, J, Delta)).real
    assert E == pytest.approx(-J * (2 + Delta))


def test_vqe_n4_periodic():
    from qmb.circuits import xxz_vqe

    J, Delta = 1.0, 0.5
    energies, _ = xxz_vqe(4, J, Delta, p_max=2, pbc=True)
    assert energies[-1] == pytest.approx(ground_state(xxz_hamiltonian(4, J, Delta, pbc=True))[0], abs=1e-6)
