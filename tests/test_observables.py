"""模型によらない観測量（src/qmb/observables.py）のテスト。"""
import numpy as np
import pytest

from qmb.ed import ground_state
from qmb.models import tfim_hamiltonian
from qmb.observables import energy_gap, entanglement_entropy, two_point_correlation
from qmb.operators import X, Z, plus_state


def test_energy_gap_n2():
    """n=2 の TFIM: 固有値は -sqrt(J^2+4h^2), -J, J, sqrt(J^2+4h^2)"""
    J, h = 1.0, 0.5
    H = tfim_hamiltonian(2, J, h)
    assert energy_gap(H) == pytest.approx(np.sqrt(J**2 + 4 * h**2) - J)


def test_entanglement_entropy_product_and_bell():
    # 積状態 |+...+> はもつれていない
    assert entanglement_entropy(plus_state(4), 4, 2) == pytest.approx(0.0, abs=1e-12)
    # ベル状態 (|00> + |11>)/sqrt(2) は S = log 2
    bell = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    assert entanglement_entropy(bell, 2, 1) == pytest.approx(np.log(2))


def test_two_point_correlation():
    # |+...+> では <X_i X_j> = 1、<Z_i Z_j> = 0
    psi = plus_state(4)
    assert two_point_correlation(psi, X, 0, X, 3, 4) == pytest.approx(1.0)
    assert two_point_correlation(psi, Z, 0, Z, 3, 4) == pytest.approx(0.0, abs=1e-12)
    # h=0 に近い強磁性の基底状態では <Z_0 Z_3> ≈ 1
    _, g = ground_state(tfim_hamiltonian(4, 1.0, 0.01))
    assert two_point_correlation(g, Z, 0, Z, 3, 4) > 0.99
