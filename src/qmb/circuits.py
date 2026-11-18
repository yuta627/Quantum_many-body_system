"""量子回路（Qiskit）による Trotter 時間発展と VQE（TFIM、開放境界条件）。"""
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
from qiskit.quantum_info import SparsePauliOp, Statevector
from scipy.optimize import minimize


def tfim_op(n, J, h):
    """H = -J Σ Z_j Z_{j+1} - h Σ X_j を SparsePauliOp で返す。"""
    zz = [("ZZ", [j, j + 1], -J) for j in range(n - 1)]
    x = [("X", [j], -h) for j in range(n)]
    return SparsePauliOp.from_sparse_list(zz + x, num_qubits=n)


def mx_op(n):
    """M_x = (1/n) Σ X_j"""
    return SparsePauliOp.from_sparse_list([("X", [j], 1 / n) for j in range(n)], num_qubits=n)


def plus_circuit(n):
    """|+...+> を作る回路。"""
    qc = QuantumCircuit(n)
    qc.h(range(n))
    return qc


def trotter_step(n, J, h, dt, order=1):
    """Trotter 分解の1ステップ。

    order=1: e^{-iH_X dt} e^{-iH_ZZ dt}
    order=2: e^{-iH_X dt/2} e^{-iH_ZZ dt} e^{-iH_X dt/2}
    """
    qc = QuantumCircuit(n)

    def zz_layer(tau):  # e^{-iH_ZZ tau} = Π R_ZZ(-2J tau)
        for start in (0, 1):
            for j in range(start, n - 1, 2):
                qc.rzz(-2 * J * tau, j, j + 1)

    def x_layer(tau):  # e^{-iH_X tau} = Π R_X(-2h tau)
        for j in range(n):
            qc.rx(-2 * h * tau, j)

    if order == 1:
        zz_layer(dt)
        x_layer(dt)
    else:
        x_layer(dt / 2)
        zz_layer(dt)
        x_layer(dt / 2)
    return qc


def mx_trotter(n, J, h, t_max, r, order=1):
    """Trotter 分解で m_x(t) を計算する（t = 0, δt, ..., t_max）。"""
    dt = t_max / r
    step = trotter_step(n, J, h, dt, order)
    sv = Statevector(plus_circuit(n))
    M = mx_op(n)
    out = [sv.expectation_value(M).real]
    for _ in range(r):
        sv = sv.evolve(step)
        out.append(sv.expectation_value(M).real)
    return np.linspace(0, t_max, r + 1), np.array(out)


def hva_circuit(n, J, h, p, local=False):
    """Hamiltonian Variational Ansatz（e^{-iγ H_X} e^{-iβ H_ZZ} を p 層）。

    local=False: 各層で β, γ を全サイト共通（パラメータ 2p 個）
    local=True : ゲートごとに別のパラメータ（パラメータ p(2n-1) 個）
    """
    nb, ns = (n - 1, n) if local else (1, 1)
    theta = ParameterVector("θ", p * (nb + ns))
    qc = plus_circuit(n)
    k = 0
    for _ in range(p):
        beta = theta[k:k + nb]
        k += nb
        gamma = theta[k:k + ns]
        k += ns
        for start in (0, 1):
            for j in range(start, n - 1, 2):
                qc.rzz(-2 * J * beta[j if local else 0], j, j + 1)
        for j in range(n):
            qc.rx(-2 * h * gamma[j if local else 0], j)
    return qc


def energy(qc, H_op, params):
    """エネルギー期待値 E(θ)。"""
    return Statevector(qc.assign_parameters(params)).expectation_value(H_op).real


def vqe(n, J, h, p_max, local=False):
    """層を1つずつ増やし、前の層の最適値を初期値にして VQE を行う。

    戻り値: (各 p での最適エネルギーの配列, 最後のパラメータ)
    """
    H_op = tfim_op(n, J, h)
    n_per_layer = (2 * n - 1) if local else 2
    params = np.array([])
    energies = []
    for p in range(1, p_max + 1):
        qc = hva_circuit(n, J, h, p, local)
        x0 = np.concatenate([params, np.full(n_per_layer, 0.01)])
        res = minimize(lambda x: energy(qc, H_op, x), x0, method="BFGS")
        params = res.x
        energies.append(res.fun)
    return np.array(energies), params
