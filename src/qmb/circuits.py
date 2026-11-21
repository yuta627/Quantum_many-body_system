"""量子回路（Qiskit）による Trotter 時間発展と VQE（TFIM・XXZ 模型）。

TFIM: Trotter 時間発展は開放境界条件。VQE（tfim_op, hva_circuit, vqe）は pbc=True で周期境界条件も扱える。
XXZ : xxz_ で始まる関数。Trotter 時間発展は開放境界条件、VQE は pbc=True で周期境界条件も扱える。

量子ビット j はサイト j に対応させる（Qiskit の状態ベクトルの並び順は np.kron と逆だが、
演算子もゲートもサイトの番号で指定しているので、期待値はそのまま厳密対角化と比較できる）。
"""
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
from qiskit.quantum_info import SparsePauliOp, Statevector
from scipy.optimize import minimize


def _bonds(n, pbc=False):
    """ZZ 相互作用のボンドの一覧。pbc=True なら (n-1, 0) も含む。"""
    bonds = [(j, j + 1) for j in range(n - 1)]
    if pbc:
        bonds.append((n - 1, 0))
    return bonds


def tfim_op(n, J, h, pbc=False):
    """H = -J Σ Z_j Z_{j+1} - h Σ X_j を SparsePauliOp で返す。

    pbc=True: 周期境界条件（Z_{n-1} Z_0 のボンドも含む）
    """
    zz = [("ZZ", [a, b], -J) for a, b in _bonds(n, pbc)]
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


def hva_circuit(n, J, h, p, local=False, pbc=False):
    """Hamiltonian Variational Ansatz（e^{-iγ H_X} e^{-iβ H_ZZ} を p 層）。

    local=False: 各層で β, γ を全サイト共通（パラメータ 2p 個）
    local=True : ゲートごとに別のパラメータ（パラメータ p(n_bonds + n) 個。
                 開放境界では n_bonds = n-1 なので p(2n-1) 個）
    pbc=True   : 周期境界条件（ボンド (n-1, 0) の R_ZZ も入れる）
    """
    bonds = _bonds(n, pbc)
    nb, ns = (len(bonds), n) if local else (1, 1)
    theta = ParameterVector("θ", p * (nb + ns))
    qc = plus_circuit(n)
    k = 0
    for _ in range(p):
        beta = theta[k:k + nb]
        k += nb
        gamma = theta[k:k + ns]
        k += ns
        # 偶数番目のボンド → 奇数番目のボンド（周期境界の (n-1, 0) は最後）の順にかける
        order = list(range(0, len(bonds), 2)) + list(range(1, len(bonds), 2))
        for i in order:
            a, b = bonds[i]
            qc.rzz(-2 * J * beta[i if local else 0], a, b)
        for j in range(n):
            qc.rx(-2 * h * gamma[j if local else 0], j)
    return qc


def energy(qc, H_op, params):
    """エネルギー期待値 E(θ)。"""
    return Statevector(qc.assign_parameters(params)).expectation_value(H_op).real


def vqe(n, J, h, p_max, local=False, pbc=False):
    """層を1つずつ増やし、前の層の最適値を初期値にして VQE を行う。

    pbc=True: 周期境界条件
    戻り値: (各 p での最適エネルギーの配列, 最後のパラメータ)
    """
    H_op = tfim_op(n, J, h, pbc)
    n_per_layer = (len(_bonds(n, pbc)) + n) if local else 2
    params = np.array([])
    energies = []
    for p in range(1, p_max + 1):
        qc = hva_circuit(n, J, h, p, local, pbc)
        x0 = np.concatenate([params, np.full(n_per_layer, 0.01)])
        res = minimize(lambda x: energy(qc, H_op, x), x0, method="BFGS")
        params = res.x
        energies.append(res.fun)
    return np.array(energies), params


# ---------------------------------------------------------------------
# XXZ 模型: H = J Σ (X_j X_{j+1} + Y_j Y_{j+1} + Δ Z_j Z_{j+1})
# ---------------------------------------------------------------------


def _xxz_bonds(n, pbc=False):
    """ボンドの一覧。pbc=True なら (n-1, 0) も含む（n=2 では同じボンドを2回数える。qmb.models と同じ）。"""
    bonds = [(j, j + 1) for j in range(n - 1)]
    if pbc:
        bonds.append((n - 1, 0))
    return bonds


def _even_odd_bonds(n, pbc=False):
    """偶数ボンド (0,1), (2,3), ... と奇数ボンド (1,2), (3,4), ...（周期境界なら (n-1,0) も奇数側）に分ける。"""
    bonds = _xxz_bonds(n, pbc)
    even = [b for b in bonds if b[0] % 2 == 0 and b[1] == b[0] + 1]
    odd = [b for b in bonds if b not in even]
    return even, odd


def xxz_op(n, J, Delta, pbc=False):
    """H = J Σ (X_j X_{j+1} + Y_j Y_{j+1} + Δ Z_j Z_{j+1}) を SparsePauliOp で返す。"""
    terms = []
    for a, b in _xxz_bonds(n, pbc):
        terms += [("XX", [a, b], J), ("YY", [a, b], J), ("ZZ", [a, b], J * Delta)]
    return SparsePauliOp.from_sparse_list(terms, num_qubits=n)


def staggered_op(n):
    """M_s = (1/n) Σ (-1)^j Z_j（量子ビット j = サイト j）"""
    return SparsePauliOp.from_sparse_list([("Z", [j], (-1) ** j / n) for j in range(n)], num_qubits=n)


def neel_circuit(n):
    """Néel 状態 |0101...>（奇数サイトに X をかける）を作る回路。"""
    qc = QuantumCircuit(n)
    for j in range(1, n, 2):
        qc.x(j)
    return qc


def _xxz_bond_gate(qc, a, b, J, Delta, tau):
    """e^{-iτJ(X_aX_b + Y_aY_b + ΔZ_aZ_b)} = R_XX(2Jτ) R_YY(2Jτ) R_ZZ(2JΔτ)（3つの項は交換するので厳密）。"""
    qc.rxx(2 * J * tau, a, b)
    qc.ryy(2 * J * tau, a, b)
    qc.rzz(2 * J * Delta * tau, a, b)


def xxz_trotter_step(n, J, Delta, dt, order=1):
    """XXZ 模型（開放境界）の Trotter 分解の1ステップ。偶数ボンドと奇数ボンドの2つに分ける。

    order=1: e^{-iH_odd dt} e^{-iH_even dt}
    order=2: e^{-iH_even dt/2} e^{-iH_odd dt} e^{-iH_even dt/2}
    """
    qc = QuantumCircuit(n)
    even, odd = _even_odd_bonds(n)

    def layer(bonds, tau):
        for a, b in bonds:
            _xxz_bond_gate(qc, a, b, J, Delta, tau)

    if order == 1:
        layer(even, dt)
        layer(odd, dt)
    else:
        layer(even, dt / 2)
        layer(odd, dt)
        layer(even, dt / 2)
    return qc


def ms_trotter(n, J, Delta, t_max, r, order=1):
    """Trotter 分解で、Néel 状態からの m_s(t) を計算する（t = 0, δt, ..., t_max）。"""
    dt = t_max / r
    step = xxz_trotter_step(n, J, Delta, dt, order)
    sv = Statevector(neel_circuit(n))
    M = staggered_op(n)
    out = [sv.expectation_value(M).real]
    for _ in range(r):
        sv = sv.evolve(step)
        out.append(sv.expectation_value(M).real)
    return np.linspace(0, t_max, r + 1), np.array(out)


def singlet_circuit(n):
    """偶数ボンド (0,1), (2,3), ... に一重項 (|01> - |10>)/√2 を並べた状態を作る回路（n は偶数）。

    各ボンドで X_a → H_a → CX(a, b) → X_b とかける：|00> → |10> → (|00> - |10>)/√2
    → (|00> - |11>)/√2 → (|01> - |10>)/√2。
    """
    if n % 2:
        raise ValueError("n must be even")
    qc = QuantumCircuit(n)
    for a in range(0, n, 2):
        qc.x(a)
        qc.h(a)
        qc.cx(a, a + 1)
        qc.x(a + 1)
    return qc


def xxz_hva_circuit(n, J, Delta, p, local=False, pbc=False):
    """XXZ 模型の Hamiltonian Variational Ansatz。初期状態は偶数ボンドの一重項の積。

    1層は e^{-iβ4 H_even^XY} e^{-iβ3 H_even^ZZ} e^{-iβ2 H_odd^XY} e^{-iβ1 H_odd^ZZ}
    （H^XY = J Σ (XX + YY)、H^ZZ = JΔ Σ ZZ）。これを p 層重ねる。
    local=False: 各層で4つの角度を全ボンド共通（パラメータ 4p 個）
    local=True : ゲートごとに別のパラメータ（パラメータ 2p × ボンドの数 個）
    pbc=True   : 周期境界条件（ボンド (n-1, 0) を奇数側に入れる）
    """
    even, odd = _even_odd_bonds(n, pbc)
    groups = [(odd, "zz"), (odd, "xy"), (even, "zz"), (even, "xy")]
    sizes = [len(g) if local else 1 for g, _ in groups]
    theta = ParameterVector("θ", p * sum(sizes))
    qc = singlet_circuit(n)
    k = 0
    for _ in range(p):
        for (bonds, kind), size in zip(groups, sizes):
            th = theta[k:k + size]
            k += size
            for i, (a, b) in enumerate(bonds):
                angle = th[i if local else 0]
                if kind == "zz":
                    qc.rzz(2 * J * Delta * angle, a, b)
                else:
                    qc.rxx(2 * J * angle, a, b)
                    qc.ryy(2 * J * angle, a, b)
    return qc


def xxz_vqe(n, J, Delta, p_max, local=False, pbc=False):
    """XXZ 模型の VQE。層を1つずつ増やし、前の層の最適値を初期値にする。

    戻り値: (各 p での最適エネルギーの配列, 最後のパラメータ)
    """
    H_op = xxz_op(n, J, Delta, pbc)
    even, odd = _even_odd_bonds(n, pbc)
    n_per_layer = 2 * (len(even) + len(odd)) if local else 4
    params = np.array([])
    energies = []
    for p in range(1, p_max + 1):
        qc = xxz_hva_circuit(n, J, Delta, p, local, pbc)
        x0 = np.concatenate([params, np.full(n_per_layer, 0.01)])
        res = minimize(lambda x: energy(qc, H_op, x), x0, method="BFGS")
        params = res.x
        energies.append(res.fun)
    return np.array(energies), params
