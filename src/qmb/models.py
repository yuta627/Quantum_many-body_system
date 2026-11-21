"""模型ごとの定義：ハミルトニアンの行列表示、観測量、自由フェルミオン表示の行列 A, B など。"""
import numpy as np
from scipy import sparse
from scipy.integrate import quad

from . import bethe
from . import free_fermion as ff
from .ed import sector_basis
from .operators import X, Y, Z, site_op

# ---------------------------------------------------------------------
# 横磁場イジングモデル（TFIM）: H = -J Σ Z_i Z_{i+1} - h Σ X_i
# ---------------------------------------------------------------------


def tfim_hamiltonian(n, J, h, pbc=False):
    """TFIM のハミルトニアンの行列（2^n x 2^n）。

    pbc=True: 周期境界条件（Z_{n-1} Z_0 のボンドも含む）
    pbc=False: 開放境界条件
    """
    dim = 2**n
    H = np.zeros((dim, dim), dtype=complex)
    n_bonds = n if pbc else n - 1
    for i in range(n_bonds):
        H -= J * site_op(Z, i, n) @ site_op(Z, (i + 1) % n, n)
    for i in range(n):
        H -= h * site_op(X, i, n)
    return H


def magnetization_x(n):
    """M_x = (1/n) Σ X_j"""
    return sum(site_op(X, j, n) for j in range(n)) / n


def tfim_free_fermion(n, J, h):
    """TFIM（開放境界）の Jordan-Wigner 変換後の行列 A, B と定数項を返す。

    規約: |+> = 空、|-> = 占有（X_j = 1 - 2 c_j† c_j）。
    H = Σ c† A c + (1/2)(c† B c† + h.c.) + const で
      A_jj = 2h, A_{j,j+1} = A_{j+1,j} = -J, B_{j,j+1} = -J, B_{j+1,j} = +J, const = -hn。
    """
    A = np.zeros((n, n))
    B = np.zeros((n, n))
    for j in range(n):
        A[j, j] = 2 * h
    for j in range(n - 1):
        A[j, j + 1] = A[j + 1, j] = -J
        B[j, j + 1] = -J
        B[j + 1, j] = J
    return A, B, -h * n


def tfim_ground_energy_ff(n, J, h):
    """TFIM（開放境界）の基底エネルギーを自由フェルミオン表示で求める。"""
    return ff.ground_energy(*tfim_free_fermion(n, J, h))


def tfim_mx_quench_ff(n, J, h, ts):
    """|+...+>（フェルミオンの真空）からのクエンチ後の m_x(t) = 1 - (2/n) Σ <n_j(t)>。"""
    A, B, _ = tfim_free_fermion(n, J, h)
    occ = ff.vacuum_occupations(ff.bdg_matrix(A, B), ts)
    return 1 - 2 / n * occ.sum(axis=1)


# ---------------------------------------------------------------------
# XXZ 模型: H = J Σ (X_i X_{i+1} + Y_i Y_{i+1} + Δ Z_i Z_{i+1})
# ---------------------------------------------------------------------


def _xxz_bonds(n, pbc):
    """ボンドの一覧。pbc=True なら (n-1, 0) も含む（n=2 では同じボンドを2回数える）。"""
    return [(i, (i + 1) % n) for i in range(n if pbc else n - 1)]


def xxz_hamiltonian(n, J, Delta, pbc=False):
    """XXZ 模型のハミルトニアンの行列（2^n x 2^n）。

    pbc=True: 周期境界条件（ボンド (n-1, 0) も含む）
    pbc=False: 開放境界条件
    """
    dim = 2**n
    H = np.zeros((dim, dim), dtype=complex)
    bond = J * (np.kron(X, X) + np.kron(Y, Y) + Delta * np.kron(Z, Z))   # 隣り合う2サイトの 4x4 行列
    for a, b in _xxz_bonds(n, pbc):
        if b == a + 1:
            # I ⊗ ... ⊗ (bond) ⊗ ... ⊗ I（行列の積を使わずクロネッカー積だけで作る）
            H += np.kron(np.kron(np.eye(2**a), bond), np.eye(2 ** (n - a - 2)))
        else:  # 周期境界のボンド (n-1, 0)
            H += J * (site_op(X, a, n) @ site_op(X, b, n) + site_op(Y, a, n) @ site_op(Y, b, n))
            H += J * Delta * site_op(Z, a, n) @ site_op(Z, b, n)
    return H


def staggered_magnetization(n):
    """スタッガード磁化 M_s = (1/n) Σ_j (-1)^j Z_j"""
    return sum((-1) ** j * site_op(Z, j, n) for j in range(n)) / n


def _z_values(basis, n):
    """各基底状態での Z_j の値（+1 または -1）を並べた配列。形は (len(basis), n)。"""
    bits = (basis[:, None] >> (n - 1 - np.arange(n))[None, :]) & 1
    return 1 - 2 * bits


def xxz_hamiltonian_sector(n, J, Delta, m, pbc=False):
    """XXZ 模型のハミルトニアンを、|1> が m 個の部分空間（Σ Z_j = n - 2m）に制限した疎行列。

    基底は qmb.ed.sector_basis(n, m) の順。XX + YY = 2(σ+σ- + σ-σ+) なので、
    ボンドの両端が異なる（01 または 10）ときだけ、両方を反転させる行列要素 2J がある。
    """
    basis = sector_basis(n, m)
    z = _z_values(basis, n)
    rows, cols, vals = [], [], []
    diag = np.zeros(len(basis))
    for a, b in _xxz_bonds(n, pbc):
        diag += J * Delta * z[:, a] * z[:, b]
        flip = z[:, a] != z[:, b]
        src = np.nonzero(flip)[0]
        dst = np.searchsorted(basis, basis[src] ^ ((1 << (n - 1 - a)) | (1 << (n - 1 - b))))
        rows.append(dst)
        cols.append(src)
        vals.append(np.full(len(src), 2.0 * J))
    rows = np.concatenate(rows + [np.arange(len(basis))])
    cols = np.concatenate(cols + [np.arange(len(basis))])
    vals = np.concatenate(vals + [diag])
    return sparse.csr_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)))


def staggered_magnetization_sector(n, m):
    """M_s = (1/n) Σ_j (-1)^j Z_j を、|1> が m 個の部分空間に制限した対角行列（疎行列）。"""
    z = _z_values(sector_basis(n, m), n)
    return sparse.diags(z @ ((-1.0) ** np.arange(n)) / n)


def xxz_free_fermion(n, J):
    """Δ=0（XX 模型、開放境界）の Jordan-Wigner 変換後の行列 A, B と定数項を返す。

    規約: |0> = 空、|1> = 占有（Z_j = 1 - 2 c_j† c_j）。
    X_j X_{j+1} + Y_j Y_{j+1} = 2(c_j† c_{j+1} + c_{j+1}† c_j) より
      A_{j,j+1} = A_{j+1,j} = 2J, B = 0, const = 0。
    Δ≠0 では Z_j Z_{j+1} が4次の項になるので、この形には書けない。
    """
    A = np.zeros((n, n))
    for j in range(n - 1):
        A[j, j + 1] = A[j + 1, j] = 2 * J
    return A, np.zeros((n, n)), 0.0


def xxz_ms_quench_ff(n, J, ts):
    """Δ=0 で Néel 状態 |0101...> から時間発展したときの m_s(t)（Jordan-Wigner 変換）。

    奇数サイトだけが占有された Fock 状態から出発し、Z_j = 1 - 2<n_j(t)> を使う。
    """
    A, _, _ = xxz_free_fermion(n, J)
    n0 = np.arange(n) % 2
    occ = ff.number_conserving_occupations(A, n0, ts)
    signs = (-1.0) ** np.arange(n)
    return (1 - 2 * occ) @ signs / n


def xxz_bethe_functions(J, Delta):
    """XXZ 模型（周期境界、J>0）の Bethe 方程式に使う関数 (θ1, θ2, θ1 の逆関数, 1粒子エネルギー ε)。

    -1 < Δ < 1 : Δ = cos γ、θ_m(λ) = 2 arctan[cot(mγ/2) tanh λ]、ε(λ) = -4J sin^2γ / (cosh 2λ - cos γ)
    Δ = 1      : θ_m(λ) = 2 arctan(2λ/m)、ε(λ) = -2J / (λ^2 + 1/4)
    Δ > 1      : Δ = cosh η、θ_m(λ) = 2 arctan[coth(mη/2) tan λ]、ε(λ) = -4J sinh^2η / (cosh η - cos 2λ)
    マグノンの運動量は k = π - θ1(λ)。導出は models/02_xxz/02_jordan_wigner.ipynb を参照。
    """
    if not -1 < Delta:
        raise ValueError("Delta must satisfy Delta > -1")
    if Delta < 1:
        g = np.arccos(Delta)

        def theta(lam, m):
            return 2 * np.arctan(np.tanh(lam) / np.tan(m * g / 2))

        def theta1_inv(t):
            t = np.clip(t, -(np.pi - g) + 1e-12, (np.pi - g) - 1e-12)
            return np.arctanh(np.tan(g / 2) * np.tan(t / 2))

        def eps(lam):
            return -4 * J * np.sin(g) ** 2 / (np.cosh(2 * lam) - np.cos(g))

    elif Delta == 1:

        def theta(lam, m):
            return 2 * np.arctan(2 * lam / m)

        def theta1_inv(t):
            return np.tan(np.clip(t, -np.pi + 1e-12, np.pi - 1e-12) / 2) / 2

        def eps(lam):
            return -2 * J / (lam**2 + 0.25)

    else:
        e = np.arccosh(Delta)

        def theta(lam, m):
            # tan λ の不連続（λ = ±π/2）をまたいでも連続になるよう 2π の整数倍を足す
            return 2 * np.arctan(np.tan(lam) / np.tanh(m * e / 2)) + 2 * np.pi * np.floor(lam / np.pi + 0.5)

        def theta1_inv(t):
            return np.arctan(np.tanh(e / 2) * np.tan(t / 2))

        def eps(lam):
            return -4 * J * np.sinh(e) ** 2 / (np.cosh(e) - np.cos(2 * lam))

    return (lambda lam: theta(lam, 1)), (lambda lam: theta(lam, 2)), theta1_inv, eps


def xxz_bethe_ground_state(n, J, Delta):
    """XXZ 模型（周期境界、n は偶数、J>0、Δ>-1）の基底エネルギーを Bethe 仮説で求める。

    基底状態は |1> が m = n/2 個の部分空間にあり、量子数 I_j = -(m-1)/2, ..., (m-1)/2 の
    実数のラピディティで表される。E0 = JΔn + Σ_j ε(λ_j)。
    戻り値: (E0, ラピディティ λ_j の配列, Bethe 方程式の残差の最大値)
    """
    if n % 2:
        raise ValueError("n must be even")
    m = n // 2
    theta1, theta2, theta1_inv, eps = xxz_bethe_functions(J, Delta)
    I = bethe.ground_state_quantum_numbers(m)
    if Delta <= 1:
        # 熱力学極限のラピディティの分布 ρ(λ) = 1/(2γ cosh(πλ/γ)) から作った初期値（Δ=1 は γ→π の形）
        scale = np.arccos(Delta) / np.pi if Delta < 1 else 1 / np.pi
        x0 = scale * np.arcsinh(np.tan(np.pi * I / n))
        lam, res = bethe.solve_bethe_equations(n, I, theta1, theta2, x0)
    else:
        x0 = theta1_inv(2 * np.pi * I / n)
        lam, res = bethe.solve_bethe_equations(n, I, theta1, theta2, x0, theta1_inv=theta1_inv, n_iter=5000)
    return J * Delta * n + np.sum(eps(lam)), lam, res


def xxz_e0_thermodynamic(J, Delta):
    """XXZ 模型（J>0）の熱力学極限での1サイトあたりの基底エネルギー e0 = lim E0/n。

    -1 < Δ < 1 (Δ = cos γ): e0 = J[Δ - 2 sin γ ∫ dx sinh((π-γ)x) / (sinh(πx) cosh(γx))]
    Δ = 1                 : e0 = J(1 - 4 ln 2)
    Δ > 1 (Δ = cosh η)    : e0 = J[Δ - 4 sinh η (1/2 + 2 Σ_{m≥1} 1/(e^{2mη} + 1))]
    """
    if Delta < 1:
        g = np.arccos(Delta)

        def f(x):  # sinh((π-γ)x) / (sinh(πx) cosh(γx)) を x>0 でオーバーフローしない形に書いたもの
            if x == 0:
                return (np.pi - g) / np.pi
            return 2 * np.exp(-2 * g * x) * (1 - np.exp(-2 * (np.pi - g) * x)) / (
                (1 - np.exp(-2 * np.pi * x)) * (1 + np.exp(-2 * g * x)))

        return J * (Delta - 4 * np.sin(g) * quad(f, 0, np.inf, limit=200)[0])
    if Delta == 1:
        return J * (1 - 4 * np.log(2))
    e = np.arccosh(Delta)
    m = np.arange(1, int(40 / e) + 2)
    return J * (Delta - 4 * np.sinh(e) * (0.5 + 2 * np.sum(1 / (np.exp(2 * m * e) + 1))))
