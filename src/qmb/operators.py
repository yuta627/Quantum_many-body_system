"""Pauli 演算子と、1サイトだけに作用する演算子の構成。

量子ビットの並び順は np.kron の順（サイト0が一番左 = 最上位ビット）。
Qiskit は逆順なので注意（docs/conventions.md を参照）。
"""
from functools import reduce

import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op, i, n):
    """i番目のサイトだけに op が作用する 2^n 次元の行列 I⊗...⊗op⊗...⊗I を返す。"""
    if i < 0 or i >= n:
        raise ValueError("i must be in the range [0, n-1]")
    op_list = [I2] * n
    op_list[i] = op
    return reduce(np.kron, op_list)


def plus_state(n):
    """|+ + ... +> （全成分が 1/sqrt(2^n)）"""
    return np.ones(2**n, dtype=complex) / np.sqrt(2**n)
