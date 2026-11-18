"""qmb: 1次元量子多体系シミュレーションの共通コード。

- operators    : Pauli 演算子、site_op
- models       : ハミルトニアン・観測量の行列表示
- ed           : 厳密対角化
- free_fermion : Jordan-Wigner 変換による自由フェルミオン解
- circuits     : 量子回路（Trotter 分解・VQE、Qiskit が必要）
- observables  : ギャップ、相関関数、もつれエントロピー（模型によらない）
"""
__version__ = "0.1.0"
