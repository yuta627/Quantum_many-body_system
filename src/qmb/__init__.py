"""qmb: 1次元量子多体系シミュレーションの共通コード。

- operators    : Pauli 演算子、site_op、|+...+> や Néel 状態
- models       : ハミルトニアン・観測量の行列表示
- ed           : 厳密対角化、保存量によるブロック対角化
- free_fermion : Jordan-Wigner 変換による自由フェルミオン解
- bethe        : Bethe 方程式（対数形）の数値解法
- circuits     : 量子回路（Trotter 分解・VQE、Qiskit が必要）
- observables  : ギャップ、相関関数、もつれエントロピー（模型によらない）
"""
__version__ = "0.1.0"
