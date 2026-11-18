# 01 横磁場イジングモデル（TFIM）

$$H = -J\sum_i Z_iZ_{i+1} - h\sum_i X_i$$

基底エネルギー $E_0$ と、初期状態 $\ket{+\cdots+}$ からのクエンチ後の横磁化 $m_x(t)$ を3つの手法で求め、互いに比較する。

| 手法 | 扱う行列・回路 | 結果 |
|---|---|---|
| 1. 厳密対角化 | $2^n\times2^n$ 行列 | $n=2$ の手計算と一致 |
| 2. Jordan-Wigner変換 | $2n\times2n$ の BdG 行列 | 厳密対角化と $10^{-14}$ の精度で一致、$n=100$ も計算可能 |
| 3. 量子回路（Trotter分解・VQE） | $n$ 量子ビットの回路 | Trotter 誤差は $O(\delta t^2)$、VQE は $p=n/2$ 層で厳密値に到達 |

ノートブック（この順に読む）:
1. [01_exact_diag.ipynb](01_exact_diag.ipynb) — モデルの定義、厳密対角化
2. [02_jordan_wigner.ipynb](02_jordan_wigner.ipynb) — Jordan-Wigner変換
3. [03_quantum_circuit.ipynb](03_quantum_circuit.ipynb) — 量子回路（Trotter分解・VQE）

2と3は共通パッケージ `qmb` を使うので、先にリポジトリ直下で `pip install -e .` を実行しておく。

対応する共通コードは `src/qmb/`（`models.py`, `ed.py`, `free_fermion.py`, `circuits.py`）、テストは `tests/test_tfim.py`。
