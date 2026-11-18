# Quantum_many-body_system

1次元量子多体系を、**厳密対角化・自由フェルミオン（Jordan-Wigner変換）・量子回路（Trotter分解 / VQE）** などの複数の手法で解き、互いに比較しながら学ぶリポジトリです。

*Learning repository for 1D quantum many-body physics: each model is solved with several methods (exact diagonalization, free fermions via Jordan-Wigner, and quantum circuits with Trotterization / VQE), and the results are cross-checked against each other.*

## 構成の考え方

- **`methods/`：手法の解説**。一般論と、横磁場イジングモデルの $n=2$ の手計算・コードによる確認。手法の説明はここに1回だけ書く。
- **`models/`：模型ごとの「物理の問い」**。1つの問い（基底エネルギー、クエンチ後の時間発展など）につき1つのノートブックで、複数の手法を比較する。

## 手法

| # | ノートブック | 内容 | 状態 |
|---|---|---|---|
| 01 | [厳密対角化](methods/01_exact_diag.ipynb) | $2^n\times2^n$ 行列の対角化 | 完成 |
| 02 | [Jordan-Wigner変換](methods/02_jordan_wigner.ipynb) | 自由フェルミオンへの帰着、BdG 行列 | 完成 |
| 03 | [Trotter分解](methods/03_trotter.ipynb) | 量子回路による時間発展 | 完成 |
| 04 | [VQE](methods/04_vqe.ipynb) | 量子回路による基底状態の変分計算 | 完成 |
| 05 | [TEBD](methods/05_tebd.ipynb) | 行列積状態による時間発展 | 準備中 |

## 模型

| # | 模型 | 状態 |
|---|---|---|
| 01 | [横磁場イジングモデル](models/01_tfim/) | 基底エネルギー・クエンチ完成、相転移・もつれは準備中 |
| 02 | [キタエフ鎖・XY模型](models/02_kitaev_chain/) | 予定 |
| 03 | [横磁場イジング＋縦磁場](models/03_tfim_longitudinal/) | 予定 |

## フォルダ構成

```
├── methods/        手法の解説（一般論 + n=2 の手計算と確認）
├── models/         模型ごとの物理の問い（複数の手法の比較）
│   └── 01_tfim/
│       ├── README.md               TFIM の一般的な説明
│       ├── 00_model.ipynb          今回扱う模型の定義
│       ├── 01_ground_state.ipynb   E0
│       ├── 02_quench_dynamics.ipynb  m_x(t)
│       ├── 03_phase_transition.ipynb
│       └── 04_entanglement.ipynb
├── src/qmb/        共通パッケージ（ノートブックから import して使う）
│   ├── operators.py      Pauli 演算子、site_op
│   ├── models.py         ハミルトニアン・観測量の行列表示
│   ├── ed.py             厳密対角化、時間発展
│   ├── free_fermion.py   Jordan-Wigner 変換による解（BdG 行列）
│   ├── circuits.py       量子回路（Trotter 分解・VQE）
│   └── observables.py    ギャップ、相関関数、もつれエントロピー
├── tests/          手法どうしの一致を確かめるテスト（pytest）
└── docs/           記号と規約（conventions.md）
```

## セットアップ

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"   # qmb を編集可能モードでインストール
python -m pytest                    # テストを実行
```

インストール後は、どのフォルダのノートブックからでも次のように使えます。

```python
from qmb.models import tfim_hamiltonian
from qmb.free_fermion import ground_energy
```

## ロードマップ

1. 横磁場イジングモデル：基底エネルギー・クエンチ（完成）
2. TFIM の量子相転移（ギャップ、もつれエントロピー、有限サイズスケーリング）
3. XY模型・キタエフ鎖（マヨラナ端状態）
4. TFIM＋縦磁場（非可積分系、熱平衡化）
5. XXZ / ハイゼンベルク模型とテンソルネットワーク（MPS, TEBD）
6. 1次元ハバード模型（フェルミオン→量子ビット変換と VQE）
7. ノイズとエラー緩和
