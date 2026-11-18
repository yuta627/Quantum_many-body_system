# Quantum_many-body_system

1次元量子多体系を、**厳密対角化・自由フェルミオン（Jordan-Wigner変換）・量子回路（Trotter分解 / VQE）** などの複数の手法で解き、互いに比較しながら学ぶリポジトリです。

*Learning repository for 1D quantum many-body physics: each model is solved with several methods (exact diagonalization, free fermions via Jordan-Wigner, and quantum circuits with Trotterization / VQE), and the results are cross-checked against each other.*

## 模型×手法の対応表

| 模型 | 厳密対角化 | Jordan-Wigner変換 | Trotter分解 | VQE | TEBD |
|---|:-:|:-:|:-:|:-:|:-:|
| [01 横磁場イジングモデル](models/01_tfim/) | [✓](models/01_tfim/01_exact_diag.ipynb) | [✓](models/01_tfim/02_jordan_wigner.ipynb) | [✓](models/01_tfim/03_trotter.ipynb) | [✓](models/01_tfim/04_vqe.ipynb) | |
| [02 キタエフ鎖・XY模型](models/02_kitaev_chain/) | | | | | |
| [03 横磁場イジング＋縦磁場](models/03_tfim_longitudinal/) | | [✗](models/03_tfim_longitudinal/README.md#jordan-wigner変換が使えない理由) | | | |

✓：実装済み（リンク先はその手法を模型に適用したノートブック）　✗：原理的に使えない（リンク先に理由）　空欄：未実装

## 構成の考え方

- **`methods/`：手法の一般論**。模型によらない部分だけを書く。
- **`models/<模型>/`：模型ごとの内容**。
  - `00_model`：模型の定義
  - `01`〜`04`（**解き方**）：各手法をその模型に適用したもの（`methods/` と同じ番号）。模型ごとの導出と $n=2$ の手計算、コードによる確認。
  - `11`〜（**物理の問い**）：1つの問い（基底エネルギー、クエンチ後の時間発展など）につき1つのノートブックで、複数の手法を比較する。

## 手法

| # | ノートブック | 内容 | 状態 |
|---|---|---|---|
| 01 | [厳密対角化](methods/01_exact_diag.ipynb) | 行列表示、対角化と時間発展、計算コスト | 完成 |
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
├── methods/        手法の一般論（模型によらない）
│   ├── 01_exact_diag  02_jordan_wigner  03_trotter  04_vqe  05_tebd
├── models/         模型ごとの内容
│   └── 01_tfim/
│       ├── README.md                 TFIM の一般的な説明
│       ├── 00_model.ipynb            今回扱う模型の定義
│       ├── 01_exact_diag.ipynb       ┐
│       ├── 02_jordan_wigner.ipynb    │ 解き方（methods/ と同じ番号）
│       ├── 03_trotter.ipynb          │
│       ├── 04_vqe.ipynb              ┘
│       ├── 11_ground_state.ipynb     ┐
│       ├── 12_quench_dynamics.ipynb  │ 物理の問い
│       ├── 13_phase_transition.ipynb │
│       └── 14_entanglement.ipynb     ┘
├── src/qmb/        共通パッケージ（ノートブックから import して使う）
│   ├── operators.py      Pauli 演算子、site_op
│   ├── models.py         模型ごとの定義（行列表示、JW 変換後の A, B など）
│   ├── ed.py             厳密対角化、時間発展
│   ├── free_fermion.py   自由フェルミオン（BdG 行列）の一般的な計算
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
from qmb.free_fermion import ground_energy, bdg_matrix
```

## ロードマップ

1. 横磁場イジングモデル：基底エネルギー・クエンチ（完成）
2. TFIM の量子相転移（ギャップ、もつれエントロピー、有限サイズスケーリング）
3. XY模型・キタエフ鎖（マヨラナ端状態）
4. TFIM＋縦磁場（非可積分系、熱平衡化）
5. XXZ / ハイゼンベルク模型とテンソルネットワーク（MPS, TEBD）
6. 1次元ハバード模型（フェルミオン→量子ビット変換と VQE）
7. ノイズとエラー緩和
