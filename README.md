# Quantum_many-body_system

1次元量子多体系を、**厳密対角化・自由フェルミオン（Jordan-Wigner変換）・量子回路（Trotter分解 / VQE）** などの複数の手法で解き、互いに比較しながら学ぶリポジトリです。各模型のノートブックは「一般論 → n=2 の手計算 → コードで確認 → 大きい n」という流れで、読んで学べる形にしています。

*Learning repository for 1D quantum many-body physics: each model is solved with several methods (exact diagonalization, free fermions via Jordan-Wigner, and quantum circuits with Trotterization / VQE), and the results are cross-checked against each other.*

## 模型一覧

| # | 模型 | 手法 | 状態 |
|---|---|---|---|
| 01 | [横磁場イジングモデル](models/01_tfim/) | 厳密対角化 / Jordan-Wigner変換 / Trotter分解・VQE | 完成 |
| 02 | [キタエフ鎖・XY模型](models/02_kitaev_chain/) | Jordan-Wigner変換（マヨラナ端状態） | 予定 |
| 03 | [横磁場イジング＋縦磁場](models/03_tfim_longitudinal/) | 厳密対角化 / 量子回路（非可積分系） | 予定 |

## フォルダ構成

```
├── src/qmb/        共通パッケージ（ノートブックから import して使う）
│   ├── operators.py      Pauli 演算子、site_op
│   ├── models.py         ハミルトニアン・観測量の行列表示
│   ├── ed.py             厳密対角化、時間発展
│   ├── free_fermion.py   Jordan-Wigner 変換による解（BdG 行列）
│   └── circuits.py       量子回路（Trotter 分解・VQE）
├── models/         模型ごとのノートブックと要約
├── tests/          手法どうしの一致を確かめるテスト（pytest）
└── docs/           記号と規約（conventions.md）
```

## セットアップ

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"     # qmb を編集可能モードでインストール
pytest                       # テストを実行
```

インストール後は、どのフォルダのノートブックからでも次のように使えます。

```python
from qmb.models import tfim_hamiltonian
from qmb.free_fermion import ground_energy
```

## ロードマップ

1. 横磁場イジングモデル（完成）
2. TFIM の量子相転移（ギャップ、もつれエントロピー、有限サイズスケーリング）
3. XY模型・キタエフ鎖（マヨラナ端状態）
4. TFIM＋縦磁場（非可積分系、熱平衡化）
5. XXZ / ハイゼンベルク模型とテンソルネットワーク（MPS, TEBD）
6. 1次元ハバード模型（フェルミオン→量子ビット変換と VQE）
7. ノイズとエラー緩和
