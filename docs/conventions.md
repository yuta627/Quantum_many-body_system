# 記号と規約

リポジトリ全体で共通の記号・規約をまとめる。新しい模型を追加するときもこれに従う。

## Pauli 演算子とハミルトニアン

- $X, Y, Z$ は Pauli 行列。$Z\ket{0}=\ket{0},\ Z\ket{1}=-\ket{1}$、$X\ket{\pm}=\pm\ket{\pm}$。
- 横磁場イジングモデル: $H = -J\sum_i Z_iZ_{i+1} - h\sum_i X_i$（$J>0$ で強磁性）。
- 境界条件は特に断らない限り**開放境界条件**（`pbc=False`）。

## 量子ビットの並び順

- `qmb.operators.site_op` などの行列表示は **`np.kron` の順**: サイト0が一番左（状態ベクトルの添字の最上位ビット）。
- **Qiskit は逆順**（qubit 0 が最下位ビット）。`SparsePauliOp` のラベル文字列も右端が qubit 0。
- 左右反転で対称な量（今回の $H$ や $M_x$）の期待値は並び順によらないが、サイトごとの量を比べるときは注意する。

## Jordan-Wigner 変換

- $\ket{+}$ をフェルミオンなし、$\ket{-}$ をフェルミオンありとする（$X_j = 1 - 2c_j^\dagger c_j$）。
- $c_j = \left(\prod_{k<j} X_k\right) a_j$、$a_j = \ket{+}\bra{-}_j = (Z_j - iY_j)/2$。
- このとき $\ket{+\cdots+}$ はフェルミオンの真空。
- Nambu 表示 $\Psi = (c_1,\dots,c_n,c_1^\dagger,\dots,c_n^\dagger)^T$、$H = \frac12\Psi^\dagger H_{\mathrm{BdG}}\Psi$。

## 量子回路のゲート

- $R_{ZZ}(\theta) = e^{-i\frac{\theta}{2}Z\otimes Z}$、$R_X(\theta) = e^{-i\frac{\theta}{2}X}$（Qiskit の定義）。
- $e^{-iH_{ZZ}\delta t} = \prod_j R_{ZZ}(-2J\delta t)$、$e^{-iH_X\delta t} = \prod_j R_X(-2h\delta t)$。
- 1次の Trotter ステップは $e^{-iH_X\delta t}e^{-iH_{ZZ}\delta t}$（$ZZ$ 層を先にかける）。

## ノートブックの書き方

**`methods/`（手法の解説）** は次の型にそろえる。

1. 一般論
2. 具体例（横磁場イジングモデルの $n=2$ の手計算）
3. コードで確かめる（一般の $n$ で書いたコードに $n=2$ を代入し、手計算と一致するか確認）

**`models/`（物理の問い）** は次の型にそろえる。

1. 何を調べるか（手法の説明は `methods/` へのリンクで済ませる）
2. 複数の手法による計算と比較（大きい $n$ を含む）
3. 物理的な考察とまとめ

どちらも、関数はノートブックの中で定義せず `src/qmb/` から import する。
