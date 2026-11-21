# 02 XXZ模型（XXZ chain）

$$H = J\sum_i \left(X_iX_{i+1} + Y_iY_{i+1} + \Delta\,Z_iZ_{i+1}\right)$$

## どんな模型か

隣り合うスピンの間に、$xy$面内の交換相互作用（第1項・第2項）と、強さが$\Delta$倍の$z$方向の相互作用（第3項）が働く模型である。$\Delta$は**異方性パラメータ**で、$\Delta$の値によって次のよく知られた模型を含む。

- $\Delta=0$：**XX模型**。Jordan-Wigner変換で自由フェルミオンになり、簡単に解ける。
- $\Delta=1$：**ハイゼンベルク模型**（XXX模型）。スピン回転の対称性 SU(2) をもつ。
- $|\Delta|\to\infty$：**イジング模型**に近づく。

**特徴**

- **全磁化の保存**：$\sum_i Z_i$ がハミルトニアンと交換する（U(1) 対称性）。$H$ は磁化ごとのブロックに分かれる。
- **相図**（$J>0$）：$\Delta>1$ ではギャップのある反強磁性（Néel）相、$-1<\Delta\le1$ ではギャップのない臨界相（朝永-Luttinger液体）、$\Delta<-1$ では強磁性相になる。$\Delta=1$ の転移は Berezinskii-Kosterlitz-Thouless（BKT）型である。
- **Jordan-Wigner変換では解けない**：変換すると、$\Delta Z_iZ_{i+1}$ が隣り合うフェルミオンの密度どうしの相互作用（4次の項）になる。$\Delta\neq0$ では自由フェルミオンにならない。
- **Bethe仮説で厳密に解ける**：相互作用があるのに可積分で、Bethe仮説で厳密に解ける（H. Bethe, Z. Phys. 71, 205 (1931)、C. N. Yang and C. P. Yang, Phys. Rev. 150, 321 (1966)）。
- **実験との関係**：擬1次元の磁性体（KCuF₃ など）や、冷却原子・超伝導量子ビットを使った量子シミュレータで実現されている。Néel状態やドメインウォール状態からの時間発展は、スピン輸送の研究の代表的な題材になっている。

横磁場イジングモデル（[01_tfim](../01_tfim/)）は自由フェルミオンに帰着する「解ける模型」の基準だったが、XXZ模型は**相互作用があるのに解ける模型**の基準である。Jordan-Wigner変換が使えないかわりに Bethe仮説を使う点が、横磁場イジングモデルとの一番の違いである。

## ノートブック

**模型の定義**

| ノートブック | 内容 |
|---|---|
| [00_model.ipynb](00_model.ipynb) | 今回扱う模型の定義、前提、対称性 |

**解き方**（各手法をこの模型に適用したもの。手法の一般論は [methods/](../../methods/) の同じ番号のノートブック）

| ノートブック | 内容 |
|---|---|
| [01_exact_diag.ipynb](01_exact_diag.ipynb) | 厳密対角化：$n=2$ の行列と手計算、全磁化によるブロック対角化 |
| [02_jordan_wigner.ipynb](02_jordan_wigner.ipynb) | Jordan-Wigner変換とBethe仮説：4次の項が残る理由、$\Delta=0$ の自由フェルミオン、Bethe方程式の導出、$n=2$ の手計算 |
| [03_trotter.ipynb](03_trotter.ipynb) | Trotter分解：偶奇ボンドへの分け方と回路、誤差の特徴、$n=2$ の手計算 |
| [04_vqe.ipynb](04_vqe.ipynb) | VQE：一重項から始めるアンザッツ（HVA）、測定、$n=2$ の手計算 |

Bethe仮説の一般論は [methods/06_bethe_ansatz.ipynb](../../methods/06_bethe_ansatz.ipynb) にある。

**物理の問い**（複数の手法の比較）

| ノートブック | 内容 |
|---|---|
| [11_ground_state.ipynb](11_ground_state.ipynb) | 基底エネルギー $E_0$ を厳密対角化・Bethe仮説・VQE で比較、熱力学極限と有限サイズ補正 |
| [12_staggered_magnetization.ipynb](12_staggered_magnetization.ipynb) | Néel状態からのクエンチ後のスタッガード磁化 $m_s(t)$ を厳密対角化・Jordan-Wigner変換（$\Delta=0$）・Trotter分解で比較 |

## 手法の比較

| | 厳密対角化 | Jordan-Wigner変換 | Bethe仮説 | 量子回路（Trotter分解・VQE） |
|---|---|---|---|---|
| 扱う行列・方程式 | $2^n\times2^n$ 行列（磁化のブロックでは $\binom{n}{n/2}$ 次元） | $2n\times2n$ 行列 | $n/2$ 変数の連立方程式 | $n$ 量子ビットの回路 |
| 使える範囲 | すべて（小さい $n$ のみ） | $\Delta=0$ のみ | 固有エネルギー（周期境界条件） | すべて |
| 精度 | 厳密 | 厳密 | 厳密 | Trotter誤差・VQEの変分誤差・測定の統計誤差がある |

XXZ模型では、基底エネルギーは Bethe仮説で大きな $n$ まで厳密に求まるが、$\Delta\neq0$ の時間発展を大きな $n$ で厳密に計算する簡単な方法はない。古典計算機では厳密対角化（またはテンソルネットワーク）に頼ることになり、量子回路による方法の意義がより大きくなる。
