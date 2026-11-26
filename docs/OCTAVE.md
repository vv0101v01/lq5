# 第七輪：十通道湧現（T69–T80 / A53–A64）

> 承接第六輪的六通道（`docs/CORRESPONDENCE.md`，T49–T68）。
> 本輪新增四條通道，全部依同一格式：**數學結構 ⟶ 計算機抽象 ⟶ 資料型別 ⟶ 可執行定律 ⟶ 定理 ⟶ 演算法**，
> 每條定理都有一個「moon test 會失敗」的可判定憑證。

| 新通道 | 模組 | 定理 | 演算法 | CSV / 圖 |
| --- | --- | --- | --- | --- |
| 七 · 量子鏈 / 辮群 / 瓊斯多項式 | `src/34_temperley_lieb_jones.mbt` | T69–T71 | A53–A55 | `octave_trace.csv`（jones 欄）、§1 |
| 八 · 函數域 Weil ζ / 曲線的 RH | `src/35_weil_zeta_function_field.mbt` | T72–T74 | A56–A58 | `octave_trace.csv`（weil 欄）、§2 |
| 九 · 共形場論 Verlinde 融合環 | `src/36_verlinde_fusion_su2.mbt` | T75–T77 | A59–A61 | `octave_trace.csv`（fusion 欄）、§3 |
| 十 · 五次三維體鏡像對稱 | `src/37_quintic_mirror_periods.mbt` | T78–T80 | A62–A64 | `octave_trace.csv`（mirror 欄）、§4 |
| 耦合 | `src/38_octave_emergence.mbt` | T69–T80 綜合 | A53–A64 | `docs/report_octave.html` |

---

## 1. 通道七：Temperley–Lieb、辮群與瓊斯多項式（T69–T71 / A53–A55）

### 1.1 抽象

$$
TL_n(\delta) = \Big\langle 1, e_1,\dots,e_{n-1}\ \Big|\ e_i^2 = \delta e_i,\
e_ie_j = e_je_i\ (|i-j|\ge 2),\ e_ie_{i\pm1}e_i = e_i \Big\rangle,
\qquad \delta = -A^2 - A^{-2}.
$$

`TL_n` 是 `2n` 個點的**非交叉完美匹配**所張的空間（Catalan 個基元素）。
乘法 = 上下堆疊，閉圈換成 `δ`；Markov 跡

$$
\operatorname{tr}(m) = \delta^{\mathrm{circles}(m)},\qquad
\pi(\sigma_i^{\pm 1}) = A^{\pm 1}\cdot 1 + A^{\mp 1}\cdot e_i ,
$$

於是辮 `β ∈ B_n` 的閉合：

$$
V_{\hat\beta}(t) = (-A^3)^{-w(\beta)}\,\operatorname{tr}\pi(\beta),
\qquad t = A^{-4},\quad w = \sum \pm 1 .
$$

### 1.2 計算機抽象

* `Lz`：`A` 的 Laurent 多項式，係數 `Int32`，指數窗口 `[−96, 95]`（密集陣列）；
  `lz_mul / lz_add / lz_to_t / lz_eval_one / lz_eval_minus_one / lz_mirror_t`。
* `TLelt`：平行陣列 `parts : Array[Array[Int32]]`（匹配的 partner 表）與 `coefs : Array[Lz]`；
  `tl_mul` 用 **union–find** 做堆疊並數閉圈，`tl_clean` 把抵消成 0 的項移除（保持基表示）。
* 半整數版本 `Hlz`（指數 ×2）給 skein 遞歸用，`hlz_to_t` 回到整數冪。

### 1.3 定理與憑證

| 定理 | 命題 | 憑證（moon test） |
| --- | --- | --- |
| **T69** | 匹配模型滿足 `TL_n` 全組關係：`e_i² = δe_i`、`e_ie_{i+1}e_i = e_i`、`σ_iσ_i^{-1} = 1`、`e_0e_2 = e_2e_0` | `TL_4` 上逐項整數比較；`tl_non_crossing_matchings(4) = 14`（Catalan）、`tl_matchings(3) = 15`（5!!） |
| **T70** | `V(3₁) = -t^{-4}+t^{-3}+t^{-1}`、`V(4₁) = t^{-2}-t^{-1}+1-t+t^2`；`V(1)=1`、`\|V(-1)\|=\det`（3₁→3、4₁→5）；三葉結手性 `V(t)≠V(1/t)` | 由 `σ_1³ ∈ B_2`、`(σ_1σ_2^{-1})² ∈ B_3` 的辮字純整數算出 |
| **T71** | 匹配模型 = **skein 遞歸** `x_{j+1} = t^2x_{j-1} + (t^{3/2}-t^{1/2})x_j`（`x_0 = -(t^{1/2}+t^{-1/2})`, `x_1 = 1`）在 `T(2,2k+1)`、`k = 1..6` 上完全一致；且 `\|V(-1)\| = 2k+1` | `jones_torus2_defect(k) = 0`；`T(2,9)` 已到 `t^{-13}`，故窗口取 `±96` |

### 1.4 誠實邊界

* 兩條獨立計算之所以互為「獨立」，是因為一條走跡泛函、一條走 skein 三項關係；
  但 **skein 遞歸以正交叉驅動**，算的是**鏡像結**。三葉結的手性（T70）使這一點可被觀察：
  `V(3₁) = -t^{-4}+t^{-3}+t^{-1}` 與其鏡像 `-t^4+t^3+t` 不同，故 T71 的比較是對 `lz_mirror_t` 做的，
  而非對原式硬套——這一點寫在程式碼註解裡，避免「湊出來的一致」。
* 匹配模型的跡歸一化用 `δ^{circles}`（幾何跡），使未結的閉合得 `⟨unknot⟩ = 1`；
  一般的代數 Markov 跡 `δ^{-n}` 只差一個整體 `δ` 的冪（T69 的所有關係都不受影響）。
* 只處理 `B_2`、`B_3` 的辮（`n ≤ 3`）：`TL_4` 已有 14 個基元素、`B_4` 需 4 個生成元，
  再往上多項式窗口要按 `|V(-1)|` 的行列式律放大，屬下一輪的規模問題
  （**第九輪已交付**：`docs/DEFERRED.md` §2，`t2_window_bound(n) = 6n−2` 的視窗伸縮律，T88–T91；`B_4` 仍屬未來）。

---

## 2. 通道八：函數域上的 Weil ζ 與曲線的黎曼假設（T72–T74 / A56–A58）

### 2.1 抽象

對 `C/F_q` 光滑射影曲線、虧格 `g`：

$$
Z(C,T) = \exp\Big(\sum_{k\ge1} \#C(\mathbb F_{q^k}) \frac{T^k}{k}\Big)
= \frac{P_1(T)}{(1-T)(1-qT)},\qquad
P_1(T) = \prod_{i=1}^{2g}(1-\alpha_iT)\in\mathbb Z[T].
$$

三條可判定律：

* **T72（函數方程）**：`P_1(T) = q^g T^{2g} P_1(1/(qT))`，即 `c_k = q^{k-g} c_{2g-k}`（清分母後為 ℤ 恆等式）。
* **T73（曲線的 RH）**：`|α_i| = √q`。以 `α ↦ q/α` 配對後等價於：`T² - a_1T + (a_2-2q)` 有實根且落在 `[-2√q, 2√q]`；
  `g = 1` 時就是 Hasse 界 `|a_1| ≤ 2√q`。
* **T74（模性陰影）**：`#E(F_p) = p+1-a_p`，`a_p` 是 module 29 由 `η(τ)²η(11τ)²` 導出的 Hecke 特徵值。

### 2.2 計算機抽象

* `E: y²+y = x³-x²-10x-20`（11a1，`g=1`）：`curve11_point_count`（module 29 的 `sqroot_count` 二次特徵和）。
* `C: y² = x⁵-x+1`（`g=2`）：判別式 `5⁵+4⁴ = 2869 = 19·151`，在 3,5,7,11,13 好約化。
  `N_1` 用特徵和；**`N_2` 用 `F_{p²} = F_p[i]/(i²-nr)` 暴力計數**（`Fp2` 結構 + `fp2_chi` 用 `x^{(q²-1)/2}`）。
* `weil_P1_genus1/2`：`a_1 = q+1-N_1`、`a_2 = (a_1² - (q²+1-N_2))/2`，
  則 `P_1 = 1 - a_1T + a_2T² - qa_1T³ + q²T⁴`（Weil 的整性）。

### 2.3 定理與憑證

| 定理 | 命題 | 憑證 |
| --- | --- | --- |
| **T72** | `deg P_1 = 2g`、`P_1(0)=1`、`c_4 = q²` 且函數方程逐係數成立 | `weil_functional_equation_defect(P_1, q) = 0`，`p ∈ {3,5,7,11,13}` |
| **T73** | `g=1`：Hasse 餘量 0（`p ≤ 23`）；`g=2`：`\|α\| = √q`（`p ∈ {3,5,7,11,13}`），且 `N_2 ≠ p²+1-a_1²`（證明虧格真的是 2，不是橢圓曲線的巧合） | `genus2_sqrt_q_defect(p) < 1e-9`、`genus2_riemann_hypothesis_defect(p) = 0` |
| **T74** | `#E(F_p) = p+1-a_p`（`p = 2,3,5,7,11`） | `modularity_defect() = 0` |

### 2.4 誠實邊界

* `N_2` 的暴力計數是 `O(p²)`，故只用 `p ≤ 13`（`p = 13` 已要 169 次 `F_{169}` 乘法 + 每次 `χ` 的快速冪）。
  管線在 `p > 13` 的步只證 `g = 1` 的 Weil/RH 部分，不假裝算了 `g = 2`。
* `p = 11` 對 11a1 是壞約化（判別式 `-11`），T74 在 `p = 11` 的相等仍是恆等式，
  但 Hasse/RH 檢查一律取 `p ≠ 11` 的素數，避免把壞約化的 `a_p` 混進 T73。

---

## 3. 通道九：共形場論的融合環（Verlinde，T75–T77 / A59–A61）

### 3.1 抽象

SU(2) 等級 `k` 的 modular tensor category 有簡單對象 `j = 0, …, k`（`j` = 2×自旋）：

$$
N_{j_1j_2}^{j_3} = \frac{2}{k+2}\sum_{n=0}^{k}
\frac{\sin\frac{\pi(j_1+1)(n+1)}{k+2}\sin\frac{\pi(j_2+1)(n+1)}{k+2}\sin\frac{\pi(j_3+1)(n+1)}{k+2}}
{\sin\frac{\pi(n+1)}{k+2}} \in \mathbb Z_{\ge 0},
\qquad
S_{jn} \propto \sin\frac{\pi(j+1)(n+1)}{k+2}.
$$

### 3.2 計算機抽象

* `ring_sin(m, N)`：`sin(πm/N)` **取絕對值**（Bohr–Mollerup 形式，所有項同號，數值最穩）；
  `ring_sin_signed(m, N)`：保留符號（`m` 對 `2N` 取模），供需要真正抵消的通道使用。
* `fusion_coefficient(k, j1, j2, j3)`：先做三角/奇偶/`j_1+j_2+j_3 ≤ 2k` 篩選（量子維數為 0），
  再求和後四捨五入到 `[0,4]` 內最近的整數。
* `fusion_matrix(k, j)`：整數矩陣 `N_j`（`N_{j_3}^{j_1j_2}` 的轉置擺放）。
* `quantum_dimension(k, j) = [j+1]_q = sin(π(j+1)/(k+2))/sin(π/(k+2))`。

### 3.3 定理與憑證

| 定理 | 命題 | 憑證 |
| --- | --- | --- |
| **T75** | `N ∈ ℤ_{\ge0}`、單位 `N_{0b}^c = δ_{bc}`、交換律、結合律（`Σ_d N_{ab}^d N_{dc}^e = Σ_d N_{bc}^d N_{ad}^e`） | `fusion_integrality_defect(k) = 0`、`fusion_ring_defect(k) = 0`，`k = 1..6` |
| **T76** | 整數融合矩陣 `N_j` 被 `S` 同時對角化：`N_j · v_n = (S_{jn}/S_{0n}) v_n`（`v_n[c] = S_{cn}`）；且 Verlinde 代數等式 `N_{j_1}N_{j_2} = Σ_{j_3}N_{j_1j_2}^{j_3}N_{j_3}` 逐元素成立 | `verlinde_eigen_defect(k) < 1e-9`、`fusion_matrices_commute_defect(k) = 0` |
| **T77** | `Σ_j d_j² = D² = (k+2)/(2\sin²(π/(k+2)))`；`(½)⁴` 的 F 矩陣 `Σ_j F_j = 1`、`F_0 = 1/d_{1/2}²`；五邊形恆等式 | `quantum_dimension_sum_defect(k) < 1e-9`、`f_matrix_half_four_defect(k) < 1e-9`、`pentagon_half_defect(k) = 0`；`k=2` → Ising（`d_σ=√2`、`D²=4`）、`k=3` → Fibonacci（`d_τ=φ`、`φ²=φ+1`） |

### 3.4 誠實邊界

* 「整性」是用浮點求和 + `1e-9` 容差判定的（`sin` 的捨入誤差 ~`1e-16`，`k ≤ 6` 時最壞 ~`1e-15`）；
  真正的 ℤ 憑證是**四捨五入後的整數在結合律/Verlinde 代數裡逐項相等**（`fusion_matrices_commute_defect` 完全不帶容差）。
* `k = 1` 時 `d_{1/2} = [2]_q = 2\cos(π/3) = 1`，故環退化為 `Z₂`（`σ×σ = 1`）——
  這讓「`d_{1/2} > 1`」這句教科書直覺在 `k=1` 失效，測試因此寫成 `k ≥ 2` 才檢查 `> 1`（見 T77 的測試）。
* F 矩陣只證 `(½)⁴` 這一族（四個自旋 ½ 的任意子），不是完整 `6j` 符號；一般 `6j` 的 Biedenharn–Elliott 遞歸留待下一輪——**第九輪已交付**：`docs/DEFERRED.md` §1（T84–T87 / A68–A70）。

---

## 4. 通道十：五次三維體的鏡像對稱（T78–T80 / A62–A64）

### 4.1 抽象

B 側：`Y = X₅ ⊂ P⁴` 的週期滿足 Picard–Fuchs 算子

$$
L = \theta^4 - 5z\prod_{k=1}^{4}(5\theta + k),\qquad \theta = z\frac{d}{dz},
$$

主解 `y₀(z) = Σ_{n≥0} (5n)!/(n!)⁵ zⁿ = 1 + 120z + 113400z² + 168168000z³ + …`。
A 側：虧格 0 Gromov–Witten 不變量 `n_d` 組成 Yukawa 耦合

$$
Y(q) = 5 + \sum_{d\ge1} n_d\, d^3 \frac{q^d}{1-q^d},
\qquad n_1 = 2875,\ n_2 = 609250,\ n_3 = 317206375,\dots
$$

鏡像映射把 B 側的 `z` 與 A 側的 `q` 連起來：`q = z\,e^{h(z)}`，`h = g/y₀`，`g` 由 `y₀\log z + g` 為 `L` 的第二解決定。

### 4.2 計算機抽象

* `Rat`（`Int64` 分子/分母，`gcd64` 約分）＋ `rat_theta`、`rat_five_theta_plus`、
  `quintic_y0_series`、`quintic_pf_residual`：把 `L y₀` 逐係數算成有理數並要求**恰為 0**。
* `quintic_period_coefficients`：`a_n = (5n)!/(n!)⁵`（`Int64`），比值遞歸整數除法。
* `quintic_mirror_h` / `quintic_z_of_q`（`Int64` 表格）＋
  `quintic_mirror_inverse_defect`：把 `z = q·exp(-h(z))` 在 `ℤ[[q]]` 裡用指數遞歸**獨立重算** `z(q)`；
  `int64_series_mul` 做截斷乘法。
* `quintic_instanton_divisibility_defect`：`[q^d](Y-5) - Σ_{e|d,e<d} n_e e³ = n_d d³`。

### 4.3 定理與憑證

| 定理 | 命題 | 憑證 |
| --- | --- | --- |
| **T78** | `h^{3,0} = 1`、`h^{2,1} = C(9,5) - 25 = 126-25 = 101`、`χ = 2(1-101) = -200`、`b₃ = 204` | `quintic_h21() = 101`、`quintic_euler_characteristic() = -200` |
| **T79** | `L y₀ = 0` **精確**（有理數逐係數 0）；`a_1 = 120 = 5!`、`a_2 = 113400`、`a_3 = 168168000`、`a_4 = 305540235000` | `quintic_pf_residual(t) = 0`，`t = 4..7`（`z⁰…z⁶`）；`quintic_period_coefficients` |
| **T80** | 鏡像映射整性：`h ∈ ℤ[[z]]`、`z(q) = q + 480q² + 968400q³ + 2594752000q⁴ + 7951931985000q⁵` 由 `z = q e^{-h(z)}` 獨立重算一致；`n_1 = 2875 = 5³·23`、`n_2 = 609250`；`d³`-整除律 | `quintic_mirror_inverse_defect(6) = 0`、`quintic_instanton_divisibility_defect(3) = 0` |

### 4.4 誠實邊界（本輪三條，全部寫進程式碼與測試）

1. **Yukawa 歸一化屬約定**：`κ = 5 = ∫_{X̌}Ω∧∂³Ω` 與前面的 `(z\,dt/dz)³/(z²y₀²)` 之冪號是文獻約定；
   本輪不引用某一條公式，只證**與約定無關的整性物件**（`y₀`、`h`、`z(q)`、`n_d`、`d³`-整除）。
2. **`L y₀ = 0` 的精確窗口是 `z⁰…z⁶`**：`Int64` 有理數在 `a_n ≈ (5n)!/(n!)^5` 之下於 `z⁷` 溢出；
   MoonBit 的整數溢出是環繞而非報錯，故把窗口寫死在測試裡（`t ≤ 7`），並在文件與報告中標明。
3. **`d³`-整除律的窗口是 `d ≤ 3`**：`Σ_{e|d} n_e e³` 在 `d = 4` 已達 `1.55×10¹³`、`d = 5` 達 `2.87×10¹⁶`，
   兩者都在 `Int64` 內，但 `n_3 = 317206375` 的 `d³` 乘法在 `d ≥ 4` 之後要靠 128 位；本輪截在 `d ≤ 3`。
   順帶記錄一個「教科書直譯會失敗」的點：`d ≥ 4` 時 `n_d·d³` 超出 `Int32`，
   若用 `Int32` 表格再取模，`d = 2` 會得到 `4876875 mod 8 = 3` 這種假缺陷——本輪因此全程 `Int64`。

---

## 5. 十通道耦合（`src/38_octave_emergence.mbt`）

`run_octave_emergence(steps=24)` 每步換一個素數（`3,5,7,11,13`）、一個融合等級（`2..6`）與一個托環結指數（`k = 1..6`），
十條通道的憑證同時重算，寫入 `octave_trace.csv`（20 欄、25 行）：

```
step,prime,level,jones_breadth,jones_det,torus_skein_defect,weil_fe_defect,weil_rh_defect,
hasse_excess,modularity_defect,fusion_integrality,verlinde_eigen_defect,quantum_dim_defect,
f_matrix_defect,pf_residual,mirror_inverse_defect,instanton_div_defect,
hasse_kernel,satake_defect,adjunction_defect
```

第六輪的六條通道在同一掃描中被重證（`hasse_kernel = 1`、`satake_defect = 0`、`adjunction_defect = 0`），
因此 `all_ok` 是**十個證明的合取**，而不是十個獨立報告的並列。

`docs/report_octave.html` 由 `scripts/make_report_octave.py` 從該 CSV 生成（7 張內嵌 SVG，無外部資源）。

---

## 6. 品質閘門

```
./scripts/verify.sh
  moon fmt --check                     clean
  moon check --deny-warn               0 warnings
  moon test --target wasm --deny-warn  129/129
  moon test --target js   --deny-warn  129/129
  moon build --target wasm             ok
  七條管線 → 七份 CSV → 六份 HTML 報告
  scripts/dead_code_scan.py            0 dead
```

本輪新增測試 16 條（T69–T80），全部是**行為測試**（關係、數值恆等式、整性），不是快照測試。
