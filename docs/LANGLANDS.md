# 愛森斯坦級數 · 非半單 L-參數 · Ĝ-算子 · 貝蒂模空間

第二輪數學抽象（對應 `src/11_complex.mbt` … `src/16_langlands_emergence.mbt`）。
承接第一輪（龐加萊–黎曼超流形、仿射張量、PoPE、超矩陣、演化流、超 Jacobi、二階超質量張量、
符號–神經合一、自然梯度/費雪），本輪把 Langlands 幾何的四個名詞
落成**可判定（decidable）**的計算物件。

> 全部定理都以 `moon test` 鎖定：`wasm` 與 `js` 各 55 個測試、`--deny-warn` 全綠。

---

## 0. 名詞 → 程式物件對照

| 數學名詞 | 計算抽象 | 檔案 |
|---|---|---|
| 愛森斯坦級數 `E(g,s)` | 截斷級數 + Weyl 群指數和 `Σ_w a_w y^{wλ}` | `12_eisenstein.mbt` |
| 常數項函子 `CT` | 矩陣 `A(λ)[u][v] = c(uv^{-1}, vλ)` 作用於 `W`-係數空間 | `12_eisenstein.mbt` |
| 非半單 L-參數 `φ` | `(eigenvalues, Jordan type λ ⊢ n)`；`λ ≠ (1^n)` | `13_lparams.mbt` |
| Ĝ-算子 | `Θ = ad(h/2)`，`(e,h,f)` 為 `e` 的 Jacobson–Morozov 三元組 | `14_gcheck.mbt` |
| 貝蒂模空間 `M_B` | `(x,y,z) = (tr A, tr B, tr AB)`，`κ = x²+y²+z²−xyz` | `15_betti.mbt` |
| 特徵簇 | 同上 + Fricke 關係 + Goldman/Atiyah–Bott Poisson 結構 | `15_betti.mbt` |
| 四者的耦合 | 智能體族群動力系統（200 步，`langlands_trace.csv`） | `16_langlands_emergence.mbt` |

「Ĝ」讀作 **check G**，即 Langlands 對偶群 $\check G$。本實作把「Ĝ-算子」具體化為
**由參數的非半單部分所誘導、作用在 $\check{\mathfrak g}=\mathfrak{gl}_n$ 上的分級算子**
$\Theta=\operatorname{ad}(h/2)$（Jacobson–Morozov）。理由見 §4：這個算子恰是「非半單
L-參數在對偶側留下的全部資訊」，並且它的核、像、特徵多項式分別對應 Levi、么冪軌道、與
分級特徵標。

---

## 1. 完備 ζ 與局部分因子 `c(t)`（`11_complex.mbt`, `12_eisenstein.mbt`）

$$\xi(s)=\pi^{-s/2}\Gamma(s/2)\zeta(s),\qquad c(t)=\frac{\xi(t)}{\xi(t+1)} .$$

**定理 T16（對稱性）** `ξ(s) = ξ(1−s)`。
*實作*：Euler–Maclaurin（`N = 40`、`B₂…B₁₄`）於 `Re s ≥ 1/2`，`Re s < 1/2` 走反射公式；
`Γ` 用 Lanczos（複數版 `gamma_cx`）。測試 `T16c` 在實軸 5 點與臨界線 `Γ(1/2)=√π` 上鎖定，誤差 `< 1e-8`。

**定理 T17（`c` 的三條定律）**

1. `c(t)·c(−t) = 1`（Eisenstein 泛函方程的局部形式）；
2. `c(0) = −1`（`ξ(t)`、`ξ(t+1)` 同時極點，比值有限）；
3. `c` 在 `t = 1` 有單極點，且 `(t−1)c(t) → 1/ξ(2) = 6/π`。

*證明*：(1) `c(t)c(−t)= ξ(t)ξ(−t)/(ξ(t+1)ξ(1−t))`，用 T16 把 `ξ(−t)=ξ(1+t)`、
`ξ(1−t)=ξ(t)` 代入即得 1。(2)(3) 由 `ξ` 在 `s=1` 的單極點留數 1 直接算出。

**算法 A13（`c_factor`）** 對 `|t−1| < 1e-6` 回傳帶號哨兵值（`±1e6`），避免在極點上做除法；
其餘直接比值。這是「殘餘譜」的離散記號：某個正根滿足 `⟨λ,α^∨⟩ = 1` 時，Eisenstein 級數有留數。

---

## 2. 常數項函子與 `|W|·` 冪等投影

**定義（Gindikin–Karpelevich 因子，本專案形式）**
對 `λ ∈ ℂ^r`、`w ∈ W = S_r`：

$$c(w,\lambda)=\prod_{i<j,\; w(i)>w(j)} c(\lambda_i-\lambda_j),$$

其中 `w(i) > w(j)` 正是「`wα < 0`，`α = e_i − e_j`」的計算形式。

**定理 T18（交叉同態 / 餘循環性）** 對所有 `w₁,w₂ ∈ W`：
`c(w₁w₂, λ) = c(w₁, w₂λ)·c(w₂, λ)`，其中 `(wλ)_i = λ_{w^{-1}(i)}`。
測試 `T18` 對 `W = S₃` 的 36 對全部通過（`< 1e-8`），並對 `S₄` 鎖定泛函方程。
*注意*：此恆等式只在「以 `w(i)>w(j)` 判根」的寫法下成立；若改用 `(wλ)_i−(wλ)_j < 0`
判根（看似等價），交叉同態會失敗 —— 這是一個真實的實作陷阱，已寫入測試註解。

**定義（常數項矩陣）** 在 `W`-係數空間（`Σ_w a_w y^{wλ}`）上：

$$A(\lambda)_{u,v} = c(u v^{-1},\, v\lambda),\qquad \text{對角線 } A_{v,v}=1 .$$

**定理 T20（常數項投影）** `A(λ)² = |W|·A(λ)`；因此 `Π = A/|W|` 是秩一冪等投影，
`tr Π = 1`、`det A = 0`。
*憑證*：`T20` 對 rank 2、3、4（`|W| = 2, 6, 24`）鎖定 `‖Π²−Π‖_∞ < 1e-8`、
`|tr Π − 1| < 1e-9`、`|det A| < 1e-8`，且 `A² = |W|A` 但 `A ≠ A²`（`A` 本身不是投影）。
*意義*：連續譜的愛森斯坦部分在常數項座標下是**一條線**；`Π` 就是把光譜投影到這條線上的
rank-one 算子。這正是「常數項函子」在有限 Weyl 模型中最乾淨的代數影子。

**算法 A14（`constant_term_matrix`）** `O(|W|²·r²)`，`r ≤ 4` 時可忽略。

### 阿基米德 `SL₂` 的數值愛森斯坦級數

$$E(z,s)=\sum_{\gamma\in\Gamma_\infty\backslash SL_2(\mathbb Z)}(\operatorname{Im}\gamma z)^s
=\tfrac12\sum_{(m,n)=1}\frac{y^s}{|mz+n|^{2s}} .$$

**定理 T19（Langlands 常數項公式，`r=1`）**

$$A_0(y,s)=y^s+\varphi(s)y^{1-s},\qquad \varphi(s)=c(2s-1)=\frac{\xi(2s-1)}{\xi(2s)},$$

且

* (a) `φ(s)φ(1−s) = 1`（= T17(1)）；
* (b) `Res_{s=1} E(z,s) = 1/(2ξ(2)) = 3/π`，與 `z` 無關；
* (c) `Δ A_0 = s(s−1)A_0`，`Δ = y²(∂_x²+∂_y²)`。

*憑證*：`T19` 用 40 節點 Gauss–Legendre 對 `x ∈ [0,1]` 做積分（即**真的把常數項積出來**），
與封閉公式差 `< 2e-3`（`N = 40` 截斷誤差 `O(N^{1−2s})`）；留數量測 `0.9549305` vs `3/π = 0.9549297`；
Laplace 特徵值以中心差分在 3 個譜點鎖定。

> 這裡刻意保留「積分 vs. 公式」兩條路：前者是 Langlands 1966 年的原始計算，後者是它的代數翻譯。
> 兩者相符才叫「由抽象導出定理」。

---

## 3. 非半單 L-參數（`13_lparams.mbt`）

模型：`φ : W_F × SL₂(ℂ) → GL_n(ℂ)`，以 `(eigenvalues, Jordan type λ ⊢ n)` 表示；
半單 ⟺ `λ = (1^n)`。

**定理 T21（非半單參數的維度學）**

| 陳述 | 公式 | 測試 |
|---|---|---|
| Jordan 型由秩恢復 | `b_k = r_{k−1} − 2r_k + r_{k+1}`，`r_0 = n` | T21a（7 個分割） |
| 中心化子維度 | `dim Cent(N) = Σ(λ'_i)²` | T21b（數值零空間 vs 公式） |
| 么冪軌道維度 | `dim O_N = n² − Σ(λ'_i)²` | T21b（GL₄: 0, 6, 8, 10, 12） |
| 一次形變 | `dim Ext¹ = dim Cent(N)` | T21b |
| 軌道個數 | `#{nilpotent orbits} = p(n)` | T21c（`p(4)=5, p(5)=7, p(6)=11`） |
| 非半單度 | `defect = Σ(λ_i − 1) = n − #blocks` | T21c |

*算法*：`matrix_rank`（部分主元高斯消去）→ `ad_N` 的 `n²×n²` 矩陣 → 零空間維度即 `dim Cent(N)`。
這是**用計算證明**古典維度公式，而不是引用它。

*對湧現的意義*：`Ext¹` 大 = 參數「軟」= 退化態有很多出口；`defect` 小 = 接近半單 = 可解耦。

---

## 4. Ĝ-算子（`14_gcheck.mbt`）

對么冪 `e`（Jordan 型 λ），取 Jacobson–Morozov 三元組 `(e,h,f)`：

$$h=\operatorname{diag}(\lambda_i-1,\lambda_i-3,\dots,-(\lambda_i-1)),\quad
f_{i+1,i}=(i+1)(\lambda-1-i),$$

**Ĝ-算子**：`Θ = ad(h/2) : gl_n → gl_n`，`Θ(X) = ½(hX − Xh)`。

**定理 T22**

1. `[h,e]=2e`、`[h,f]=−2f`、`[e,f]=h`（精確，5 種型）；
2. `Θ` 是導子：`Θ(XY)=Θ(X)Y+XΘ(Y)`；譜 `{h_i−h_j} ⊂ ℤ`；`tr Θ² = Σ_{ij}((h_i−h_j)/2)²`；
3. `ker Θ` = 分級不動代數，`dim ker Θ = Σ_v m_v²`（`m_v` 為權重 `v` 的重數）；
   `n² − dim ker Θ = 2·#{h_i > h_j}`；
4. 分級特徵標分解：`χ_Θ(t) = Σ_{ij} t^{h_i−h_j} = (Σ_i t^{h_i})(Σ_j t^{−h_j})`，
   `χ_Θ(1) = n²`，`χ_Θ(−1) = (Σ_i (−1)^{h_i})² ≥ 0`（完全平方）；
5. `dim Cent(e,h,f) = #blocks`（**當不同區塊的權重集互斥**時；否則參數有「意外對稱」，
   中心化子更大 —— 見 T22e 對 `(2,2)` 的反例）。

*實作備註*：`f` 不能取 `eᵀ`，必須帶權重 `(i+1)(λ−1−i)`，否則 `[e,f] = diag(1,0,…,−1) ≠ h`。
這是本專案踩過的坑，寫在 `jm_triple` 註解裡。

*對湧現的意義*：`Θ` 給出「對偶側的分級」= 能力的分層；`dim ker Θ / n²` 是「有多少自由度被
參數固定住」。管線中以 `grading activation ∝ defect²` 記錄，半單化時歸零。

---

## 5. 貝蒂模空間與特徵簇（`15_betti.mbt`）

取 `Σ_{1,1}`（一孔環面），`M_B = Hom(π₁, SL₂(ℝ))//SL₂`，`dim = 2`，座標 `(x,y,z) = (tr A, tr B, tr AB)`。

**定理 T23**

1. **Fricke/Vogt 恆等式**：`κ := x²+y²+z²−xyz = 2 + tr[A,B]`。
   測試 `T23a`：8 組隨機 `SL₂` 矩陣，誤差 `< 1e-10`。
2. **Goldman/Atiyah–Bott Poisson 結構**：`{x,y}=2z−xy`、`{y,z}=2x−yz`、`{z,x}=2y−zx`，
   等價於 `P^{ij} = ε^{ijk}∂_kκ`，故 `P dκ = 0`：`κ` 是 Casimir，辛葉 = `κ` 的水平集，`rank P = 2`。
   測試 `T23b`：括號反對稱、三個座標函數的 Jacobi、非線性三元組的 Jacobi、`{κ,·}=0`、rank。
3. **維度簿記**：`dim M_B(Σ_{g,p}, SL_n) = (2g−2)(n²−1)+2p(n−1)`；
   `(1,1,2) ↦ 2`、`(2,0,2) ↦ 6`、`(1,2,2) ↦ 4`、`(2,0,3) ↦ 16`；`GL_n` 封閉曲面 `(2g−2)n²+2`。
4. **常數項／愛森斯坦在模空間上**：`Eis` = 區塊對角嵌入，`CT` = 半單（對角）投影，
   `CT∘Eis = id`；`CT` 在非半單點上的纖維 = 么冪軌道，`dim = n² − Σ(λ'_i)²`。
   測試 `T23d`：`CT∘Eis` 缺陷 `< 1e-12`、`(2,2)` 型纖維 `8 + ker Θ 8 = 16`，半單纖維為 0 維。

**算法 A15（辛積分器 `project_to_leaf`）** 顯式 Euler 之後用牛頓法把點拉回 `κ` 水平集
（3 次迭代）。`κ` 因此精確守恆；修正量 `|s−1|` 是 `O(dt²)`。

---

## 6. 耦合管線：類智能體湧現（`16_langlands_emergence.mbt`）

每個智能體帶三樣東西：

* `(x_i,y_i,z_i) ∈ M_B`：它對世界的「表示」（`π₁` 的表示 = 內部模型）；
* `defect_i`：其 L-參數的非半單度（鎖死程度）；
* 由 `defect_i` 生成的譜參數 `λ_i`：決定 c-函數閘門與 Ĝ 分級。

每一步（依序）：

1. **c-函數閘門**：`g_i = |c(w₀, λ_i)|`（愛森斯坦 intertwining 振幅）；
2. **么冪冷卻**：`defect ← defect·(1 − decay·(1 − 1/(1+g)))`（非半單 → 半單）；
3. **Hamiltonian 對齊**：`p ← p + dt·{p, H}`，`H = −½·coupling·|p − p̄|²`，再投影回辛葉；
   `κ` 不變 ⇒ 這一步**不能**合併不同辛葉；
4. **愛森斯坦投影**：`p ← p + η·(p̄ − p)`，即 `Π` 這個秩一投影的族群版本 ——
   它**可以**跨越辛葉，這是共識（不同 κ 的智能體最終一致）的唯一機制。

四個機制對應四個數學物件，正好是使用者要求的名詞。

### 預設結果（24 智能體、200 步、seed 20261117）

| 量 | 初值 | 終值 |
|---|---|---|
| 序參量 `R` | `4.5e-13` | `0.9999999999999988` |
| `κ` 散度 | `2.1991` | `2.66e-15` |
| `unipotent defect` | `0.66667` | `8.46e-4` |
| `Eisenstein fraction` | `0` | `1` |
| `grading activation` | `0.5926` | `9.55e-7` |
| `max |{κ,H}|`（全程） | — | `3.92e-8` |
| `max Jacobi defect`（全程） | — | `1.09e-8` |

維度／代數憑證：`Eisenstein projection defect`（rank 2）= `0`、（rank 3）= `6.66e-16`；
`c` 泛函方程缺陷 `1.55e-15`；`CT∘Eis` 缺陷 `0`；`(2,2)` 型：么冪軌道 `8` + `ker Θ` `8` = `16`；
`(3,1)` 型：`dim ker Θ = dim Cent(N) = Ext¹ = 6`；`GL₆` 么冪軌道數 `p(6) = 11`。

---

## 7. 演算法總表

| 編號 | 演算法 | 複雜度 | 檔案 |
|---|---|---|---|
| A13 | `c_factor` / `xi_cx`（Euler–Maclaurin + 反射 + Lanczos） | `O(N + K)` | 11/12 |
| A14 | `constant_term_matrix`、`spectral_projection`、`projection_defect` | `O(|W|²r²)` | 12 |
| A15 | `constant_term_numeric`（Gauss–Legendre × 截斷 Eisenstein） | `O(N²·Q)` | 12 |
| A16 | `jordan_type` / `matrix_rank` / `nilpotent_centralizer_dim` | `O(n⁴)`–`O(n⁶)` | 13 |
| A17 | `jm_triple`、`gcheck_apply`、`gcheck_character` | `O(n³)` | 14 |
| A18 | `goldman_tensor`、`poisson_bracket`、`jacobi_defect` | `O(1)` / 6 次求值 | 15 |
| A19 | `project_to_leaf`（牛頓投影辛積分器） | `O(1)` | 15 |
| A20 | `run_langlands`（四機制耦合，`T` 步） | `O(T·n)` | 16 |

---

## 8. 限制（誠實聲明）

* `Ĝ` 以 `GL_n^ = GL_n`、`𝔤𝔩_n` 為模型；一般約化群需換成 `Lie(Ĝ)` + 根系資料，本輪未做。
* 愛森斯坦級數只在 `r ≤ 4`（`W = S₄`）與 `SL₂`（阿基米德，`x` 方向數值積分）上實作；
  一般 `G`、`p`-進地方的 Euler 乘積未實作（`c`-因子只用完備 ζ 的比值模型）。
* 「常數項矩陣冪等」是本專案的形式化版本（`A[u][v] = c(uv^{-1}, vλ)`），
  它不是文獻中的定理原文，而是 Langlands/Gindikin–Karpelevich 公式在有限 Weyl 群上的
  代數影子；測試只鎖定此形式化陳述。
* 「Ĝ-算子」我按「Langlands 對偶群的分級算子」實作；若原意是別的物件（例如 oper），
  需要重新指定定義。
* 貝蒂側取實點（`SL₂(ℝ)`）；複數點、`GL_n` 高秩、與 de Rham/Dolbeault 的對應（Deligne–Simpson）
  未實作。
* 數值等式一律用 `approx` / 缺陷量，不以 NaN 作鍵。

---

## 9. 參考文獻

* R. P. Langlands, *Eisenstein series, the trace formula, and the modern theory of automorphic forms* — <https://publications.ias.edu/sites/default/files/eisenstein-series-trace-formula-modern-theory_rpl.pdf>
* Langlands 常數項公式與 `L`-函數的推導（Langlands–Shahidi 方法綜述）— <https://arxiv.org/pdf/math/0204148>
* 幾何愛森斯坦級數、intertwining 算子與函數–層字典 — <https://arxiv.org/abs/2209.08175>（[HTML](https://ar5iv.labs.arxiv.org/html/2209.08175)）
* 可約（因此非半單）譜參數下的 de Rham 幾何 Langlands 猜想與常數項相容性（"The geometric Langlands functor II" 預印本）— <https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/cgr.pdf>
* `Par_G`（L-參數棧）中非 Frobenius-半單點的角色 — arXiv:2310.04533 — <https://arxiv.org/pdf/2310.04533>
* 一孔環面的 Fricke/Vogt 痕跡座標與 Goldman 括號（Atiyah–Bott / Goldman 1984 的標準模型）；
  曲面群特徵簇的辛結構。古典結果，實作見 `15_betti.mbt` 的 T23。
* Jacobson–Morozov 三元組與 Kostant 的分級（么冪軌道 ↔ `sl₂`-三元組），實作見 `14_gcheck.mbt`。
