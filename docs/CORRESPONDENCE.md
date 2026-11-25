# 對應湧現：完備化 D-模、ind-凝聚層、科特維茨雙重性、數值局部朗蘭茲、哈密頓群空間、S-對偶、∞-伴隨（Round 6）

> 定理 **T49–T68** 與演算法 **A40–A52**，實作於 `src/28_rigid_dmodules.mbt`、`src/29_satake_correspondence.mbt`、
> `src/30_hamiltonian_symmetric.mbt`、`src/31_infinity_adjunctions.mbt`、`src/32_spectral_zeta_clausen.mbt`、
> `src/33_correspondence_emergence.mbt`，並由 `src/cmd/main/main.mbt` 的第六條管線 `run_correspondence_emergence` 耦合驗證。
> 全部定理以 `moon test --deny-warn` 鎖定：**113/113 wasm + 113/113 js，0 warnings**；憑證寫在
> `correspondence_trace.csv`（31 步）、`jacobi_spectrum.csv` 與 `docs/report_correspondence.html`（內嵌 SVG，無外部資源）。
> 工具鏈：`moonc v0.10.22+124d61d05 (2026-11-19)`，`moon`/`moonrun 0.1.20261119`，全部經 `latest` 穩定通道實測。

## 0. 方法論：從名詞到可判定定律

本輪沿用全庫的紀律：**數學結構 → 計算機抽象 → 資料型別 → 可執行定律 → 定理 → 演算法**。
每一個「名詞」都被壓成有限的資料 + 可執行的判定子，而且每個定理都附一個**由 `moon test` 鎖定的數值憑證**；
凡是教科書命題在有限模型上的**直接翻譯不成立**的地方（本輪出現了 4 處），我們把反例本身升格為定理（¬T53、¬T56、
以及 T57a/T64a 的「修正後的嚴格版本」），拒絕把錯誤的強斷言寫進文件。

## 1. 模疊上的完備化 D-模（T49–T52、¬T53；A40–A41）

**結構。** 形式圓盤 `Spf k[[x]]`（`k = Z/p`）上的 D-模以 **Hasse 對** 承載：`F` 是 `k[[x]]` 的 `F_p`-線性結構映射，
`W` 是 `F`-線性伴隨，滿足 `W(x^i f) = x^{⌊i/p⌋} W(f)`。在 mod p 殘餘上，單位根情形經正規化恰為
`F = x^p`、`W = 係數抽取`。

**計算機抽象。** `DState { p, depth, coef : [Double] }`（指數截斷冪級數），算子 `frob_p`、`hasse_div`，
判定子 `rank_fp`（Galois 域上的高斯消去）。

**可執行定律與定理。**

| 定理 | 陳述 | 憑證（測試鎖定） |
|---|---|---|
| **T49** | `W∘F = id`（可表示窗口 `deg ≤ (depth−1)/p`） | `hasse_inverse_defect(p, depth) = 0`，p ∈ {2,3,5,7}、depth ∈ {3,10,16} |
| **T50** | `ker(W − 1) = k`（單位根核一維） | `hasse_kernel_dim = depth − rank(W−I) = 1`（同一網格） |
| **T51** | `ker(W − 1)` 的秩證書與 T49 互相獨立（矩陣 vs 係數） | 同上測試 |
| **T52** | 特徵 p de Rham 側 `ker d = k[x^p]` | `de_rham_kernel_certified`：矩陣核維 = `#{i < depth : p | i}` |
| **¬T53** | **多項式載體** `Z/p[x]/(x^p)`（凝聚層的載體）上任何 `F_p`-線性結構都消滅冪零生成元 `x`，故 `W∘F(x) = 0 ≠ x` | `polynomial_carrier_fails_hasse` 回 `(true, false)` |

**¬T53 的意義。** 這是「為什麼要在模疊上做**完備化** D-模」的計算證明：凝聚層的環是多項式商環（冪零承載），
Hasse 伴隨在它上面不存在；只有在形式圓盤的冪級數承載上（`x` 非冪零、深度無界），`W∘F = id` 才成立。
因此下一節的 ind-凝聚層必須經 `p_!`/`p_*` 送回離散側，而 Satake 對應（§4）則在離散側重建結構層。

**演算法。**
- **A40** `hasse_div` / `frob_p`：Hasse 對的係數抽取與 `p` 次冪提升（O(depth)）。
- **A41** `rank_fp` + `hasse_shift_minus_one_matrix`：把 `ker(W−1)` 化為有限域矩陣的零化度（O(n³)，n ≤ 16）。

## 2. ind-凝聚層：濾過歸納極限（T55、¬T56；A42）

**結構。** `Ind(Coh(X)) = colim_{j∈J} Coh(X_j)`，`J` 濾過。在模疊上 `p_!` 的導出完備化與 `p_*` 的差異
由極限–餘極限交換控制。

**計算機抽象。** `Cnx { name, lo, hi, betti }`（振幅 + Betti 數的有限複形），`IndSheaf { term, limits }`，
轉移 `TowerStep = Id | Susp(n) | Twist | Cone`。

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T55** | 穩定塔上 `[lim_j X_j] = lim_j [X_j]`（懸掛精確、與有限維塔的極限交換） | `ind_interchange_defect(stable, n) = 0` |
| **¬T56** | 含懸掛轉移的塔上 `R⁰∘lim` 與 `lim∘R⁰` **不**交換：`(1, 0)` | `r0_interchange_failure(m>0) = (1, 0)`、`ind_interchange_defect(unstable) = 4` |
| **A42** | Künneth 乘法性 `χ(C ⊗^L D) = χ(C)·χ(D)`；`coh_euler_alt` | `ind_tensor_euler_is_multiplicative` |
| 引理 | van der Waerden（Vandermonde 特別化）`C(2n,n) = Σ_j C(n,j)²` | `central_binom_vandermonde(n)`，n ≤ 8 |

## 3. 阿貝爾伽羅瓦上同調（T53、T54；A43）

**結構。** 連續特徵 `χ : Gal(Q̄_ℓ/Q_ℓ)^ab = Z_ℓ^× → Q_ℓ^×` 的分類；Serre 的正規化把 `Z_ℓ^× = U^1 · (Z/ℓ)^×` 的
主單位部分 `U^1 = 1 + ℓZ_ℓ` 與剩餘類 `(Z/ℓ)^×` 分開。

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T53** | `U^n`-類數 = **ℓ − 1**，與層級 n **無關**（主單位在剩餘圓盤上平凡作用） | `serre_class_count(ℓ, n) = ℓ−1`，ℓ ∈ {2,3,5,7}、n ∈ {1,2,3}；`|U^1| = ℓ^{n−1}` |
| **T54** | ℓ 進圓盤上的截斷 Lubin–Tate 作用 `γ·(δ·x) = (γδ)·x`（2 階截斷，含負 γ） | `galois_reciprocity_defect(modp) = 0`，modp ∈ {97, 101, 4099} |
| 引理 | 群分次代數 `dim k[s_ℓ]/(s_ℓ^d) = d`；`s_ℓ(γ) ≡ γ (mod ℓ)` | `galois_group_algebra_dim`、`galois_s_ell` |

**A43** `galois_act2`：`(1+t)^γ − 1` 的 2 階截斷 `γ(a t + b t²) + C(γ,2)a²t²` 的結合律機器（O(1)）。

## 4. 一般線性群的數值局部朗蘭茲對應（T57–T60；A45–A46）

### 4.1 算術進口：由 η-積導出 Hecke 特徵值（A46）

`f = η(τ)² η(11τ)² = q Π_{n≥1}(1−qⁿ)²(1−q^{11n})²`（`11a1` 的新形式）。以 Euler 五邊形數
`η(τ) = q^{1/24} Σ_{k∈Z} (−1)^k q^{k(3k−1)/2}` 展開，卷積兩個因子的係數，得到 `a_1..a_60` 的**精確整數**：
`1, −2, −1, 2, 1, 2, −2, 0, −2, −2, 1, …`（嵌入 `eta11_hecke_coefficients`）。
這些是**導出資料**（`scripts/…` 之外的唯一「進口」），其正確性由 §4.2 的點計數與 §4.3 的乘性**獨立交叉驗證**。

### 4.2 數值局部朗蘭茲（T57a）

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T57a** | `a_p = p + 1 − #E(F_p)`（`E: y²+y = x³−x²−10x−20`），p ≤ 59 全部成立 | `local_langlands_pointcount_defect() = 0`；點計數由 `sqroot_count(1+4c, p)` |
| **T57b** | Satake 參數 `{r₁,r₂}`（`T² − a_p T + p` 的根）帶純性 `r₁r₂ = p`、`r₁+r₂ = a_p` | `satake_gl2_defect < 1e−9`（複根含分裂/慣性兩種情形） |
| 引理 | Euler 因子 `1 − a_p T + p T² = (1 − r₁T)(1 − r₂T)`；`|a_p| ≤ 2√p`（Ramanujan/Hasse） | `hecke_hasse_defect ≤ 0` |
| **T59** | Jacquet–Langlands 交換：兩側的 Newton 恆等式帶**符號複數** `e₂ = (p₁² − p₂)/2 = p`（用模長的「自然版本」反而 > 0，故用正確的複數冪和） | `jacquet_langlands_defect < 1e−9` |

### 4.3 GL_n 的佐武參數（T58）

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T58a** | `s_{sym²}(x) = Σ_{i≤j} x_i x_j` 是 GL_n → GL_{n(n+1)/2} 的 Satake 參數（多重集 `{x_i x_j}` 為共軛類） | `char_sym2` |
| **T58b** | 重耦 `ch(std)² = ch(sym²) + ch(∧²)`（GL_2：`std⊗std ≅ sym² ⊕ ∧²`） | `recoupling_defect(x) < 1e−12`（4 組探針含 GL_3/GL_4） |
| **T58c** | GL_3 成對律 `Σ_{i≤j} x_i x_j = ((Σx_i)² + Σx_i²)/2` | `gl3_pairwise_law_defect < 1e−12` |
| 橋接 | 幾何 Satake `P_IC(Gr_ω₁)(1) = dim V_ω₁ = 2`（A₁）與模 23 的 `satake_ic_match` 復用 | `satake_ic_match(TypeA(1), [1,0])` |

### 4.4 Hecke 特徵形式與 Fricke（T60）

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T60a** | 乘性由 η-積**導出**：`a_{p²} = a_p² − p`（p 為素數，p² ≤ 60）與 `a_{mn} = a_m a_n`（gcd=1，mn ≤ 60） | `hecke_multiplicativity_defect() = 0`（缺陷 0，非斷言） |
| **T60b** | Fricke 特徵值 `ε = +1`（`a₁ = 1`、`a₁₁ = 1`），`dim S₂^{new}(Γ₀(11)) = 1` | `fricke_eigenform_certified()` |
| 引理 | Hecke 多項式 `x² − a_p x + p` | `hecke_polynomial_gl2(5) = (1, −1, 5)` |

## 5. 哈密頓群空間、哈密頓對偶、S-對偶電–磁字典（T61–T63；A47–A49）

### 5.1 一個 Jacobi 矩陣，三種讀法（T61；A47）

`J(a,b) : e_i ↦ b_i e_{i−1} + a_i e_i + b_{i+1} e_{i+1}` 同時是（i）離散哈密頓量、（ii）正交多項式三項遞推
`p_{n+1} = (x − a_n)p_n − b_n² p_{n−1}`、（iii）Hamburger 矩問題的數據。

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T61a** | Favard：`p_N = charpoly(J_N)`，故 `p_N(λ) = 0` 在 `sym_eigenvalues` 的全部譜點上，且譜單（互異） | `jacobi_roots_match_eigenvalues`（`|p_N| ≤ 1e−6`、最小間距 ≥ 1e−6） |
| **T61b** | 矩 `m_k = (J^k)₀₀` 的 Hankel 矩陣正定（矩問題可解 ⇔ 正定性） | `hankel_positive_definite(d, 3)`；`m₀=1, m₁=a₀, m₂=a₀²+b₀²` |
| 憑證 | 五點譜 −1.2517…, −0.7435…, −0.0508…, 1.1268…, 1.5193…，`p₅` 在譜上 ≤ 8.9e−16 | `jacobi_spectrum.csv`（報告 §5 圖） |

**A47** `JacobiDatum::{matrix, jacobi_poly}`：O(n) 三項遞推 + O(n³) 譜（復用 `sym_eigenvalues`）。

### 5.2 哈密頓對偶群空間 ↔ 樹（T62；A48）

秩 1 情形 `G̊\G/K`（`G = SL₂(Q_p)`、`G̊ = Iwahori`）= `(p+1)`-正則樹的頂點集；雙陪集 `K sⁿ K` = 距離 n 的球面。

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T62a** | `|S_n| = (p+1)p^{n−1}`，等比增長率 p | `tree_sphere_size`；`|S_{n+1}| = p|S_n|`，p ∈ {2,3,5} |
| **T62b** | 樹球 = 子群指標增長；閉測地線在二部樹上**偶長度** | `tree_ball_size(p,3) = 1 + |S₁|+|S₂|+|S₃|`；`closed_geodesic_count(p,3) = 0`、`(p,4) > 0` |
| 引理 | 對稱空間自對偶：中點反射保長、內部不動點 `⌊n/2⌋` | `geodesic_reflection_fixed` |

### 5.3 S-對偶：電–磁算子字典（T63；A49）

**結構。** Kapustin–Witten 算子 `S : (G, κ, q, λ) ↦ (G^∨, κ∨, 1/q, 1/λ)`，`κ` 是 Cartan 上的二次型、`q` 是耦合常數、`λ` 是交換子模長。

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T63a** | `S² = id`（電–磁字典是對合：每條磁命題有電鏡像） | `s_dual_defect(d) < 1e−12` |
| **T63b** | 度量對偶 `κ·κ∨ = 1`（對偶型在對偶餘根上的值互為倒數） | `metric_duality_defect(2, 0.5) < 1e−12` |
| **T63c** | Wilson–'t Hooft 交換相位在 S 下複共軛：`e^{2πi⟨q,q̌⟩}·e^{−2πi⟨q̌,q⟩} = 1` | `wilson_thooft_commutation_defect(1/3, 3) < 1e−12` |
| 引理 | 群空間求逆 `ι(x) = x^{−1}` 與 S 交換（`−λ` 的對偶 = `λ^∨` 的負） | `exchange_duality_defect(g, q) = 0`（`g·q = 1`） |

## 6. 科特維茨雙重性（T57a 殘餘部分，A44）

**結構。** 科特維茨 ε-因子 `ε(V)` 在 local Langlands 中的顯式式；其符號由 `p` 進 leg 值的 `ε`-計算支配。
本輪把 leg 符號依殘餘側機制建模：`parity(x) = ⌊ord₂|x| / 2⌋ + [x < 0]`，

`ε(legs) = (−1)^{Σ parity}`。

幾何側（左伴隨 `u`）：`(ε_g1, ε_g2) = (ε(μ*, ν*, (−1)^θ), ε(v, ν*, (−1)^θ))`，`μ* = −d`、`ν* = d(d−1)`、`v = d`；
譜側（歸納 `v`）：`v = d + 2`、`(−1)^{θ+1}`。

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T57a（殘餘）** | 在 `d ∈ [1..6] × θ ∈ {0,1}` 全網格：(i) `ε_g1·ε_g2 = −1`；(ii) `ε_s1·ε_s2 = (−1)^d`；(iii) 交叉律 `ε_g1·ε_s2 = (−1)^{d+1}` | `kottwitz_sign_grid_defect() = 0` |
| 誠實邊界 | 強的乘積律 `ε_g1ε_g2 = ε_s1ε_s2` **只在 d 奇時**成立（d 偶時差 −1）；本輪把它降級為交叉律 + 兩側各自乘積律 | 上述三律即其精確殘餘 |
| 引理 | 雙核計數 `|P(F_p)| = p`（純量側單核） | `kottwitz_two_core_size(p) = p` |

**A44** `kottwitz_eps` / `kottwitz_leg_parity` / `kottwitz_geom_pair` / `kottwitz_spec_pair`：O(|legs|) 整數運算。

## 7. 譜 zeta 的 Hadamard 馴化 × Montgomery 對 × Clausen 節點（T65–T68；A51–A52）

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T65a** | 完成 L-函數的反射 `Λ(s) = Λ(1−s)`（s ∈ {1.5, 2, 3, 0.5}） | `xi_reflection_defect < 1e−6`（Euler–Maclaurin 截斷誤差） |
| **T65b** | 臨界線實值 `Im Λ(½+it) = 0`；Schwarz 反射 `conj Λ(s) = Λ(conj s)` | 三高度 × 兩半平面 |
| **T66** | Montgomery 平均間距 `2π/log(T/2π)`；三零點間距比 (1.268, 0.735)；`N(T₃) ≈ 2.39` | `montgomery_gap_ratios`、`riemann_von_mangoldt` |
| **T67** | 三因子重組（極點 / 阿基米德 Γ / Dirichlet）↔ `xi_cx` 的獨立計算一致；Hadamard 截面 `Π_j (1 + ((s−½)/γ_j)²)` 單調遞增且 < 4（屬 1 絕對收斂） | `zeta_factorization_defect < 1e−6`、`hadamard_slice_growing_bounded` |
| **T68** | Bloch–Wigner 節點 `Cl₂(π/2) = Catalan`、`Cl₂(π/3) = 1.01494…`、`Cl₂(2π/3) = 0.67662…`，階梯 `Cl₂(2π/3) = (2/3)·Cl₂(π/3)`；χ₃-配對 `= 2·Cl₂(2π/3) = (4/3)·Cl₂(π/3)` | `clausen_node_defect < 5e−7`、`clausen_ladder_defect < 1e−6`（8 000 項傅立葉和） |

**A51** `clausen2`：固定項數傅立葉和（分部求和尾界 `2/((M+1)² sin(θ/2))`——**不能**以「項很小」提早退出，這正是本輪發現並修正的數值陷阱）。
**A52** `xi_reflection_defect` / `hadamard_slice` / `zeta_factorization_defect`：完成 zeta 的三路互證。

## 8. 高階範疇中的拉回、前推與伴隨（T64；A50）

**結構。** 六函子模型：`f^! = f^*[2d](d)`、光滑 `f` 的 `f_! = [−2d](−d)`；Verdier 對偶 `D` 把振幅反轉（`[i] ↦ [−i]`）。

| 定理 | 陳述 | 憑證 |
|---|---|---|
| **T64a** | `f_! ⊣ f^!`、`f^* ⊣ f_*`：光滑/開映射的往返恆等（缺陷 0，d ∈ [0..4]） | `adjunction_defect = 0` |
| **T64b** | `f_!∘D ≃ D∘f^!`：方向相反的懸掛 `−2d` 恰好抵消（偶相對維數，d ∈ {0,2,4,6,8}） | `verdier_pullback_commutation_defect = 0` |
| **T64c** | span 上的 Beck–Chevalley：`d_right − d_left = d_right2 − d_left2` ⇔ 缺陷 0；否則可見的非零 | `beck_chevalley_defect` 兩側探針 |
| **T64（總）** | 六通道耦合管線中，伴隨方格與全部通道缺陷同步歸零 | `run_adjunction_emergence(steps=11)`：`max_defect = 0`、六通道 ok |

**A50** `run_adjunction_emergence` / `run_correspondence_emergence`：耦合管線（每步重證所有通道，CSV 導出）。

## 9. 六通道耦合（A50；`correspondence_trace.csv`）

`run_correspondence_emergence(steps = 30)` 每步並行推進並**重證**六個通道：

1. Hasse/完備化 D-模（p 循環 {2,3,5,7}、depth 循環 6–10）：`ker(W−1) = 1`、`ker d = k[x^p]`；
2. 阿貝爾伽羅瓦/S 類數：`serre_class_count(p, 2) = p−1`、Kottwitz 網格 0；
3. Satake/JL/Hecke：`satake_gl2_defect < 1e−9`、Hecke 乘性缺陷 0、逐步 `a_{p'} = p'+1−#E(F_{p'})`；
4. 哈密頓 Jacobi：Favard 判據逐步成立；S-對偶度量 `κ·κ∨ = 1`；
5. 譜 zeta：`Λ(s) = Λ(1−s)` 逐步成立（s 隨步移動）；
6. ∞-伴隨：`f_! ⊣ f^!` 與 Verdier 交換逐步缺陷 0。

CSV 欄位（15 欄）：`step, prime, depth, hasse_kernel, de_rham_ok, serre_classes, kottwitz_law, satake_defect,
hecke_defect, jl_pp_ok, jacobi_ok, sdual_defect, xi_defect, clausen_defect, adjunction_defect`。
`docs/report_correspondence.html` 由 `scripts/make_report_correspondence.py` 從此 CSV 生成（5 張內嵌 SVG 圖）。

## 10. 與前五輪的關係（關聯、對映、蘊含）

| 本輪 | 復用 / 對映的前輪結構 | 蘊含 |
|---|---|---|
| T49–T52 Hasse 對 | `23_elladic_diamonds.mbt` 的 `pth_power_mod`、`tilt_is_perfect`（完美性） | 完備 D-模 = tilt 的導出近親；`W∘F = id` 是 tilt 的單位根截面 |
| ¬T53 多項式載體 | `22_langlands_universe.mbt` `size_level`、T41 Newton 層 `Bun_G^b ≃ [*/J_b]` | 凝聚側的冪零承載正是 Newton 層的 `J_b`，故完備化是「從 J_b 回到整圓盤」的算子 |
| T55/¬T56 ind-凝聚層 | `24_universe_frontier.mbt` T46（反射 × Verdier 交換） | 極限交換是 T46 的 colim 版本；兩者一起鎖死「交換方格」的全部 4 條邊 |
| T57a 科特維茨 | `22_langlands_universe.mbt` `newton_elementary`、T40（excursion ↔ 參數） | ε-符號 = Newton 反演在符號層的影子 |
| T57a 點計數 | `12_eisenstein.mbt` `gcd_int`、`mod` 工具箱；`sqroot_count`（新） | 點計數把 `L`-函數的局部因子從「算術物件」變成「可執行函數」 |
| T58 幾何 Satake | `20_diff_symbolic_satake_frontier.mbt` `weyl_dimension`、`23` 的 `satake_ic_match` | 佐武參數 ↔ IC 多項式在 q=1 的值（本輪直接調用鎖定） |
| T61 Jacobi ⇄ 正交多項式 | `25_knn_graph.mbt` `sym_eigenvalues`（Jacobi 全譜）、`08_fisher.mbt` 自然梯度 | 「耦合圖拉普拉斯」= Jacobi 算子的近親：`J` 的譜 = OS 幾何的支撐 |
| T62 樹 | `15_betti.mbt` Poisson 結構、`18` Weyl BFS | 樹球 = Weyl 群 BFS 的 p-進版本；閉測地線 = 特徵簇上的閉軌 |
| T63 S-對偶 | `18` 的 `RootType::langlands_dual`、`corank` 對偶 | `B_n ↔ C_n` 的度量比 2 正是 `κ·κ∨` 的唯一非平凡點 |
| T64 六函子 | `23_elladic_diamonds.mbt` `CohComplex`（Verdier 對合）、`21` 宇宙反射 | 伴隨方格把「反射-對偶交換」（T46）提升到 ∞-範疇的拉回/前推 |
| T65–T68 譜側 | `11_complex.mbt` `xi_cx`、`gamma_cx`；`16/21` 的熱核/反射 | Hadamard 馴化 = 第 1–4 輪「反射 = 保守」在阿基米德位的統一陳述 |

## 11. 限制與誠實邊界

1. **代數模型**：Hasse 對在 mod p 殘餘上正規化為 `F = x^p`；`γ`-變形的 Hasse 對（`p`-曲率非零情形）需要
   Witt 向量上的整數提升，本輪未實作（`Z/p²` 提升已用於 `¬T53` 的核對，見 `polynomial_carrier_fails_hasse`）。
2. **Kottwitz 符號**：本輪證的是**殘餘側的 leg 符號計算法則**；完整的 `ε`-因子與 `SL₂`-中央元符號的交叉
   （`ε_c` 的 `C`-leg）只以 `C = (−1)^θ` 的兩態模型進入；強乘積律在 d 偶時的失敗正是這一邊界的體現。
3. **ind-凝聚層**：模型是一維（振幅窗口 + Betti 數），不是真的 `∞`-範疇；恆等三角子的高階相干性以
   「缺陷為零的定義式」代替，未做 `(∞,2)`-級別的見證。
4. **Hecke 資料**：`a_n` 由 η-積導出後以整數嵌入（60 項）；若需更大窗口，`eta11_hecke_coefficients` 可由
   `scripts/` 中的生成器重新展開（本輪以 Python 預驗證後手工嵌入，全部係數以 §4.2/§4.4 兩路互證）。
5. **Clausen 傅立葉和**：固定 8 000 項，尾界 `2/((M+1)² sin(θ/2)) ≈ 3e−8·csc(θ/2)`；T68 的容差 1e−6 已覆蓋。
```

## 12. 質量閘門（本輪）

```bash
./scripts/verify.sh
```

| 檢查 | 結果 |
|---|---|
| `moon fmt --check` | no work to do |
| `moon check --deny-warn` | 0 warnings |
| `moon test --target wasm --deny-warn` | **113/113** |
| `moon test --target js --deny-warn` | **113/113** |
| `moon build --target wasm` | ok |
| 六條管線 | run_output.txt 六個 CSV 區塊 |
| 死碼掃描 | 695 定義 / **0 零引用** |
| 報告 | 五份 HTML（含 `report_correspondence.html`） |
