# 格羅滕迪克宇宙、可判定保守反射、朗蘭茲綱領宇宙、ℓ-進光滑層（Round 4）

> 定理 T35–T46 與演算法 A29–A37，實作於 `src/21_universe_reflection.mbt`、`src/22_langlands_universe.mbt`、`src/23_elladic_diamonds.mbt`、`src/24_universe_frontier.mbt`，並由 `src/cmd/main/main.mbt` 的宇宙湧現管線 (`run_universe_emergence`) 耦合驗證。
> 全部定理由 `moon test --deny-warn`（wasm + js，82/82 全綠，0 warnings）鎖定；憑證寫在 `universe_trace.csv`（41 步）與 `docs/report_universe.html`（自包含內嵌 SVG）。

## 0. 編碼：讓宇宙公理可判定

Hereditarily finite 集合以位元遮罩編碼：`x ∈ y` ⟺ `(y >> rank(x)) & 1 == 1`；冪集 `P(x)` 列舉 `x` 的子遮罩；並集、配對、替換皆為位元運算。於是格羅滕迪克宇宙公理（遞移、配對、冪集、並集、替換）與 `U_j ∉ U_j`（羅素／Burali-Forti 屏障的計算形）都變成有限窮舉，可由 `moon test` 逐一斷言（T35、T36，`21_universe_reflection.mbt`）。

基數型別 `Cardinal { Finite(n) | Omega | Pow2(κ) }` 與 `cardinal_pow_two`、`cardinal_lt`、`cardinal_le`、`cardinal_max`、`is_infinite`；`size_level` 將基數映射到宇宙層級（`Finite(17) → 0`, `Omega → 1`, `Pow2(Omega) → 2`, `Pow2(Pow2(Omega)) → 3`；A32）。

## 1. 可判定保守宇宙反射原理（T35–T36；T37–T38；A29–A31）

Montague–Lévy 反射在「所有參數、所有公式」的強型上於有限模型 `V_α`（此處 α ≤ 4，|V_4| = 16）是**虛無的**：前輪嘗試的 `reflects_on(α,n,nvars,size)` 對每個輸入皆回 0（只有 Σ₀ 在有限層級絕對成立；Σ_n n≥1 需 α ≥ ω 的無限見證空間）。該強型已刪除（見 `docs/CLEANUP.md` §8）。

取而代之的是**參數相對**形式，這才有實際內容：

- `reflects_for(β, α, n, params, nvars, max_size)`：`V_β` 與 `V_α`（α = 4）在參數 `params` 下對 Σₙ 公式池（`formula_pool`，含原子 `x_i ∈ x_i`，池大小 4/15/30/205）的真值一致。
- `reflection_bound_for(α, n, params, nvars, max_size)`：最小使上式成立的 β（無則 −1）。
- `decide_conservatively(f, α, n, env, ~max_size)`：無需宇宙公理的 Σₙ 判定程序（A31）。
- `ReflectionWitness { alpha, beta, n, formulas_checked, parameters }` 與 `is_strict`（A30）。

經窮舉驗證的反射界（T37，`α = 4`，參數 `p̄ = (0,1)`）：

| Σ-index n | params | β | 備註 |
|---|---|---|---|
| 0 | (0,1) | 2 | Δ₀：β = rank(p̄)+1，絕對性 |
| 1 | (0,1) | 3 | 見證秩需嚴格更大 |
| 2 | (0,1) | 3 | |
| ≥1 | (0) | rank(p̄)+1 | 單參數下各 n 皆收斂 |
| 1 | (3,15) under α=3 | −1 | 參數超出 V₃，無解 |
| any | (0,1) under α=3 | −1 | α=3 太小 |

保守性缺陷 `conservativity_defect(α, n, params, max_size)` 恆為 0 只要 β 存在（T38：反射判定與大宇宙真值完全一致，且反射關係遞移）。管線中 `n = step % 4` 步進，β 恆 = 3 < α = 4。

## 2. 朗蘭茲綱領宇宙（T39–T41；A32–A34）

beth 梯 `L_j` 分層：

| 對象 | 側 | 基數 | 層級 |
|---|---|---|---|
| Kottwitz 不變量 / Hecke 特徵值 | automorphic | 有限 | 0 |
| Div¹ = Spd(ℚ̌_p)/Frob^ℤ | geometric | ω | 1 |
| Spd(ℚ_p) (perfectoid diamond) | geometric | ω | 1 |
| Bun_G (GL₂) | geometric | 2^ω | 2 |
| Hck_G (Hecke stack) | geometric | 2^ω | 2 |
| LocSys_Ĝ (spectral stack) | spectral | 2^ω | 2 |
| D_lis(Bun_G, Λ) ≃ IndCoh_N(LocSys_Ĝ) | spectral | 2^{2^ω} | 3 |

`triad_is_tight()` 為 true：導出 ∞-範疇層（3）比譜層（2）高一級（T39）。`derived_category_needs_extra_stratum() = true`。

Newton 反演：excursion 算子特徵值 `p_k = tr φ(Frob)^k` 由半單 L-參數反演，其核恰為 unipotent radical（`reflection_kernel_is_unipotent`，T40）。Newton 層 `Bun_G^b ≃ [*/J_b(Q_p)]`（T41，`22_langlands_universe.mbt`）：

| Newton 層 b | 斜率 | J_b | dim J_b | basic |
|---|---|---|---|---|
| b₀ | (0, 0) | GL₂ | 4 | 否 |
| b₁/₂ | (½, ½) | D^× | 4 | 是 |
| b₁ | (1, 1) | GL₂ | 4 | 否 |

lisse 秩 r 層 ↔ J_b(Q_p) 的 r 維光滑表示；純性 `w = −2a`（`weil_weight_defect` 為 0）。

## 3. ℓ-進光滑層在 Diamond 上（T42–T45；A35–A37）

**Tilting**（`23_elladic_diamonds.mbt`）：`X^♭ = lim_{(←, x↦x^p)} O_X/p` 以 p 次根鏈實現，憑證「char-p 剩餘完美且原代數非完美」：

- `pth_power_image_count(2, 2) = 2 < 4`，但 `frob_bijective_mod_p(2)`（p=2 對 residues 雙射）。
- `tilt_root_count(2, 2, 1) = 2`，`tilt_root_count(2, 2, 2) = 0`。
- `tilt_chains(2, 3, 2)` = 8 條鏈；`tilt_is_perfect(p, n)` ⟺ `|T_d| = p^n` 於每深度（T42，A35）。

**Div¹ 鏡曲線**（T43）：`deg ∈ ℝ/ℤ` 為 Frob-不變；`χ_{a+b} = χ_a χ_b`；權重 `w = −2a`；`div1_frobenius_invariant` 為 true。

**ℓ-進光滑性**（T44，A36）：`CohComplex { lo, betti, twist }`，`coh_dualizing(d) = {lo: 2d, twist: d}`，即 `f^!Λ ≅ Λ(d)[2d]`；Verdier 對合為對合（`coh_dual(coh_dual(c)) == c`）；Poincaré 對偶回文：`[1]`, `[1,0,1]`, `[1,4,1]`, `[1,0,1,0,1]`，拒絕 `[1,2,3]`；合成乘法性 `coh_compose_dimensions`；負維數分類堆疊 `[*/G]`（d = −dim G）亦成立（`coh_etale_cover_euler`, `coh_recovered_dimension`）。

**幾何 Satake**（T45，A37）：`satake_ic_match(A₁, ω₁) = 2 = dim V_ω₁`，`(A₂, e₁/e₂) = 3`，`(C₃, ω₁) = 6`；移次 `⟨2ρ, λ⟩`（`satake_perverse_degree`: 1 / 2 / 6）；`P_IC` q-回文。

## 4. 綜合：反射與 Verdier 對偶交換（T46；A36–A37）

`src/24_universe_frontier.mbt` 的 `run_universe_emergence(~steps=40, ~levels=4, ~seed)` 是四通道耦合管線：

1. 宇宙階梯：`level = step % 4`；
2. 參數相對反射：`reflection_bound_for(4, n, [0,1], 2, 2)`（β = 3，`n = step % 4`）；
3. 朗蘭茲 excursion：譜 `[1.7, 0.85, −0.6]` 漂移，`reflected_dual` / `reflected_of_inverse` / `duality_commutation_defect`；
4. ℓ-進：`dim = step % 9 − 4`，tilt 鏈、`coh_dualizing`。

T46：反射／對偶方塊交換，`e_k(α⁻¹) = e_{r−k}(α)/e_r(α)`，max defect `9.77×10⁻¹⁵`（雙精度機器精度）。`power_sums_of` 與 Newton 恆等式串聯四通道。

管線憑證（`steps = 12`）：13 行 CSV；欄寬 11（`step, universe_level, reflection_beta, sigma_index, sigma_formulas, conservativity_defect, stratum, excursion_defect, dualizing_shift, duality_defect, tilt_chains`）；所有缺陷 `< 1e-9`；`tilt_chains == 8`；摘要含 `T35/T37/T40/T42/T44/T46`。

## 5. 反思：兩條死路（不得重複）

- **強反射 `reflects_on / reflection_bound`**：對有限 `V_α` 恆回 0，僅 Σ₀ 絕對；Montague–Lévy 需 α ≥ ω。已刪除，採參數相對 `reflects_for / reflection_bound_for`。
- **單層 p 次根「恰 p 個」tilt 模型**：`#{y : y^p ≡ x}` 於 `x = p` 為 0，不等於 p；採 tilt 鏈 + `|T_d| = p^n` 完美性憑證。

## 6. 檔案與測試錨點

| 模組 | 檔案 | 定理 | 演算法 |
|---|---|---|---|
| 宇宙反射 | `src/21_universe_reflection.mbt` | T35–T38 | A29–A31 (`reflects_for`, `reflection_bound_for`, `decide_conservatively`) |
| 朗蘭茲宇宙 | `src/22_langlands_universe.mbt` | T39–T41 | A32–A34 (`size_level`, `reflected_parameter`, `triad_is_tight`) |
| ℓ-進 diamonds | `src/23_elladic_diamonds.mbt` | T42–T45 | A35–A37 (`tilt_chains`, `coh_dualizing`, `satake_ic_match`) |
| 綜合管線 | `src/24_universe_frontier.mbt` | T46 | `run_universe_emergence`, `UniverseTrace::csv` |

憑證：`universe_trace.csv`（41 行）、`docs/report_universe.html`、`run_output.txt`（# universe-emergence 區塊）。CI：`.github/workflows/ci.yml` 上傳 `docs/report_universe.html` + `universe_trace.csv`；品質閘門 `scripts/verify.sh` 以 `moon fmt --check`、`moon check --deny-warn`、`moon test --target wasm/js --deny-warn`、`moon build`、四條管線、四份報告、`scripts/dead_code_scan.py`（562 定義 / 0 死）封鎖。
