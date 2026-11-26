# 類智能體湧現的超幾何抽象 · MoonBit 0.10

[![gate](https://github.com/vv0101v01/frame1/actions/workflows/ci.yml/badge.svg)](https://github.com/vv0101v01/frame1/actions/workflows/ci.yml)
[![tests](https://img.shields.io/badge/tests-146%2F146%20wasm%20%2B%20146%2F146%20js-brightgreen)](#2-品質閘門全部通過一條指令)
[![HF Space](https://img.shields.io/badge/%F0%9F%A4%97%20Space-frame1--hexad-yellow)](https://huggingface.co/spaces/Transynergy/frame1-hexad)
[![MoonBit](https://img.shields.io/badge/MoonBit-0.10.22-1f6feb)](https://www.moonbitlang.com/download)

以 **數學結構 ⟶ 計算機抽象 ⟶ 資料型別 ⟶ 可執行定律 ⟶ 定理 ⟶ 演算法** 的順序，把下列數學名詞實作成一個可編譯、可測試、可執行的 MoonBit 專案，並把抽象過程中浮現的定理與演算法整理在 [`docs/THEORY.md`](docs/THEORY.md)（第一輪，T1–T15 / A1–A12）、[`docs/LANGLANDS.md`](docs/LANGLANDS.md)（第二輪，T16–T24 / A13–A20）、[`docs/FRONTIER.md`](docs/FRONTIER.md)（第三輪前沿擴充與限制消除，T25–T34 / A21–A28）、
[`docs/UNIVERSE.md`](docs/UNIVERSE.md)（第四輪，T35–T46 / A29–A37）、`docs/HEXAD_6IN1_PLAN.md`（第五輪，T47–T48）
、[`docs/CORRESPONDENCE.md`](docs/CORRESPONDENCE.md)（第六輪對應湧現，T49–T68 / A40–A52）
與 [`docs/OCTAVE.md`](docs/OCTAVE.md)（第七輪十通道，T69–T80 / A53–A64）、
[`docs/DEFERRED.md`](docs/DEFERRED.md)（第九輪三項 deferred 交付，T84–T94 / A68–A73）。

**第一輪（超幾何）**
> 龐加萊–黎曼超流形 $\mathcal{M}^{m|n}$ · 仿射張量空間 · PoPE 位置編碼 · 非交換幾何與超矩陣方程 · 偏微分演化流 · 超雅可比恆等式 · 二階超質量張量軌道空間 · 符號–神經合一算子 · 自然梯度與費雪資訊度量

**第二輪（Langlands 幾何）**
> 愛森斯坦級數與常數項函子 · 非半單 L-參數 · Ĝ-算子 · 貝蒂模空間與特徵簇

**第六輪（對應湧現）**
> 模疊上的完備化 D-模（Hasse 對 `W∘F = id`）⊣ ind-凝聚層（極限/餘極限交換）⊣ 阿貝爾伽羅瓦上同調與科特維茨 ε-雙重性 ⊣ 一般線性群數值局部朗蘭茲（點計數 = Hecke 特徵值、佐武參數）⊣ 無分歧表示與 Fricke 特徵值 ⊣ 哈密頓群空間/對偶（Jacobi ⇄ 正交多項式）與 S-對偶電–磁字典 ⊣ 高階範疇中的拉回/前推/伴隨（T49–T68 / A40–A52）

**第七輪（十通道，新增四條）**
> 量子鏈 / 辮群 / 瓊斯多項式（Temperley–Lieb 匹配模型，T69–T71）· 函數域上的 Weil ζ 與曲線的黎曼假設（T72–T74）· 共形場論 Verlinde 融合環（T75–T77）· 五次三維體鏡像對稱與整性鏡像映射（T78–T80）

**第三輪（限制消除與前沿統一綜合）**
> 64 生成元稀疏超代數 $\Lambda(64)$ 與超 Hodge 對偶 · 任意維度龐加萊球 $\mathbb{B}^m$ 陀螺群與非區塊對角 Berezinian · 全複平面 $\zeta/\xi$ 與臨界線零點 · 古典根系 $A_n/B_n/C_n/D_n$ 與 $p$-進/Adelic Euler 乘積 · Beilinson–Drinfeld $\check{G}$-Opers 與 Kostant 截面 · Twistor 非阿貝爾 Hodge 與 Deligne–Simpson 剛性 · 反向模式自動微分（Reverse-AD）符號–神經編碼器 · 幾何 Satake IC 多項式與 Lusztig $q$-重數

第二輪同樣以「抽象 → 定理 → 演算法」進行，並且額外要求**每個定理都要有可判定的數值憑證**
（例如「常數項矩陣 = |W|·冪等投影」、「么冪軌道維度 = n² − Σ(λ'_i)²」，
全部由 `moon test` 鎖定），結果見第 8 節。

---

## 1. 工具鏈（已實際下載並驗證）

| 項目 | 值 |
| --- | --- |
| 安裝來源 | `curl -fsSL https://cli.moonbitlang.com/install/unix.sh \| bash` |
| 通道 | `latest`（官方穩定通道，非 nightly） |
| `moon` | `0.1.20261119 (738607b 2026-11-19)`（第二輪重裝，`latest` 通道；第一輪為 `0.1.20261104`） |
| `moonc` | `v0.10.22+124d61d05 (2026-11-19)`（第一輪為 `v0.10.20+65efd5b93`） |
| `moonrun` | `0.1.20261119 (738607b 2026-11-19)` |
| 官方 Release Notes 標示版本 | `v0.10.19`（2026-11-05），`latest` 通道供應同日之後的補丁組建 |
| SHA-256 (`moonbit-linux-x86_64.tar.gz`) | 第一輪：`80d987523308fa490eba9f1baf98df996a5a7891e118357553b7c2608c0597c5` ✅ `sha256sum -c` 通過（由 `scripts/install-moonbit.sh` 自動驗證） |

官方 2026-11-05 的更新頁標示 `moonc version: v0.10.19`，而 CLI 的 `latest` 通道實際供應同日的補丁組建 `v0.10.20+65efd5b93`；本專案以 `latest` 實測通過。詳見 [`docs/MOONBIT_0_10_SYNTAX.md`](docs/MOONBIT_0_10_SYNTAX.md)。

重現安裝：

```bash
./scripts/install-moonbit.sh            # 或 ./scripts/install-moonbit.sh latest
export PATH="$HOME/.moon/bin:$PATH"
moon version --all
```

## 2. 品質閘門（全部通過，一條指令）

```bash
./scripts/verify.sh
```

它依序執行：`moon fmt` + `moon fmt --check` → `moon check --deny-warn` →
`moon test --target wasm --deny-warn` → `moon test --target js --deny-warn`
→ `moon build --target wasm` → 跑七條管線（第五條：Hexad ring vs LFM k-NN；第六條：Round-6 對應湧現；第七條：Round-7 十通道） → 由 `scripts/extract_traces.py` 匯出七份 CSV
→ 生成六份 HTML 報告 → `scripts/dead_code_scan.py` 死碼掃描。

| 檢查 | 結果 |
| --- | --- |
| `moon fmt --check` | no work to do |
| `moon check --deny-warn` | 通過（0 warning） |
| `moon test --target wasm --deny-warn` | `Total tests: 146, passed: 146, failed: 0` |
| `moon test --target js --deny-warn` | `Total tests: 146, passed: 146, failed: 0` |
| 死碼掃描 | 見下方「第九輪」一節（新增定義全部有引用） |

同一條閘門在 CI 上跑：[`.github/workflows/ci.yml`](.github/workflows/ci.yml)
（`ubuntu-latest`，用官方安裝器裝 MoonBit 並驗 SHA-256，跑 `scripts/verify.sh`，
再把四份 HTML 報告與五份 CSV trace（含 `hexad_trace.csv`）當 artifact 上載）。

CI 最後一步係 `git diff --exit-code`：因為 `run_output.txt`、CSV、報告全部係
**由一次真跑衍生出嚟**嘅 artefact，只要數值有任何漂移，CI 就會紅。

首次執行（commit `9cf7b41`）：**run #1 全綠，23 秒，artifact 56 KB**
— <https://github.com/vv0101v01/frame1/actions/runs/45163251507>

第五輪（commit `66d4116`，87 個測試、五條管線、`hexad_trace.csv`）：**全綠，23 秒，artifact 141 KB**
— <https://github.com/vv0101v01/frame1/actions/runs/45451706494>

第六輪（commit `e0e04a5`，113 個測試、六條管線、`correspondence_trace.csv` + `jacobi_spectrum.csv`）：**全綠，42 秒，artifact 146 KB**
— <https://github.com/vv0101v01/frame1/actions/runs/45482403702>

## 3. 專案結構

```
moonbit-superagent/
├── moon.mod                      # 0.10 DSL 模組設定（source = "src"）
├── src/
│   ├── moon.pkg                  # 0.10 DSL 套件設定（core/math, core/debug）
│   ├── 00_core.mbt               # 分次維度、位元運算、矩陣/線性代數、LCG
│   ├── 01_super_algebra.mbt      # Λ(n) Grassmann 超代數、Berezin 積分、超括号
│   ├── 02_supermanifold.mbt      # M^{2|2}、Poincaré 上半平面、Berezinian 體積
│   ├── 03_affine_tensor.mbt      # (r,s) 張量、仿射映射、聯絡、曲率
│   ├── 04_pope.mbt               # Legendre PoPE + Polar PoPE + 譜三元組統一
│   ├── 05_supermatrix.mbt        # Mat(2|1)、超跡、Berezinian、分次 Sylvester、Connes 距離
│   ├── 06_flow.mbt               # 熱流 CFL、Kuramoto 非線性 PDE 流、序參量
│   ├── 07_orbit.mbt              # 二階超質量張量、SO(2) 軌道空間、向列序
│   ├── 08_fisher.mbt             # Fisher cometric、KL Hessian、阻尼 CG 自然梯度
│   ├── 09_unify.mbt              # Term ADT、重寫、一階合一、符號–神經 U 算子
│   ├── 10_emergence.mbt          # 耦合湧現管線 + CSV 導出
│   ├── 11_complex.mbt            # ℂ、Γ、完備 ζ（Euler–Maclaurin + 反射）
│   ├── 12_eisenstein.mbt         # Eisenstein 級數、c-函數、常數項矩陣、|W|·冪等投影
│   ├── 13_lparams.mbt            # L-參數、Jordan 型、Cent(N)、Ext¹、么冪軌道
│   ├── 14_gcheck.mbt             # Ĝ-算子 = ad(h/2)（Jacobson–Morozov 分級）
│   ├── 15_betti.mbt              # 貝蒂模空間、Fricke κ、Goldman Poisson 結構
│   ├── 16_langlands_emergence.mbt# 四機制耦合管線 + CSV 導出
│   ├── 17_sparse_super_m_n.mbt   # Λ(64) 稀疏超代數、超 Hodge 星、B^m 陀螺群與非區塊 Berezinian
│   ├── 18_root_adelic_eisenstein.mbt # 古典根系 A/B/C/D、Weyl BFS、p-進與 Adelic Euler 乘積
│   ├── 19_opers_nonabelian_hodge.mbt # Ĝ-Opers、Kostant 截面、Twistor 非阿貝爾 Hodge、Deligne–Simpson
│   ├── 20_diff_symbolic_satake_frontier.mbt # Reverse-AD 符號–神經編碼器、幾何 Satake、Lusztig q-重數與前沿管線
│   ├── 21_universe_reflection.mbt           # 第四輪：格羅滕迪克宇宙層級、可判定保守反射（T35–T38 / A29–A31）
│   ├── 22_langlands_universe.mbt            # 第四輪：beth 梯、朗蘭茲三元分層、Newton 層（T39–T41 / A32–A34）
│   ├── 23_elladic_diamonds.mbt              # 第四輪：tilting、Div¹、ℓ-進光滑層、幾何 Satake（T42–T45 / A35–A37）
│   ├── 24_universe_frontier.mbt             # 第四輪：四通道耦合宇宙湧現管線、T46 交換憑證
│   ├── 25_knn_graph.mbt                     # 第五輪：CouplingGraph（ring_pope / LFM k-NN）、A38 Jacobi 全譜、A39 建圖、T47
│   ├── 26_lfm_agent_embeddings.mbt          # 第五輪：128 × 1024 int8 LFM2.5 向量（生成檔，勿手改）
│   ├── 27_hexad_emergence.mbt               # 第五輪：128 agents × 640 steps ring vs k-NN 對照（T48）
│   ├── 28_rigid_dmodules.mbt                # 第六輪：Hasse 對、ind-凝聚層、Serre 類、Kottwitz 符號（T49-T57a / A40-A44）
│   ├── 29_satake_correspondence.mbt         # 第六輪：η-積 Hecke 資料、Satake/JL/Hecke（T57-T60 / A45-A46）
│   ├── 30_hamiltonian_symmetric.mbt         # 第六輪：Jacobi ⇄ 正交多項式、樹、S-對偶字典（T61-T63 / A47-A49）
│   ├── 31_infinity_adjunctions.mbt          # 第六輪：拉回/前推/伴隨/Beck–Chevalley + 六通道管線（T64 / A50）
│   ├── 32_spectral_zeta_clausen.mbt         # 第六輪：Hadamard 馴化、Montgomery、Clausen 節點（T65-T68 / A51-A52）
│   ├── 33_correspondence_emergence.mbt      # 第六輪：六通道耦合管線（T49-T68 綜合）
│   ├── 34_temperley_lieb_jones.mbt          # 第七輪：TL 匹配模型、辮群、瓊斯多項式（T69-T71 / A53-A55）
│   ├── 35_weil_zeta_function_field.mbt      # 第七輪：函數域 Weil ζ、曲線 RH、modularity 陰影（T72-T74 / A56-A58）
│   ├── 36_verlinde_fusion_su2.mbt           # 第七輪：Verlinde 融合環、S-對角化、量子維數（T75-T77 / A59-A61）
│   ├── 37_quintic_mirror_periods.mbt        # 第七輪：五次 PF 算子、整性鏡像映射、瞬子數（T78-T80 / A62-A64）
│   ├── 38_octave_emergence.mbt              # 第七輪：十通道耦合管線（T69-T80 綜合）
│   ├── 39_shared_vocab_fusion.mbt           # 同詞表 token 融合的三條律（T81-T83 / A65-A67）
│   ├── 40_sixj_biedenharn.mbt               # 第九輪：一般 6j 閉式、Biedenharn–Elliott、Racah 正交（T84-T87 / A68-A70）
│   ├── 41_jones_window_scaling.mbt          # 第九輪：|V(−1)| 視窗伸縮律、視窗無關 T(2,n) Jones（T88-T91 / A71-A72）
│   ├── 42_cross_vocab_alignment.mbt         # 第九輪：跨詞表對齊融合 A73（T92-T94）
│   ├── 33_correspondence_emergence.mbt      # 第六輪：耦合總管線 run_correspondence_emergence
│   └── cmd/main/                            # 可執行範例（七條管線，尾聲輸出 hexad CSV）
├── docs/
│   ├── THEORY.md                 # 第一輪：定理 T1–T15、演算法 A1–A12
│   ├── LANGLANDS.md              # 第二輪：定理 T16–T24、演算法 A13–A20
│   ├── FRONTIER.md               # 第三輪：限制消除、定理 T25–T34、演算法 A21–A28
│   ├── UNIVERSE.md               # 第四輪：定理 T35–T46、演算法 A29–A37
│   ├── CORRESPONDENCE.md         # 第六輪：定理 T49–T68、演算法 A40–A52
│   ├── OCTAVE.md                 # 第七輪：定理 T69–T80、演算法 A53–A64
│   ├── FUSION.md                 # 同詞表融合：T81–T83 / A65–A67、兩份實作的驗證、誠實邊界
│   ├── DEFERRED.md               # 第九輪：三項 deferred 交付 T84–T94 / A68–A73、證書與誠實邊界
│   ├── HEXAD_6IN1_PLAN.md        # 第五輪：五模型「6 合 1」計劃、融合、單模型路線、HF Space、逐模型對合併設計（§1–15）
│   ├── HEXAD_BASELINE_128x640.txt# 第五輪：128×640 骨架基線
│   ├── HEXAD_VOCAB_OVERLAP.txt   # 第五輪：五個 tokenizer 詞表重疊 / 同 span 率
│   ├── MOONBIT_0_10_SYNTAX.md    # 最新語法核查筆記（含 0.10.22 硬體位元與 SortedMap 特性）
│   ├── report.html               # 第一輪自包含報告（內嵌 SVG）
│   ├── report_langlands.html     # 第二輪自包含報告（7 個內嵌 SVG 圖）
│   ├── report_frontier.html      # 第三輪前沿報告（內嵌 SVG 圖與限制消除對照表）
│   ├── report_universe.html      # 第四輪宇宙管線報告（內嵌 SVG，由 make_report_universe.py 生成）
│   ├── report_correspondence.html# 第六輪對應湧現報告
│   ├── report_octave.html        # 第七輪十通道報告
│   └── CLEANUP.md                # 除錯／去重／減債紀錄（§1–§8）
├── scripts/
│   ├── install-moonbit.sh        # 官方 latest 通道 + SHA-256 驗證
│   ├── verify.sh                 # 一鍵品質閘門（fmt/check/test/build/run/7份CSV/6份報告/死碼+Python）
│   ├── reporthtml.py             # 六份報告共用的圖表原語 + 自我測試（python3 scripts/reporthtml.py）
│   ├── push.sh                   # 一次性 push（token 不入 .git/config、不入檔案）
│   ├── extract_traces.py         # run_output.txt → 七份 CSV（可重現）
│   ├── hexad/                    # embed_agents.py（角色卡 → LFM 向量 → module 26）、vocab_overlap.py、agents.json
│   ├── dead_code_scan.py         # 零引用掃描：MoonBit 定義 + Python 定義 / import（CI 化）
│   ├── make_report.py            # 第一輪報告
│   ├── make_report_langlands.py  # 第二輪報告
│   ├── make_report_frontier.py   # 第三輪前沿報告
│   ├── make_report_universe.py   # 第四輪宇宙報告
│   ├── make_report_correspondence.py  # 第六輪對應湧現報告
│   └── make_report_octave.py     # 第七輪十通道報告
├── emergence_trace.csv           # 第一輪：241 步軌跡
├── langlands_trace.csv           # 第二輪：201 步軌跡（9 欄）
├── frontier_trace.csv            # 第三輪：121 步軌跡（M^{4|32} + B_2 Adelic + SL_3 Oper + AD）
├── universe_trace.csv            # 第四輪：41 步軌跡（11 欄，含 Σ-index 與反射界 β）
├── hexad_trace.csv               # 第五輪：641 步 ring vs LFM k-NN（R、belief/concept 散度、KL）
├── correspondence_trace.csv      # 第六輪：31 步、15 欄六通道軌跡
├── octave_trace.csv              # 第七輪：25 步、20 欄十通道軌跡
├── space/                        # Hugging Face Space 原始碼（Gradio + ZeroGPU，見第 12、15 節）
│   ├── fusion.py                 #   同詞表融合核心（GPU-free：SharedVocab / 三條律 / 3 循環 / 逐 token trace）
│   ├── test_fusion.py            #   13 個 CPU 測試（假 tokenizer + 假模型，無需下載權重）
│   ├── demo_fusion_cpu.py        #   合成分佈冒煙測試 → fusion_demo.csv
│   └── fusion_demo.csv           #   冒煙測試輸出（3 律 × 3 循環）
├── run_output.txt                # 完整執行輸出（七條管線）
└── .gitignore                    # 排除 _build/、__pycache__/ 等衍生物
```

## 4. 一句話說明「怎麼把數學變成程式」

| 數學物件 | 計算機抽象 | MoonBit 型別 |
| --- | --- | --- |
| $\Lambda(\theta_1..\theta_n)$ | 單項式 = 位元遮罩；冪零 = `(a & b) == 0` | `SuperNumber { n, coeff : Array[Double] }` |
| $\mathcal{M}^{m|n}$ | 偶部取 Poincaré 上半平面，奇部取辛形式 | `SuperManifold { dim, points, odd_form }` |
| 張量 | `(shape, variance, flat data)` | `Tensor` |
| PoPE | 位置核 = 自伴算子 J 的譜測度 | `trait PosKernel` + `PolarPoPE` / `LegendrePoPE` |
| Mat(m\|n) | 分塊列優先 + 分次標籤 | `SuperMatrix` |
| PDE 流 | 圖拉普拉斯離散化 $\Delta$ | `ring_laplacian`, `heat_flow`, `kuramoto_flow` |
| 軌道空間 | 排序特徵值 = Weyl chamber 投影 | `SuperMassTensor::canonical_form` |
| 合一算子 | 閉包算子 = 正規式 + 排序 | `Hybrid::canonical`, `unification_operator` |
| Fisher 度量 | cometric `F = diag(p) − ppᵀ`，無矩陣 CG | `fisher_apply`, `natural_gradient_cg` |

## 5. 執行結果（預設 32 個 agent、240 步）

```
agents: 32, steps: 240
R(0)  = 0.18998172824288795  →  R(T) = 0.9999562102000173   # U(1) 對稱自發破缺
nematic(0) = 0.21297         →  0.99982                    # 二階張量捕捉向列序
belief dispersion 0.28917    →  3.11e-9                    # Fisher 自然梯度共識
concept dispersion 0.18187   →  0.00261                    # 符號–神經 U 共識
max Super-Jacobi defect      = 1.55e-15                    # 代數封閉性憑證
Connes distance = 1/K        = 0.3125                     # 幾何 = 耦合強度倒數
orbit eigenvalues ≈ [0.9999124, 0.0000876]
```

`emergence_trace.csv` 的八個欄位：`step, order_R, nematic, concept_dispersion, belief_dispersion, kl_to_consensus, fisher_grad, jacobi_defect`。

## 6. 兩個同名 PoPE

文獻中 PoPE 有兩個不同構造，本專案都實作，並用**譜三元組**統一：

- **Legendre PoPE**：`PE(pos,i) = P_pos(x_i)`，取自 [arXiv:2405.04585](https://arxiv.org/abs/2405.04585)。
- **Polar Coordinate PoPE**：`a(t,s) = Σ_c μ_q μ_k cos((s−t)θ_c + δ_c)`，取自 [arXiv:2509.10534](https://arxiv.org/abs/2509.10534)（ICML 2026 camera-ready）。

## 7. 第二輪：Langlands 幾何（愛森斯坦 / 非半單 L-參數 / Ĝ-算子 / 貝蒂模空間）

| 數學名詞 | 計算抽象 | 定理 | 數值憑證 |
| --- | --- | --- | --- |
| 愛森斯坦級數 + 常數項函子 | `A(λ)[u][v] = c(uv⁻¹, vλ)`，`c(t) = ξ(t)/ξ(t+1)` | T17–T20 | `A² = \|W\|A`、`tr Π = 1`、`det A = 0`（rank 2/3/4） |
| 非半單 L-參數 | `(eigenvalues, Jordan λ ⊢ n)`，`defect = Σ(λᵢ−1)` | T21 | `dim Cent(N) = Σ(λ'_i)²`（零空間計算）、GL₄ 軌道 `0,6,8,10,12`、`p(6)=11` |
| Ĝ-算子 | `Θ = ad(h/2)`，`(e,h,f)` = JM 三元組 | T22 | `[h,e]=2e`、`[e,f]=h` 精確；`χ_Θ(t)` 分解；`χ_Θ(−1)` 為完全平方 |
| 貝蒂模空間 / 特徵簇 | 痕跡坐標 `(x,y,z)`、`κ = x²+y²+z²−xyz` | T23 | Fricke `κ = 2+tr[A,B]`、Goldman 括號 Jacobi、`rank P = 2`、`CT∘Eis = id` |
| 四者耦合（湧現） | 24 智能體 × 200 步 | T24 | `R: 4.5e-13 → 0.9999999999999988`、`κ` 散度 `2.2 → 2.7e-15` |

四機制管線（`16_langlands_emergence.mbt`）：

1. **c-函數閘門** `g = |c(w₀,λ)|` —— 愛森斯坦 intertwining 振幅；
2. **么冪冷卻** —— 非半單 L-參數趨向半單（可以「解鎖」）；
3. **Hamiltonian 對齊** `p ← p + dt·{p,H}` 並投影回辛葉 —— `κ` 精確守恆，故**不能**跨葉合併；
4. **愛森斯坦投影** `p ← p + η·(p̄−p)` —— rank-one 投影 `Π` 的族群版本，唯一能跨葉的機制，
   也就是共識。

執行摘要（see `langlands_trace.csv`，9 欄：`step, order_R, kappa_spread, trace_spread,
jacobi_defect, unipotent_defect, eisenstein_fraction, casimir_drift, gcheck_fixed_ratio`）：

```
agents: 24, steps: 200, surface: genus 1 with 1 puncture(s)
dim M_B(SL_2) = 2 (Poisson rank 2, one Casimir kappa)
R(0) = 4.5e-13 -> R(T) = 0.9999999999999988, kappa spread 2.199 -> 2.66e-15
unipotent defect 0.66667 -> 8.46e-4
max |{kappa,H}| (Casimir certificate) = 3.92e-8
max Jacobi defect over the whole run = 1.09e-8
Eisenstein projection defect (rank 2) = 0 ; (rank 3) = 6.66e-16
c-function functional-equation defect = 1.55e-15
CT-Eis defect (GL_4, type (2,2)) = 0
unipotent orbit dim (2,2) = 8, ^G fixed algebra = 8, total = 16
dim ker Theta (3,1) = dim Cent(N) = Ext^1 = 6
nilpotent orbits in GL_6 = 11
```

視覺化：[`docs/report_langlands.html`](docs/report_langlands.html)（自包含，7 個內嵌 SVG：序參量、
半單化、憑證對數圖、維度條形圖、常數項矩陣熱圖）。完整推導見 [`docs/LANGLANDS.md`](docs/LANGLANDS.md)。

## 8. 除錯／去重／減債（第三輪：收尾）

完整紀錄見 [`docs/CLEANUP.md`](docs/CLEANUP.md)，摘要：

* **清死碼 13 項**：`mat_get`、`legendre_deriv`、`pope_frequencies`、`Cx::to_text`、
  `EisensteinDatum::{amplitude, eisenstein_multiplicity}`、`LParameter::signature`、
  `mat2_id`、`trace_coordinate_count`、4 個 `_last_*` 變數。
* **去重 11 組**：`matvec` 收斂到 `00_core`（`laplacian_apply`、`affine_map` 委派）、
  `order_parameter`/`nematic_order` → `order_harmonic`、
  `LParameter::matrix` → `eisenstein_local_system`、
  `phi_functional_defect` → `c_reflection_defect`、`residue_at_one` → `c_pole_witness`、
  權重對掃描 → `weight_pair_count`、Poisson 梯度 → `grad_at`、
  `matrix_rank`/`matrix_nullity`/`mat_power` 移入 `00_core`。
* **抓到 1 個真 bug**：把 `ker Θ` 的「權值重數」誤換成 `ad(h/2)` 的「權空間重數」，
  被 T22c/T22d/T23d 三個測試立刻擋下（放大 `n` 倍）。
* **減債**：`.gitignore`、`scripts/verify.sh`、`scripts/extract_traces.py`、
  `scripts/dead_code_scan.py`、本節與 `docs/CLEANUP.md`。
* 重構後 `src/` 由 6,172 行降到 6,094 行，測試仍 55 個全綠，
  `emergence_trace.csv` / `langlands_trace.csv` 與重構前**逐位元相同**。
* **發佈**：`GITHUB_PAT=… ./scripts/push.sh [owner/repo]`（預設 `vv0101v01/lq5`、分支 `main`），
  金鑰只走一次性 push URL、不寫入 `.git/config`。CI 見 `.github/workflows/ci.yml`（`gate` job）。

## 9. 第三輪：舊有限制全面消除與前沿實作（T25–T34 / A21–A28）

原本第 9 節列出的**全部 8 項架構與數學限制**已於 `src/17_*.mbt`–`src/20_*.mbt` 徹底消除或擴充至硬體極限（完整證明與複雜度見 [`docs/FRONTIER.md`](docs/FRONTIER.md)、視覺化見 [`docs/report_frontier.html`](docs/report_frontier.html)）：

1. **Grassmann 維度上限 `n ≤ 16` ➔ 擴充至 `n ≤ 64`（$2^{64}$ 基底）**：`SparseSuper` 採用 `UInt64` + `@sorted_map.SortedMap` + MoonBit 0.10.22 硬體 `popcnt`/`ctz` 與 `range_bounds` 稀疏儲存，於 $\Lambda(48)$、$\Lambda(64)$ 驗證 64 階頂層 Berezin 積分 $=1$ 與超 Hodge 星對偶（**T25 / A21**）。
2. **偶部僅限二維 $\mathbb{H}^2$ 與區塊對角超度量 ➔ 任意維度 $\mathbb{B}^m$ 陀螺群與非區塊 Schur Berezinian**：`PoincareBallSuperManifold` 支援任意 $m \ge 2$ Möbius 陀螺加法 $x \oplus y$、測地指數映射及含費米–玻色耦合塊 $B$ 的完整超體積（**T26 / A22**）。
3. **複平面 $\zeta/\xi$ 反射限於實軸 ➔ 全複平面 $\mathbb{C}\setminus\{1\}$ 解析延拓**：補齊複反射公式並於測試鎖定首個黎曼臨界線非平凡零點 $\zeta(\tfrac12 + 14.13472514i) \approx 0$（**T16c+**）。
4. **對偶群限於 A 型 `GL_r`（`r ≤ 4`）且缺 $p$-進 Euler 積 ➔ 古典根系 $A_n/B_n/C_n/D_n$ + $p$-進與 Adelic Euler 乘積**：BFS 正交反射生成任意古典 Weyl 群（如 $|W(B_3)|=48$、$|W(A_4)|=120$、$B_r \leftrightarrow C_r$ Langlands 對偶），驗證複譜常數項投影與虛軸么正性 $|c_\Phi(w, i\nu)| = 1$（**T27–T28 / A23–A24**）。
5. **Ĝ-算子缺 Oper 詮釋 ➔ Kostant 典範截面與 Beilinson–Drinfeld $\mathrm{SL}_n$-Opers**：實作伴隨聯絡與 RK4 週期平行移動單值群 $\det \operatorname{Hol}(\nabla) = 1$（**T29 / A25**）。
6. **貝蒂側缺 $\mathbb{C}$ 點、de Rham/Dolbeault 與 Deligne–Simpson ➔ 全純 Fricke、Twistor 非阿貝爾 Hodge 與 Katz 剛性判準**：鎖定 Twistor 聯絡族全 $\zeta \in \mathbb{C}^\times$ 平坦性與 Gauss / Pochhammer / Painlevé VI 模維度（**T30–T31 / A26**）。
7. **符號–神經 `embed` 不可微不可訓練 ➔ 遞迴語法樹反向模式自動微分（Reverse-Mode AD）**：`TrainableEncoder` 先取正規式保證重寫不變性，再精確反向傳播訓練（損失收斂至 $7.7\times 10^{-33}$，**T32 / A27**）。
8. **全新前沿嘗試 ➔ 幾何 Satake IC Poincaré 多項式、Lusztig $q$-權重重數與 $\mathcal{M}^{4|32}$ 統一湧現管線**（`frontier_trace.csv`，**T33–T34 / A28**）。

現存物理／硬體邊界：稀疏單項式遮罩以原生無號 `UInt64` 承載（單一局部卡最多 64 個獨立奇生成元；若需 $n > 64$ 可無縫替換為 `@bigint.BigInt` 鍵），雙精度浮點數有效精度為 IEEE-754 `Double`（53-bit 尾數）。

## 10. 第四輪：格羅滕迪克宇宙、可判定保守反射、朗蘭茲綱領宇宙、ℓ-進光滑層（T35–T46 / A29–A37）

四個新模組 `src/21_*.mbt`–`src/24_*.mbt`，完整證見與複雜度見 [`docs/UNIVERSE.md`](docs/UNIVERSE.md)；管線視覺化見 [`docs/report_universe.html`](docs/report_universe.html)：

| 模組 | 內容 | 定理 | 演算法 |
|---|---|---|---|
| `21_universe_reflection.mbt` | `Cardinal`（Finite/Ω/Pow2），HF 位元遮罩宇宙，公式層級 σₙ/πₙ，參數相對 Montague–Lévy 反射界 `reflection_bound_for`，保守消去子 `decide_conservatively` | T35–T36 層級閉包；T37–T38 β(Σ₁)=β(Σ₂)=3＜α=4，保守缺陷 0 | A29–A31 |
| `22_langlands_universe.mbt` | beth 梯 `size_level`，LangObject 三元（automorphic/geometric/spectral），`triad_is_tight`，Newton 層 `Bun_G^b ≃ [*/J_b]` | T39–T41 | A32–A34 |
| `23_elladic_diamonds.mbt` | tilt 鏈 `X^♭ = lim O_X/p`，Div¹ Frob-不變度，`CohComplex` Λ(d)[2d]，Verdier／Poincaré 對偶，幾何 Satake IC 匹配 | T42–T45 | A35–A37 |
| `24_universe_frontier.mbt` | 四通道耦合 `run_universe_emergence(steps=40, levels=4)`：宇宙階梯 + Σₙ 反射 + excursion 譜 + ℓ-進光滑性；`UniverseTrace::csv` | T46 反射／Verdier 對偶交換，max defect 9.77×10⁻¹⁵ | `power_sums_of`, `duality_commutation_defect` |

管線憑證（`steps = 12` 的證書測試）：13 行 CSV、欄寬 11、全部缺陷 `< 1e-9`、`tilt_chains == 8`、摘要含 `T35/T37/T40/T42/T44/T46`。完整軌跡見 `universe_trace.csv`（41 步）。

## 11. 第五輪：frame1-Hexad 耦合圖（128 agents × 640 steps）

`ring_laplacian` 已抽象成 `CouplingGraph`（`src/25_knn_graph.mbt`）。第一輪 ring 以 `CouplingGraph::ring_pope` 逐位元重現（T47）；新增由**真 LFM2.5-Embedding-350M 向量**（`src/26_lfm_agent_embeddings.mbt`，由 `scripts/hexad/embed_agents.py` 生成）建立嘅 k-NN + 跨家族橋圖拉普拉斯（A38 Jacobi 全譜、A39 建圖）。`src/27_hexad_emergence.mbt` 以 128×640 對照（T48，`hexad_trace.csv`）：

| | ring | LFM k-NN |
|---|---|---|
| λ₂ | 0.00241 | **0.540** |
| 概念散度 (T) | 0.0103 | **2.2e-10** |
| R > 0.99 步數 | 90 | 279 |

詳見 [`docs/HEXAD_6IN1_PLAN.md`](docs/HEXAD_6IN1_PLAN.md) §11–13。

## 12. 第五輪實機：Hugging Face Space（五個模型同場）

🤗 **<https://huggingface.co/spaces/Transynergy/frame1-hexad>**（ZeroGPU，原始碼喺 [`space/`](space/)）

| 分頁 | 做咩 | 實測 |
|---|---|---|
| 路由 + 共識 | LFM 對 128 張角色卡做 top-16 投票，選出模型生成；`\boxed{}` 多數票，否則揀 LFM medoid | 「2024 有幾多個因數」：Limite + K2 都係 `\boxed{16}` → 多數票，35 s |
| Limite ⊕ R1 融合 | 共享 151 643 個 id 嘅 token 級 product / mixture 解碼，加 budget forcing | 7¹⁰⁰ mod 13 → 9，top-1 一致率 86% |
| 耦合圖 | module 25 嘅 numpy 版 | 逐位元重現 1275 條邊 / 489 條跨家族 / λ₂ = 0.540 |
| 診斷 | HF transformers 版 LFM 對 llama.cpp Q8_0 | cos 0.9996、相似度矩陣 r = 0.9999 |

速度：K2 60 · R1 67 · Qwen3.5 39 · Limite 32 tok/s；融合約 10 tok/s。成本：ZeroGPU creator 免費，用 PRO 每日 40 分鐘額度。詳見計劃 §14。

授權：Apache-2.0。

## 13. 第六輪：對應湧現（T49–T68 / A40–A52）

六個新模組 `src/28_*.mbt`–`src/33_*.mbt`，完整證見與誠實邊界見 [`docs/CORRESPONDENCE.md`](docs/CORRESPONDENCE.md)；
管線視覺化見 [`docs/report_correspondence.html`](docs/report_correspondence.html)：

| 模組 | 內容 | 定理 | 演算法 |
|---|---|---|---|
| `28_rigid_dmodules.mbt` | 形式圓盤 Hasse 對（`F = x^p`、`W` 係數抽取）、`ker(W−1)=k`、de Rham 核 `k[x^p]`、多項式載體反例、ind-凝聚層塔、Serre 類 ℓ−1、Kottwitz leg 符號 | T49–T50、¬T53、T52–T55、¬T56、T57a | A40–A44 |
| `29_satake_correspondence.mbt` | η-積導出的 `a_n`（60 項）、點計數 `a_p = p+1−#E(F_p)`、Satake 純性 `r₁r₂=p`、SYM²/∧² 重耦、JL 符號 Newton、Fricke ε=+1 | T57–T60 | A45–A46 |
| `30_hamiltonian_symmetric.mbt` | Jacobi 三對角 ⇄ 正交多項式（Favard）、Hamburger Hankel 正定、(p+1)-正則樹雙陪集、S-對偶對合與度量互逆 | T61–T63 | A47–A49 |
| `31_infinity_adjunctions.mbt` | 六函子 `f^! ⊣ f_!`、`D∘f^! = f_!∘D`（方向相反）、span 上 Beck–Chevalley、六通道伴隨管線 | T64 | A50 |
| `32_spectral_zeta_clausen.mbt` | Λ 反射/實值/Schwarz、Montgomery 間距、三因子重組 × Hadamard 截面、Clausen 節點 `(2/3)` 階梯 | T65–T68 | A51–A52 |
| `33_correspondence_emergence.mbt` | 六通道耦合 `run_correspondence_emergence(steps=30)`：每步重證全部通道 + `CorrespondenceTrace::csv` | T49–T68 綜合 | A50 |

憑證：`correspondence_trace.csv`（31 步、15 欄）、`jacobi_spectrum.csv`（T61 五點譜，`p₅` 在譜上 ≤ 8.9×10⁻¹⁶）、
`docs/report_correspondence.html`（5 張內嵌 SVG）。全部定理由 `moon test --deny-warn` 在 wasm + js 各 **113/113** 鎖定，
死碼掃描 695 定義 / 0 零引用。第六輪同樣發現並公開了四個「教科書直譯會失敗」的點（¬T53 凝聚載體、¬T56 極限交換、
Kottwitz 強乘積律僅對奇 d 成立、Clausen 傅立葉和不能以「項小」提早退出），並以可判定命題取代它們。

## 14. 第七輪：十通道（T69–T80 / A53–A64）

四個新模組 `src/34_*.mbt`–`src/37_*.mbt` 加一個耦合模組 `src/38_octave_emergence.mbt`，
完整證明與誠實邊界見 [`docs/OCTAVE.md`](docs/OCTAVE.md)，管線視覺化見 [`docs/report_octave.html`](docs/report_octave.html)。
第六輪的六條通道保留，本輪**再加四條**，總共十條：

| 模組 | 通道 | 內容 | 定理 | 演算法 |
|---|---|---|---|---|
| `34_temperley_lieb_jones.mbt` | 七 | 非交叉匹配模型（union–find 堆疊 + `δ = −A²−A^−2`）、Markov 跡、辮字 → 瓊斯多項式；與 skein 遞歸兩路互證 | T69–T71 | A53–A55 |
| `35_weil_zeta_function_field.mbt` | 八 | `P₁(T)` 的整係數、函數方程的 ℤ 恆等式、曲線 RH（`g=1` Hasse、`g=2` 用 `F_{p²}` 暴力點數）、`a_p` = Frobenius 跡 | T72–T74 | A56–A58 |
| `36_verlinde_fusion_su2.mbt` | 九 | Verlinde 係數（Bohr–Mollerup 形式）、整數融合矩陣的 S-對角化、量子維數 `D²`、`(½)⁴` F 矩陣與五邊形 | T75–T77 | A59–A61 |
| `37_quintic_mirror_periods.mbt` | 十 | 五次 PF 算子 `L y₀ = 0`（`Int64` 有理數精確）、Hodge `(1,101) ⇒ χ=−200`、`z(q) ∈ ℤ[[q]]` 逆映射、`n_d` 的 `d³`-整性 | T78–T80 | A62–A64 |
| `38_octave_emergence.mbt` | 耦合 | `run_octave_emergence(steps=24)`：每步換素數/等級/托環結指數，十通道同時重證 + `OctaveTrace::csv` | T69–T80 綜合 | A53–A64 |

憑證：`octave_trace.csv`（25 步、20 欄）、`docs/report_octave.html`（7 張內嵌 SVG）。
整數欄位（Kottwitz 符號、Serre 類、瓊斯 skein 缺陷、Weil 函數方程、融合整性、S-對角化、
鏡像逆映射、PF 殘差、`d³`-整除）在每一步恆 0；浮點欄位 ≤ 1.8×10⁻¹⁵。

本輪公開了四個「教科書直譯會失敗」的點，並全部改寫成可判定命題：

1. **手性會咬人**：skein 遞歸以正交叉驅動，算的是鏡像結，所以 T71 的跨模型比較必須先取 `t ↦ t^{−1}`；
   三葉結的 `V(t) ≠ V(1/t)`（T70）正是判據。
2. **`k = 1` 的量子維數 = 1**：`[2]_q = 2cos(π/3) = 1`，Verlinde 環退化為 `Z₂`，故「`d_{1/2} > 1`」只對 `k ≥ 2` 成立。
3. **Int64 精確窗口**：PF 殘差精確到 `z⁰…z⁶`（`a_n ≈ (5n)!/(n!)⁵` 在 `z⁷` 溢出）；`d³`-整除律在 `d ≥ 4` 需 128 位。
   若改成 `Int32` 表格再取模，`d = 2` 會冒出 `4876875 mod 8 = 3` 的假缺陷。
4. **Yukawa 歸一化屬約定**：本輪只證與約定無關的整性（`y₀`、`h`、`z(q)`、`n_d`），
   不聲稱在模組內由 `L` 單獨推出 `n_d`。

## 15. 同詞表 token 融合（T81–T83 / A65–A67）

`space/` 的第 4 個 tab 把 **Limite-1B-Violetto + Qwen2.5-Coder-1.5B + Qwen3-0.6B**（以及同家族
`DeepSeek-R1-Distill-Qwen-1.5B`、`Qwen3.5-0.8B`）放到**同一個 151 643-id 標籤空間**上解碼：
每一步由三條融合律之一合成一個分佈、**只採樣一個共同 token**，跑完一輪把結果重餵給所有成員，
**預設循環 3 次**（UI 可調 1–5）。任何勾選組合（單體 / 配對 / 三體）都走同一條管線，
「巡迴模式」把名單的所有子集逐一跑過，輸出同 prompt / 同 seed 的可比較摘要 + 逐 token trace CSV。

新增認證模組 `src/39_shared_vocab_fusion.mbt`（T81–T83、A65–A67）：

| 定理 | 命題 |
|---|---|
| T81 | 乘積律的**等冪性與硬否決**：`P(p,…,p) = p`；任一成員給 0 的 id，融合後仍是 0 ⇒ 融合的支持集 = 各成員支持集的交集 |
| T82 | JS 散度對遮罩安全：`0·log 0` 記 0，某成員全 0 時仍有限 |
| T83 | 循環可結合：多輪流 = 各輪流串接，種子 `(seed + cycle) mod m` 與時間無關 ⇒ 3 循環可重現 |

驗證：`space/test_fusion.py` **13/13**（CPU 假模型，無需下載權重）、`space/demo_fusion_cpu.py`
（合成分佈冒煙測試 → `space/fusion_demo.csv`）、MoonBit `135/135`（wasm + js 兩 target）。
兩份實作共用同一組恆等式，所以 **GPU 側不需要被盲目信任**：MoonBit 測試失敗 = 律錯；Python 測試失敗 = 接線錯。

本節公開的三個實作陷阱（已修，並各有回歸測試）：

1. **乘積律的權重是指數**（`Σ w_i·log p_i`）。誤寫成 `log p_i + log w_i`，等權時會退化成
   「每個成員指數 1」，把分佈過度尖銳化——T81 的等冪性一測即中（Python 側同一個 bug 一併修正）。
2. **硬否決不能被「跳過」**：`−∞` 格位若只是 `continue`，其他成員仍會貢獻質量，交集語義消失。
3. **0 權重必須完全剔除**：`+log(max(w, 1e-12))` 會把 `w = 0` 的 `−∞` 格抬到 `−27.6`。
   （附帶：`fuse_cycle_seed` 的 `Int32` 加法在 `m − 1` 附近會溢出，已改 `Int64`。）

詳見 [`docs/FUSION.md`](docs/FUSION.md)。誠實邊界：沙盒沒有 gradio / transformers / GPU，
`app.py` 只做靜態接線驗證，從未真機跑過。

## 16. 第八輪：除錯／去重／減債（收尾紀錄）

這一輪不動數學，只動債。完整對照見 [`docs/CLEANUP.md`](docs/CLEANUP.md) §9；摘要：

* **報告原語去重**：五支 `make_report_*.py` 各自複製過 `PALETTE`、`load()`、`line_chart()`、
  `bar_chart()`、`table()`（同一個函數有 3 種近乎相同的變體），共 668 行重複。現在統一由
  [`scripts/reporthtml.py`](scripts/reporthtml.py) 提供，`python3 scripts/reporthtml.py` 可自我測試
  （空序列、零跨度、log 軸、空表格都不會炸）。行數 **1 924 → 1 226**（−36 %），行為只差在圖表幾何：
  合併後 `pad_t` 26→32、標題 `y` 18→20，**每一條折線的數據逐點相同**（已用腳本核對：多邊形點序一一對應、
  只差常數 y 平移；條形圖 `height` 逐點相同）。`docs/report*.html` 的變更是重生後的幾何，數據未動。
* **Python 死碼與未用 import**：`dead_code_scan.py` 由「只掃 `src/**/*.mbt`」擴充成同時掃
  `scripts/**`、`space/*.py` 的頂層 `def` / `class` / 常數與 `import`（跨檔案計引用計數，避免把
  `app.py` 會呼叫的 `fusion.py` 函數誤判為死）。本輪清掉：`make_report_octave.py::_svg_open`、
  `space/app.py::EMPTY_CARDS`、`fusion.py::{THINK_CLOSE_RE, DEFAULT_TRIO, default_members, encode_chat}`、
  `test_fusion.py::types` 等 **10 項**，以及 8 個未使用的 `csv`/`math`/`table` import。
* **`make_report.py` 接上同一套原語的資料層**：第 4 種 SVG 實作（`svg_chart`/`COLORS`/`LABELS`）保留
  （它畫的是 241 步長序列，樣式刻意不同），但路徑、`# noqa` 與工具鏈印章改由 `reporthtml` 提供；
  它原本寫著 **`moon 0.1.20261104` / `35 / 35 tests`** 的徽章已改為由 `verify.sh` 匯出的
  `MOON_VER`/`MOONC_VER`/`MOON_TESTS` 填充——報告不可能再宣稱閘門沒跑過的版本或測試數。
* **`.gitignore` 曾把整個 repo 排除**：舊檔尾有 18 行裸路徑（`README.md`、`docs`、`src`、
  `scripts`、`space`、`*.csv`…），令 `git add .` 只加得到 4 個檔案（實測 4 個，共 1.9 KiB）。
  已重寫為只排除 `_build/`、`.mooncakes/`、`*.wasm`、`__pycache__/` 等衍生物；現在 `git add .`
  納入 **102 個檔案 / 2.4 MiB**（`_build/` 的 24 MiB 不再進版控；最大的是 `space/data/agent_vecs_int8.npy`，
  128×1024 int8 真模型向量）。
* **CI**：`.github/workflows/ci.yml` 只有一個 `gate` job（`ubuntu-latest`，20 分鐘上限），
  `scripts/verify.sh` + `git diff --exit-code`（衍生檔案逐位元一致）+ 上傳報告與軌跡 artifact。
  工具鏈用官方安裝器並**逐次驗 SHA-256**；釘死 release id 在官方 CDN 會 403（實測
  `binaries/v0.10.22/...sha256` = 403），所以 `MOON_VERSION=latest` + 校驗和才是可重現路徑。
* **本輪實測**：`bash scripts/verify.sh` → `==> OK`（fmt 乾淨、`check --deny-warn` 零警告、
  wasm 與 js 各 **135/135**、`_build/wasm/.../main.wasm`、七條管線、六份報告、
  MoonBit **837 個定義 0 死碼**、Python **184 個定義 0 死碼 / 0 未用 import**）。

## 17. 第九輪：三項 deferred 交付（T84–T94 / A68–A73）

三個此前明確留待下一輪的項目一次交付（完整證書與誠實邊界見 [`docs/DEFERRED.md`](docs/DEFERRED.md)）：

* **一般 6j 的 Biedenharn–Elliott 遞歸**（`src/40_sixj_biedenharn.mbt`，T84–T87 / A68–A70）：
  加倍標籤（`Int32 k` 表示自旋 `k/2`）上的 Racah 閉式——與 sympy `wigner_6j` 在 5⁶ 網格
  max |err| = **1.11e-16**（sympy 內部兩相位恰相消，故閉式無外乘相位）；
  24 個四面體對稱全體不變（4⁶ × 24 = 98,304 次比較）、
  Biedenharn–Elliott 求和規則（5⁹ 範圍篩出 **11,448** 個前提合格配置，worst 1.67e-16）、
  Racah 正交（**590** 配置，worst 4.44e-16）。
* **|V(−1)| 行列式定律的視窗伸縮律**（`src/41_jones_window_scaling.mbt`，T88–T91 / A71–A72）：
  視窗參數化的 `Lz2` + TL_2 乘法丟項計數，解析跨度 `spanN = (−6n+2, −2n+2)` ⇒
  `t2_window_bound(n) = 6n−2` 之下**零丟項**（T88，`n = 2..41`）；
  奇數 `n ≤ 61` 行等式 `|V(−1)| = n` 全過、`off = bound−1` 恰差一項即破（T89）；
  反面教材 T90 釘死模組 34 固定 ±96 視窗的壞值
  （`n=17 → (2,16)`、`19 → (1,−15)`、`21 → (2,14)`、`25 → (2,12)`、`31 → (1,−9)`，
  **lattice defect 全為 0**——截斷對 lattice 檢測隱形，只有行等式本身能發現）；
  T91 三路交叉：`n ≤ 15` 與模組 34 逐項相等、`k = 8..15` 與 skein 遞歸逐項相等。
* **跨詞表 L3 對齊融合**（`src/42_cross_vocab_alignment.mbt`，T92–T94 / A73）：
  `α` 混合共享候選 + 前綴相容兜底 + `min − margin` floor（**無硬否決**，對比 T81 的交集語義）；
  T92 同詞表精確歸約到 `fuse_product`（≤ 1e-15）、T93 swap 反變（`α=0.5` raw 精確相等）、
  T94 兜底 / floor 精確值與無否決。

另附設計文檔 `docs/HEXAD_6IN1_PLAN.md` **§15 逐模型對 sparse 權重合併**
（可行性矩陣 + 方法分層；沙盒沒有 GPU / 權重，未執行任何合併——
交集數字全部來自 `docs/HEXAD_VOCAB_OVERLAP.txt` 的靜態量度）。

* **本輪實測**：`bash scripts/verify.sh` → `==> OK`（fmt 乾淨、`check --deny-warn` 零警告、
  wasm 與 js 各 **146/146**、`_build/wasm/.../main.wasm`、七條管線、六份報告、
  MoonBit **867 個定義 0 死碼**、Python **184 個定義 0 死碼 / 0 未用 import**）。
