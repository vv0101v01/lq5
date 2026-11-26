# frame1-Hexad：六合一合成詳盡方案（128 agents × 640 steps）

> 基準記錄點：`1545a42`（round 4）。本方案**不沿用** `2581a4e` / `066ed0d` 嘅三模型骨架（qwen3-4b / phi-4-mini / gemma-4-e2b + mock 候選），改為由 `1545a42` 開新分支 `hexad-6in1` 重新做。
> 狀態：規劃文件，所有「實測」數字均有出處；「估算」會明確標示。

---

## 0. 一頁摘要

| 項目 | 決定 |
| --- | --- |
| 合成對象 | **frame1（骨骼）** + K2-Horizon-0.9B + Limite 1B-Violetto + LFM2.5-Embedding-350M + Qwen3.5-0.8B + DeepSeek-R1-Distill-Qwen-1.5B（8-bit） |
| 規模 | **128 agents、640 steps**（frame1 原預設 32 × 240） |
| 核心事實 | 五個神經模型**架構、維度、詞表全部唔同** → 權重直接合併（SLERP / TIES / DARE / Frankenmerge）**不可行** |
| 合成路線 | **L1 系統級合成**（一個入口、六個器官）→ **L2 輸出空間融合**（同詞表 logit 融合 + 跨詞表序列級共識）→ **L3（可選）蒸餾成單一權重** |
| frame1 角色 | 唔係「第六個模型」，而係**骨骼／外骨骼**：狀態空間、耦合動力學、共識、可判定憑證、驗證閘門 |
| 已實測 | frame1 四條管線喺 128 × 640 下**全部收斂**，wasm 單核共 13.4 秒（見 §3、`docs/HEXAD_BASELINE_128x640.txt`） |
| 已實作（v2） | `ring_laplacian` → **LFM k-NN 圖拉普拉斯**：λ₂ 0.0024 → **0.540**；概念散度 0.0103 → **2.2e-10**（§11.3） |

---

## 1. 六個組件盤點（已逐一核對官方 config / 模型卡）

| # | 組件 | 參數 | 架構要點 | 隱藏維 / 層 | 詞表 | 上下文 | 強項 | 已知弱點 |
|---|---|---|---|---|---|---|---|---|
| S | **frame1** | — | MoonBit 0.10.22，24 模組，10,871 行，82/82 測試 | — | — | — | 超幾何狀態空間、Kuramoto / Fisher / Eisenstein 共識、反射憑證 | 無 IO、無語言能力 |
| A | **K2-Horizon-0.9B**（IFM/MBZUAI） | 0.9B dense | `K2HorizonForCausalLM`，GQA 32/8，YaRN×16 | 1536 / 28 | **64,256** | 128K | 數學、程式、STEM、tool-use；多教師蒸餾；全開源 | 需新版 llama.cpp（b11454+）/ sglang |
| B | **Limite 1B-Violetto**（Paradigma） | 1B dense | `LimiteForCausalLM`：48 層、3:1 局部(1024窗)/全域、MUDD、value-embedding、sigmoid softcap | 1280 / 48 | **151,680**（Qwen2Tokenizer，實詞 151,667） | 131K | 競賽數學極強（AIME 2026 94.01%、BeyondAIME 74.25%） | 只適合單輪；**非數學問題會被當成數學題**；要 vLLM 0.26 + 自家 plugin、Python 3.12 |
| C | **LFM2.5-Embedding-350M**（Liquid AI） | ~354M | LFM2 雙向化 encoder，CLS pooling | 輸出 **1024-d** | LFM 詞表 | 512 tokens | 快（query p50 ≈ 1.5–7 ms）、跨語言檢索 | **官方 11 種語言唔包中文**；只出向量，唔生成 |
| D | **Qwen3.5-0.8B**（Alibaba） | 0.8B dense（+視覺） | 6×(3×Gated DeltaNet + 1×Gated Attention)，MTP | 1024 / 24 | **248,320** | 262K | 201 種語言（含中文/粵語）、原生多模態、tool calling | 細模型容易出現 thinking loop → 用 non-thinking |
| E | **DeepSeek-R1-Distill-Qwen-1.5B**（8-bit） | 1.5B | `Qwen2ForCausalLM`，GQA 12/2 | 1536 / 28 | **151,936**（Qwen2.5 BPE） | 128K | 長 CoT 推理、自我檢查 | 會語言混雜、冗長；8-bit 輕微損失 |

**關鍵觀察**

1. **B（Limite）同 E（R1-Distill）詞表已實測完全相同**：151,643 個基礎 token 嘅（字串, id）逐一相同（`docs/HEXAD_VOCAB_OVERLAP.txt`）→ logit 可以直接相加，零對齊成本。
2. **A（K2，64,000）同 E（R1）唔係同一個詞表**：共享 42,986 個 token 字串（K2 嘅 67.2%）；同一段文字切出嚟相同 token 嘅比例：英文 100%、程式 95.5%、數學 81.2%、粵語 77.5%、中文 64.7% → 可以做 **跨詞表 token 級融合**（§11.1），但要做字節對齊。
3. C 係 encoder，唔參與生成，負責**度量**：佢嘅 1024-d 餘弦相似度就係 frame1 耦合圖嘅邊權。
4. 隱藏維 1536 / 1280 / 1024 / 1536 互不相同、層數 28 / 48 / 24 / 28 → **否定權重合併**。

---

## 2. 合成總架構

```
                         ┌───────────────────────────── frame1-Hexad ─────────────────────────────┐
 使用者輸入 ──► [D] Qwen3.5-0.8B 介面器官：語言偵測、粵/中→英正規化、圖片解析、任務拆解
                         │
                         ▼
                [C] LFM2.5-Embedding 度量器官：任務向量 q ∈ R^1024、記憶檢索、路由分數
                         │
                         ▼
        ┌──────────── frame1 骨骼（MoonBit, 決定性） ────────────┐
        │ 128 個 agent 狀態 X_i = (θ_i, p_i, c_i, T_i, κ_i, ℓ_i)  │
        │  θ 相位（Kuramoto）   p 信念單形（Fisher 自然梯度）      │
        │  c 概念（符號–神經 U 算子正規式）  T 二階軌道張量         │
        │  κ Casimir 不變量（硬約束）  ℓ 宇宙層級（反射驗證等級）   │
        │ 耦合圖 W_ij = cos(e_i, e_j)（由 C 提供）→ 圖拉普拉斯 L_W   │
        └───────┬───────────────┬───────────────┬────────────────┘
                ▼               ▼               ▼
      [A] K2-Horizon     [B] Limite 1B     [E] R1-1.5B 8bit
      規劃 / 程式 / 工具   數學求解           長 CoT 批判 / 驗證
                │               │◄── 同詞表 logit 融合 ──►│
                └───────────────┴───────────────┬────────┘
                                                ▼
                    frame1 收斂：Eisenstein 投影共識 + 保守反射閘門 + moon test 憑證
                                                ▼
                    [D] Qwen3.5 以使用者語言（粵語/中文/英文）輸出最終答案 + 憑證摘要
                         └──────────────────────────────────────────────────────────────────┘
```

### 2.1 frame1 數學結構 → 合成語義對照

| frame1 模組 | 原本數學 | 在 Hexad 嘅意義 |
|---|---|---|
| `06_flow` Kuramoto 相位 θ | U(1) 對稱破缺，序參量 R | agent 嘅「立場」；R → 1 = 群體達成一致 |
| `08_fisher` 信念 p ∈ Δ^{K−1} | Fisher 度量自然梯度 | p_i = agent i 對 K 個候選答案簇嘅信心；以自然梯度向共識移動 |
| `09_unify` U 算子 | 重寫正規式 + 一階合一 | 將答案化成**符號正規式**（例如數學式化簡、程式 AST 正規化）再比較，避免「寫法唔同但係同一答案」 |
| `07_orbit` 二階超質量張量 | SO(2) 軌道、向列序 | 分歧嘅**結構**：兩派對立 vs 四散；向列序高但 R 低 = 兩極化，觸發 R1 仲裁 |
| `16_langlands` Casimir κ | Hamiltonian 流保 κ，不能跨葉 | 硬約束（題目條件、型別、單位）——局部對齊**不能**違反；只有 Eisenstein 投影（顯式共識）可以跨簇 |
| `16_langlands` 么冪冷卻 | 非半單 → 半單 | 「未完成推理缺口」數量隨步數衰減，作為收斂指標 |
| `20_frontier` Reverse-AD `TrainableEncoder` | 可微符號編碼器 | **線上學習路由權重**：邊個模型喺邊類題目可信 |
| `21_universe` 保守反射 `decide_conservatively` | Σₙ 反射界 β | **分級驗證**：ℓ=0 格式 → ℓ=1 數值檢查 → ℓ=2 執行/單元測試 → ℓ=3 跨模型共識；高層結論必須喺低層可判定先接受 |

---

## 3. 128 agents × 640 steps：已實測嘅骨架基線

以 `1545a42` 原始碼、只改入口參數，喺 sandbox 用 `moon run --target wasm` 實跑（完整輸出：`docs/HEXAD_BASELINE_128x640.txt`）：

| 管線 | 32×240（README 原值） | **128×640（實測）** |
|---|---|---|
| P1 emergence：R(0) → R(T) | 0.190 → 0.99996 | **0.0478 → 0.99996** |
| P1 belief dispersion | 0.289 → 3.11e-9 | **0.276 → 3.36e-9** |
| P1 concept dispersion | 0.182 → 0.00261 | **0.218 → 0.0103** ⚠ |
| P1 max Super-Jacobi defect | 1.55e-15 | **1.78e-15** |
| P1 λ₂(L)（環形圖） | — | **0.00241** ⚠ |
| P2 langlands：R(T)、κ spread | 24×200：0.9999999999999988、2.66e-15 | **0.9999999999999978、4.88e-15** |
| P2 unipotent defect | 0.667 → 8.46e-4 | **0.669 → 0** |
| P3 frontier：測地散度、AD loss | 16×120 | **0.486 → 0；0.281 → 7.70e-33** |
| P4 universe：T46 交換缺陷 | 40 steps：9.77e-15 | **640 steps：1.42e-14** |
| 總耗時（四條管線） | — | **13.4 s**（wasm，單核） |

**由基線得出嘅兩個設計結論**

1. **環形拉普拉斯唔夠用**：λ₂ ∝ 1/n²，128 個 agent 時只有 0.0024，概念散度只降到 0.0103（32 agent 時係 0.0026）。→ **必須將 `ring_laplacian` 換成由 LFM 向量建嘅 k-NN 圖拉普拉斯**（§5.2），目標 λ₂ ≥ 0.05。
2. 骨架動力學好平：640 步只要大約 13 秒，**成本幾乎全部喺 LLM 呼叫**。所以 640 步唔係每步都叫模型，而係每 32 步一個「神經輪」（§4）。

---

## 4. 128 個 agent 嘅分工同 640 步排程

### 4.1 Agent 分配（合共 128）

| 器官 | agent 數 | 角色 | 解碼設定（建議） |
|---|---|---|---|
| [A] K2-Horizon-0.9B | **40** | 規劃、程式、一般 STEM、工具呼叫 | T=0.6，top-p 0.95，reasoning parser 開 |
| [B] Limite 1B | **32** | 數學求解（只接數學子任務） | 單輪、無 system persona、T=0.6–0.8 多樣本 |
| [E] R1-Distill 1.5B 8-bit | **24** | 批判、逐步檢查、仲裁兩極化 | T=0.6、無 system prompt、限長 |
| [D] Qwen3.5-0.8B | **16** | 介面、翻譯、多模態解析、最終撰寫 | non-thinking、限長防 loop |
| [C] LFM2.5-Embedding | **16** | 記憶 / 檢索 / 路由 / 聚類（非生成） | normalize_embeddings=True，query/document 前綴 |

> 非數學題時，Limite 嘅 32 個 agent **自動降為 0 活躍**，名額由 K2 同 R1 補上（路由由 `TrainableEncoder` 學）。agent 總數維持 128，frame1 狀態維度唔變。

### 4.2 640 步四階段

| 階段 | 步數 | 神經輪（每 32 步一輪） | 內容 | frame1 主導機制 |
|---|---|---|---|---|
| Ⅰ 解析 | 0–63 | 2 輪 | D 正規化 + 拆題；C 檢索記憶；A 出計劃 | 宇宙層級初始化、κ（硬約束）寫入 |
| Ⅱ 發散 | 64–319 | 8 輪 | A/B 出候選；E 批判；C 聚類成 K 簇 | Kuramoto 弱耦合（保多樣性）、么冪冷卻 |
| Ⅲ 收斂 | 320–511 | 6 輪 | B+E 同詞表 logit 融合重評；兩極化時 E 仲裁 | Fisher 自然梯度、Eisenstein 投影、k-NN 拉普拉斯 |
| Ⅳ 固化 | 512–639 | 4 輪 | ℓ=0…3 分級驗證；D 撰寫輸出 | 保守反射閘門、T46 交換憑證 |

合共 **20 個神經輪**；輪與輪之間嘅 31 步係純 frame1 動力學（毫秒級）。

### 4.3 呼叫預算（估算）

| 模式 | 每題生成呼叫 | 每題輸出 token（平均 1.5K/呼叫） |
|---|---|---|
| 全員（每輪 112 個生成 agent） | 2,240 | ≈ 3.4M |
| **建議：每輪活躍上限 32**（按 Fisher 梯度範數 + 相位偏離揀） | **640** | **≈ 0.96M** |
| 經濟：每輪活躍 16 | 320 | ≈ 0.48M |

LFM 向量：128 × 20 = 2,560 次 / 題，按官方延遲屬秒級總量，可忽略。

---

## 5. 三層合成方法

### 5.1 L1：系統級合成（主線，必做）

* 單一入口 `hexad serve`（OpenAI 相容 API），外部睇落係**一個模型**。
* Python 編排器負責 IO 同推理後端；frame1 負責所有**狀態更新同決策**，保持決定性同可重現。
* 每輪 frame1 ↔ 編排器交換 JSON：
  * 入：`embeddings[128][1024]`、`answer_clusters`、`scores[128][K]`、`canonical_forms[128]`、`checks[128][ℓ]`
  * 出：`active_set`、`coupling W`、`beliefs p`、`R / nematic / dispersion / λ₂ / κ drift`、`accept|continue|escalate`

### 5.2 frame1 需要新增 / 修改嘅模組

| 檔案 | 內容 | 測試（新增） |
|---|---|---|
| `src/25_hexad_config.mbt` | `HexadConfig { agents:128, steps:640, round_every:32, roles, K_max:8 }`；角色表 | 預設值、角色總和 = 128 |
| `src/26_graph_laplacian.mbt` | 由外部相似度矩陣建 k-NN 對稱圖、`L_W = D − W`、λ₂ 冪迭代、CFL dt | 對稱、行和 0、λ₂ > 0（連通）、與 `ring_laplacian` 在環圖時一致 |
| `src/27_hexad_state.mbt` | 128 agent 狀態、外部分數注入 Fisher 更新、U 算子正規式簇合併 | 單形保持、KL 單調下降、同正規式必同簇 |
| `src/28_hexad_gate.mbt` | ℓ=0..3 分級驗證接 `decide_conservatively`；兩極化偵測（nematic 高、R 低） | 保守缺陷 = 0；人工兩極化樣本觸發仲裁 |
| `src/29_hexad_pipeline.mbt` | 640 步四階段排程、20 神經輪、CSV `hexad_trace.csv`（12 欄） | 128×640 憑證測試（見 §7） |
| `src/cmd/hexad/` | CLI：讀 JSON、走 32 步、寫 JSON | 黃金檔案 round-trip |

> IO 方式（stdin JSON / 檔案 / js target + node）要喺實作時按 MoonBit 0.10.22 實際支援核實，`docs/MOONBIT_0_10_SYNTAX.md` 同步更新。原有 82 個測試必須保持全綠、CSV 逐位元不變。

### 5.3 L2：輸出空間融合

**(a) 同詞表 token 級融合：Limite ⊕ R1（數學專線）**

* Phase 0 核實：載入兩個 tokenizer，逐一比對 id 0…151,642 嘅 token 字串；不一致嘅 id 列入遮罩。
* 解碼時：`log p(x) = α·log p_Limite(x) + (1−α)·log p_R1(x)`，只喺共同詞表上；α 由 frame1 信念 p 決定（Limite 簇信心高時 α↑）。
* 特殊 token（`<think>`、EOS 等）分別映射，唔融合。

**(b) 跨詞表序列級融合：K2 / Qwen3.5 / 其他**

* 每個候選答案 → U 算子正規式 + LFM 向量 → 歸入 K 個簇。
* 簇得分 = frame1 信念 p 嘅群體平均 × 驗證等級 ℓ 權重。
* **唔做跨詞表 token 對齊**（成本高、收益低），留待 L3 蒸餾處理。

**(c) 中文檢索缺口（LFM 唔支援中文）**

1. 預設：D（Qwen3.5）先將粵語/中文正規化成英文再交 C 嵌入；
2. 後備：以 sentence-transformers 用中英平行句對微調 LFM2.5-Embedding（官方有提供微調範例）；
3. 對照組：用 Qwen3.5 最後層 mean-pool 做中文向量，只作評估比較。

### 5.4 L3：蒸餾成單一權重（可選，最後做）

| 選項 | 學生模型 | 理由 | 代價 |
|---|---|---|---|
| L3-α | **Qwen3.5-0.8B** | 中文/粵語、多模態、262K；適合香港用戶 | 數學要由 Limite/R1 蒸入 |
| L3-β | **K2-Horizon-0.9B** | 訓練配方全開源、數理基礎較好 | 中文較弱 |

* 資料：L1 系統跑題庫，只保留**通過 ℓ≥2 驗證**嘅軌跡（frame1 憑證 = 資料過濾器）。
* 方法：序列級 KD（SFT）為主；同詞表教師（Limite、R1 → 若學生用 Qwen 系詞表）可加 token 級 KD；跨詞表用 universal logit distillation 作實驗項。
* 嵌入能力：學生加一個 1024-d 投影頭，以 LFM 向量做對比蒸餾。
* 產物：**一個約 0.8–0.9B 權重 + frame1 wasm 驗證器**；frame1 仍然保留作推理時閘門（佢本身只係幾百 KB wasm）。

---

## 6. 部署同資源

### 6.1 權重記憶體

| 組件 | 格式 | 大小 |
|---|---|---|
| K2-Horizon-0.9B | Q8_0 GGUF / BF16 | 1.15 GB / 2.16 GB（bartowski 檔案表） |
| Limite 1B | BF16 safetensors | ≈ 2.1 GB |
| LFM2.5-Embedding-350M | BF16 | ≈ 0.7 GB（估算） |
| Qwen3.5-0.8B | BF16 | ≈ 1.6 GB |
| R1-Distill-1.5B | 8-bit | ≈ 1.9 GB（估算） |
| **合計** | | **≈ 7.5–8.5 GB** |

### 6.2 KV cache（由 config 計，BF16）

| 組件 | 每 token KV | 備註 |
|---|---|---|
| K2-Horizon | **56 KB** | 28 層 × 8 KV heads × 64 |
| R1-Distill | **28 KB** | 28 層 × 2 × 128 |
| Limite | **12 KB**（全域層）+ 局部層固定 ≈ 36 MB/序列 | 36 個局部層窗口 1024 |
| Qwen3.5 | **12 KB** + DeltaNet 固定狀態 | 只有 6 層全注意力 |

例：每輪 32 個活躍 agent、每個 8K context，最壞情況（全部係 K2）≈ 14 GB；按 §4.1 比例混合約 6–8 GB。

### 6.3 建議硬件同後端

* **最低**：1 × 24 GB GPU（RTX 4090 / L4）＋ 32 GB RAM；每輪活躍上限 32、context 8K。
* **舒適**：2 × 24 GB；GPU0 = K2 + Qwen3.5 + LFM，GPU1 = Limite + R1（logit 融合要同機低延遲）。
* 後端：
  * Limite → vLLM 0.26 + Paradigma plugin（Python 3.12，獨立 venv）
  * Qwen3.5 → vLLM（≥ 0.17）
  * K2-Horizon → sglang 或 llama.cpp（≥ b11454）
  * R1 8-bit → llama.cpp Q8_0 或 vLLM W8A8
  * LFM → sentence-transformers 或 llama.cpp `--embeddings`
* 所有後端以 OpenAI 相容 HTTP 暴露畀編排器；版本、revision hash、SHA-256 全部寫入 `docs/HEXAD_RUNTIME.md`（延續 frame1「工具鏈已實際下載並驗證」傳統）。

---

## 7. 驗收憑證（全部要由 `moon test` / CI 鎖定）

| 類別 | 指標 | 門檻 |
|---|---|---|
| 骨架 | 原 82 測試 + CSV 逐位元不變 | 100% |
| 骨架 | 128×640：R(T) | > 0.999 |
| 骨架 | belief dispersion(T) | < 1e-6 |
| 骨架 | concept dispersion(T) / (0) | < 0.02（現環形圖係 0.047，要靠 k-NN 圖改善） |
| 骨架 | λ₂(L_W) | ≥ 0.05 |
| 骨架 | Super-Jacobi / Casimir / T46 缺陷 | < 1e-12 / < 1e-6 / < 1e-12 |
| 閘門 | 保守反射缺陷 | = 0 |
| 融合 | Limite⊕R1 共同詞表覆蓋 | 實測報告（預期 ≥ 151,643） |
| 端到端 | 每題 wall-clock（活躍 32） | 實測報告 |

**評測（真跑，唔准用 mock 結果當成果）**

| 能力 | 基準 | 比較組 |
|---|---|---|
| 數學 | AIME 2026、HMMT Feb 2026、MATH-500 | Limite 單獨；R1 單獨；128 票多數決（無 frame1）；Hexad |
| 推理 / 知識 | GPQA Diamond、MMLU-Pro | K2 單獨；Qwen3.5 單獨；Hexad |
| 中文 | C-Eval、CMMLU、粵語自建小集 | Qwen3.5 單獨；Hexad |
| 檢索 | 自建中英跨語言集 | LFM 原版 / 翻譯前置 / 微調版 |
| 消融 | 逐個拎走 A–E；環形 vs k-NN 圖；32/64/128 agents；240/640 steps | — |

---

## 8. 里程碑（每個都係一個 commit，CI 必須綠）

| M | 內容 | 產出 | 驗收 |
|---|---|---|---|
| M0 | 由 `1545a42` 開分支 `hexad-6in1`；入口參數化 128×640；寫入本文件同基線 | `docs/HEXAD_6IN1_PLAN.md`、`HEXAD_BASELINE_128x640.txt` | verify.sh 全綠 |
| M1 | 下載並核驗五個模型（revision + SHA-256）；各後端健康檢查；tokenizer 對比報告 | `scripts/hexad_fetch.sh`、`docs/HEXAD_RUNTIME.md` | 五個 `/health` OK；詞表對比表 |
| M2 | frame1 新模組 25–28 + 測試 | `src/25..28_*.mbt` | 新測試全綠、λ₂ ≥ 0.05 |
| M3 | 640 步四階段管線 + `hexad_trace.csv` + 報告 | `src/29_*.mbt`、`docs/report_hexad.html` | 128×640 憑證測試 |
| M4 | Python 編排器 + 真模型接入（L1） | `scripts/hexad/orchestrator.py`、`hexad serve` | 10 條真題端到端跑通，有日誌 |
| M5 | L2 融合（Limite⊕R1、序列級簇、中文檢索） | `scripts/hexad/fusion.py` | 數學子集比單模型好 |
| M6 | 全量評測 + 消融 | `docs/HEXAD_EVAL.md` + CSV | 所有數字可重現 |
| M7（可選） | L3 蒸餾成單一權重 | 學生 checkpoint + 卡片 | 對比 M6 |

---

## 9. 風險同對策

| 風險 | 影響 | 對策 |
|---|---|---|
| Limite 將非數學題當數學題 | 錯答、浪費算力 | 只經路由接數學子任務；`TrainableEncoder` 學拒收 |
| Qwen3.5-0.8B thinking loop | 卡死、超長 | 強制 non-thinking、max_tokens 上限、重複 n-gram 截斷 |
| LFM 無中文 | 中文聚類 / 路由失準 | §5.3(c) 三層方案 |
| K2 / Limite 要特定推理後端版本 | 環境衝突 | 每個後端獨立 venv / 容器；版本鎖定寫入文件 |
| 128 agent 環形圖混合慢 | 概念未收斂 | k-NN 圖拉普拉斯（M2 硬指標 λ₂ ≥ 0.05） |
| R1 8-bit 精度損失 | 批判質素下降 | 同 BF16 版做小集對照，損失 > 2 分就改用 BF16 |
| 128 票「假共識」（羊群效應） | 一致但錯 | 階段Ⅱ弱耦合保多樣性；ℓ≥2 驗證先可以接受；向列序偵測兩極化 |
| 上一手問題重演（mock 當成果） | 可信度 | 規定所有評測數字來自真跑日誌，CI 存 artifact |

---

## 10. 下一步需要確認嘅決定

1. 硬件：手上係邊張 GPU / 幾多 VRAM？（決定 §6.3 配置同每輪活躍上限）
2. 主要語言：粵語 / 中文為主？（決定 L3 學生揀 Qwen3.5 定 K2，以及中文檢索方案）
3. 「合成 1 個」要做到 L1（一個入口、系統級）定 L3（真正單一權重檔）？
4. 係咪批准喺 `1545a42` 開 `hexad-6in1` 分支並 push 去 GitHub？

---

## 11. 第二版修訂（2026-11-24）：K2⊕R1 token 級融合、4 GB VRAM、k-NN 圖拉普拉斯

### 11.1 K2 + R1 token 級 logit 融合（跨詞表，已實測對齊度）

K2（64,000 BPE）同 R1（151,643 BPE）**唔同詞表**，所以唔可以直接 logit 相加。實測（真 `tokenizer.json`，`scripts/hexad/vocab_overlap.py`）：

| 配對 | 共享 token 字串 | 備註 |
|---|---|---|
| K2 ∩ R1 | 42,986（K2 67.2% / R1 28.3%） | 要做字節對齊 |
| R1 ∩ Limite | **151,643（100%，id 都一樣）** | 直接相加 |
| R1 ∩ Qwen3.5 | 131,586（R1 86.8%） | 對齊容易 |
| K2 ∩ Qwen3.5 | 50,932（K2 79.6%） | |
| LFM | — | encoder，唔生成，唔參與 token 融合；佢負責建耦合圖 |

同一段文字，K2 同 R1 切出相同 token（同一字節區間）嘅比例：英文 100%、程式 95.5%、數學 81.2%、粵語 77.5%、中文 64.7%。

**演算法 A40：字節對齊 token 級融合（K2 ⊕ R1，可選 ⊕ Limite）**

1. 兩邊各自做一次 forward，攞 top-k（k=20）logprob（llama.cpp `n_probs`）。
2. 每個候選 token 轉成**原始字節串**（兩邊都係 ByteLevel BPE，可以無損還原）。
3. 候選集 = 兩邊候選嘅聯集。對每個候選字節串 s：
   * 兩邊詞表都有 s：`score = α·log p_K2(s) + (1−α)·log p_R1(s)`
   * 只有一邊有：另一邊攞「字節前綴相容」嘅候選（c 係 s 嘅前綴，或者 s 係 c 嘅前綴）入面最高嗰個 logprob；搵唔到就用 top-k 最低值再減 margin 做下限
4. 揀分數最高嘅 s，輸出佢嘅字節；兩個模型各自用**自己嘅分詞**將呢段字節加入 KV cache（多數情況 1 個 token，邊界唔同就 1–2 個）。
5. α 由 frame1 信念 p 決定（K2 家族 vs R1 家族嘅群體信心），數學題時再 ⊕ Limite（同 R1 共用 id，直接加權相加，冇額外對齊成本）。
6. 特殊 token（`<think>`、EOS、chat template 標記）唔融合，按角色映射。

代價：每輸出一個字要兩次 forward，兩個模型要同時放喺 VRAM。

### 11.2 「咁容量係咪就唔會疊加？」同 4 GB VRAM

**兩個意思要分開講：**

| | 會唔會疊加 | 原因 |
|---|---|---|
| **VRAM／記憶體** | **會疊加** | 輸出層融合唔共用權重，每個生成模型都要各自放入記憶體；frame1 自己幾乎唔使（wasm，128 個 1024 維向量只係 128 KB） |
| **能力** | **唔會疊加成大模型** | 唔係變成一個 4.5B 模型。融合得到嘅係互補同糾錯：錯誤唔相關嘅時候會好過最好嗰個單模型，但唔會多出新知識。要真正縮細記憶體，只有 L3 蒸餾，而容量上限就係學生模型（約 0.8–0.9B） |

**4 GB VRAM 預算（GGUF 檔案大小，來自各 HF 倉庫；R1 Q8_0 為估算）**

| 組件 | 量化 | 大小 |
|---|---|---|
| K2-Horizon-0.9B | Q4_K_M | 0.69 GB |
| R1-Distill-1.5B | **Q8_0（你指定嘅 8-bit）** | ≈ 1.9 GB |
| Limite 1B | Q4_K_M | 0.68 GiB ≈ 0.73 GB |
| Qwen3.5-0.8B | Q4_K_M | 0.53 GB |
| LFM2.5-Embedding | Q8_0 | 0.38 GB |
| **全部放入 GPU** | | **≈ 4.2 GB → 4 GB 放唔落** |

**建議 4 GB 配置（方案 C，穩陣）**

| 位置 | 內容 | 佔用 |
|---|---|---|
| GPU（常駐） | K2 Q4_K_M + R1 Q8_0（token 融合嗰對） | ≈ 2.6 GB |
| GPU（KV cache） | q8_0 KV：K2 + R1 每 token ≈ 45 KB（f16 時係 84 KB） | ≈ 0.9 GB → 約 2 萬 token，即大概 4 條融合串流 × 5K context |
| GPU（按階段換入） | Limite Q4_K_M：只喺數學輪載入（mmap 由 SSD 載入，秒級） | 0.73 GB（同其他輪錯開） |
| CPU / RAM | Qwen3.5 Q4_K_M（介面、翻譯；呼叫少）＋ LFM Q8_0（向量；sandbox 2 vCPU 實測 128 條約 21 秒） | ≈ 0.9 GB RAM |
| CUDA context | | ≈ 0.3 GB |

影響：每輪 32 個活躍 agent 要分 8 批，每批約 4 條串流，所以 4 GB 可以跑，但比 24 GB 慢好多。128 agent × 640 步本身唔受影響，因為 frame1 動力學喺 CPU 跑，每 128×640 一次只要數秒。

**「全部都做 token 級融合」嘅限制**：四個生成模型（K2、R1、Limite、Qwen3.5）權重加埋 ≈ 3.85 GB，再加 KV 就超出 4 GB。所以 4 GB 上面：K2⊕R1 做 token 級融合（數學輪 ⊕ Limite），Qwen3.5 做序列級；LFM 本身就唔可以做 token 融合（佢唔生成）。

### 11.3 已實作：`ring_laplacian` → LFM k-NN 圖拉普拉斯

**做法**

1. `scripts/hexad/embed_agents.py`：128 張 agent 角色卡（40 K2、32 Limite、24 R1、16 Qwen3.5、16 LFM）用**真嘅 LFM2.5-Embedding-350M**（官方 Q8_0 GGUF，sha256 `6ec5f8e8…2599`，llama.cpp b12676）計 1024 維向量，再以 int8 存入 `src/26_lfm_agent_embeddings.mbt`（餘弦誤差最大 1.6e-3，8-NN 一致率 99.7%）。
2. `src/25_knn_graph.mbt`：
   * `CouplingGraph` 抽象。所有共識通道（Kuramoto、Fisher 信念、符號–神經合一）統一經佢；
   * `CouplingGraph::ring_pope` 重現第一輪 ring，**逐位元相同**（T47 測試鎖定）；
   * `CouplingGraph::knn`：k=8、熱核 `exp((cos−1)/τ)`、τ=0.25，再加**跨家族橋**（每個 agent 連去其他每個模型家族最相似嗰個），最後用 Kruskal 補連通；
   * `sym_eigenvalues`：Jacobi 全譜（A38），用 ring 解析譜驗證到 1e-10。
3. `src/10_emergence.mbt`：`run_emergence(cfg)` = `run_emergence_on(cfg, ring_pope)`。新增 `run_emergence_on(cfg, graph)`。
4. `src/27_hexad_emergence.mbt`：128×640 ring vs k-NN 對照、`hexad_trace.csv`、T48 憑證測試。

**點解要跨家族橋（實測）**：LFM 向量入面，同一家族內平均餘弦 0.86–0.89，跨家族得 0.25–0.39。

| 建圖方案 | 連通分量 | λ₂ |
|---|---|---|
| ring | 1 | 0.0024 |
| 純 k-NN（k=8） | **5（斷開）** | 0 |
| k-NN + 最小生成樹補橋 | 1 | 0.0017（仲差過 ring） |
| **k-NN + 跨家族橋（採用）** | 1 | **0.540** |

**128 agents × 640 steps 實測（`hexad_trace.csv`）**

| 指標 | ring | **LFM k-NN** |
|---|---|---|
| λ₂(L) | 0.00241 | **0.540（×224）** |
| 邊數 / 跨家族邊 | 128 / 5 | 1,275 / 489 |
| 概念散度 (T) | 0.0103 | **2.23e-10** |
| 概念散度 < 0.02 嘅步數 | 227 | **17** |
| 概念散度到達 0.0103 嘅步數 | 640 都未到 | **33** |
| 信念散度 (T) | 3.36e-9 | 2.58e-9 |
| 信念散度 < 1e-6 嘅步數 | 57 | 52 |
| R(T) | 0.99996 | 0.99960 |
| R > 0.99 嘅步數 | **90** | 279 ⚠ |
| 熱流 CFL dt | 0.714 | 0.098 ⚠ |
| Super-Jacobi 最大缺陷 | 1.78e-15 | 1.33e-15 |

**代價（如實記錄）**

* Kuramoto 相位同步變慢（90 → 279 步）：原本 ring 版本嘅 Kuramoto 其實係 PoPE **全連接**，k-NN 版本係稀疏鄰居平均，而且跨家族邊權得大約 0.06。之後可以調 τ，或者 Kuramoto 用「k-NN + PoPE 全連接」混合。
* λ_max 由 4 升到 29.2，顯式熱流步長要細 7.3 倍。湧現管線本身冇用熱流，只係會報告呢個數。

**品質閘門**：`./scripts/verify.sh` 全綠。wasm 87/87、js 87/87（原 82 + 新 5），0 warning，584 個定義 0 死碼。**原本四份 CSV 同四份 HTML 報告逐位元不變**，`run_output.txt` 只係喺尾加咗第五段。

---

## 12. 第三版：收斂成「1 個模型」嘅路線（融合 → 去重疊 → 內聯自學習）

### 12.1 總流程

```
[教師層]  Limite ⊕ R1  ── token 級 logit 融合（同詞表 151,643，id 完全一樣）
          K2、Qwen3.5  ── 序列級融合（答案層，經 U 算子正規式 + LFM 聚類）
          LFM          ── 唔生成，只提供 k-NN 耦合圖（後期由學生自己嘅向量頭取代）
                │
                ▼  frame1：128 agents × 640 steps 共識 + ℓ≥2 驗證  →  只保留「已驗證」樣本
[去重疊]  ① 集成權重去相關   ② 數據去重   ③ 專家 LoRA 以 TIES/DARE 合併（drop 50%）
                │
                ▼
[學生]    單一權重（建議 R1-Distill-1.5B 架構，同 Limite 共用詞表）
                │
                ▼  內聯自學習：學生 × 128 個角色提示 → frame1 共識/驗證 → LoRA 更新 → 重複
[產物]    1 個權重檔 + frame1 wasm 骨骼（推理時 128 agents = 同一個模型嘅 128 個角色樣本）
```

### 12.2 「減去 5 成疊項」三種解讀，邊種得、邊種唔得

| 層面 | 做法 | 得唔得 | 說明 |
|---|---|---|---|
| ① 集成權重 | 喺驗證集量度各教師嘅**錯誤相關矩陣 Σ**，權重 `w = Σ⁻¹1 / (1ᵀΣ⁻¹1)` | ✅ 得 | 錯誤高度重疊嘅模型會被自動減權（減去「重複投票」）。減幾多要由數據決定，**唔應該硬性設 50%** |
| ② 數據 | 所有教師一致 **而且** 學生已經答啱嘅樣本，大部分剔走；保留 10–20% 做 replay | ✅ 得 | 訓練集中喺分歧同學生弱項；可能剔走約一半，但比例要實測 |
| ③ 權重（跨教師） | 直接將 Limite / R1 / K2 / Qwen3.5 權重減重疊再合併 | ❌ 唔得 | 架構、維度、層數全部唔同（1280×48 / 1536×28 / 1536×28 / 1024×24），冇對應關係 |
| ③′ 權重（同一學生） | 喺**同一個學生基座**上分別訓練 數學 LoRA（Limite⊕R1）、中文 LoRA（Qwen3.5）、程式/工具 LoRA（K2），再用 **TIES / DARE（drop p=0.5，重縮放 1/(1−p)）** 合併成一個 | ✅ 得 | 呢個先係「減去 5 成疊項」喺權重層面最合理嘅落點：去冗餘、解符號衝突，最後合成一個權重 |

### 12.3 學生基座選擇

| 候選 | 詞表 | 優點 | 缺點 |
|---|---|---|---|
| **R1-Distill-Qwen-1.5B 架構（建議）** | 151,643（= Limite） | 可以對 Limite⊕R1 融合分佈做 **token 級 KD**；Qwen2 架構，訓練工具最齊；8-bit ≈ 1.9 GB | 中文、多模態要靠蒸餾補 |
| Limite 1B | 151,643 | 數學最強嘅起點 | 自訂架構（MUDD / XSA / VE），訓練工具支援少；輕度指令調校，加通用能力容易蝕數學 |
| Qwen3.5-0.8B | 248,044 | 中文、多模態、262K | 同 Limite/R1 唔同詞表，只可以做序列級 KD |
| K2-Horizon-0.9B | 64,000 | 配方全開源 | 同 Limite/R1 唔同詞表 |

### 12.4 內聯自學習（ReST / STaR 類，frame1 做裁判）

1. 每題：學生用 128 個角色提示各出一個樣本（LFM 圖由學生自己嘅 1024 維向量頭建）。
2. frame1：k-NN 圖共識 + 保守反射分級驗證（ℓ=1 數值、ℓ=2 執行／單元測試）。
3. 只用**通過 ℓ≥2** 嘅答案更新 LoRA；共識高但未驗證嘅只作低權重。
4. 每輪喺固定 held-out 集評測；冇進步就停。每輪混入教師數據，防止模型崩壞（model collapse）。
5. 只喺**可驗證領域**（數學、程式、有標準答案嘅題）做自學習。開放式文字唔做，因為冇裁判就會自我強化錯誤。

### 12.5 4 GB VRAM 可行性

| 階段 | 4 GB 得唔得 | 做法 / 限制 |
|---|---|---|
| 教師數據生成 | ✅（慢） | §11.2 方案 C |
| token 級 KD 數據 | ✅ | 只存每個 token 嘅 top-32 logprob（全詞表 151K 存唔落） |
| 學生訓練 | ⚠ 勉強 | 1.5B 4-bit QLoRA、序列 1–2K、gradient checkpointing；**151K 詞表嘅 logits 好食記憶體**（2048 token × 151,936 × 4 bytes ≈ 1.2 GB），要用分塊 / fused cross-entropy。建議訓練期租 24 GB GPU，部署返 4 GB |
| 內聯自學習 | ⚠ 慢 | 生成同 LoRA 更新交替 |
| 部署 | ✅ 寬裕 | 單一模型 8-bit ≈ 1.9 GB，KV 有大約 1.5 GB → 同時跑嘅 agent 比多模型方案多好多。**記憶體唔再疊加** |

### 12.6 誠實評估

* **做唔做到「1 個模型」？做到。** 蒸餾 + LoRA 合併之後就係一個權重檔；frame1 仍然係骨骼（推理流程同驗證器），唔計入權重。
* **「創新」喺邊？** 個別技術（蒸餾、ReST、TIES/DARE、self-consistency）都有前例。新嘅係組合：用 frame1 圖拉普拉斯共識同分級反射驗證，做**蒸餾數據過濾器**同**自學習裁判**；去相關集成權重；異詞表教師分別蒸成專家 LoRA，再用 DARE 合併。佢係唔係真創新，要靠 §7 嘅基準同消融實驗證明。
* **能力上限：** 1.5B 學生**唔會**等於五個教師加埋。合理預期：數學接近 Limite⊕R1 融合（唔保證超越），中文同通用能力好過 Limite 同 R1 單獨用，但未必及得上 K2 / Qwen3.5 喺佢哋各自強項嘅表現。全部要實測。
* **主要風險：** 學數學嘅時候中文同通用能力被蝕（用 replay + 分開 LoRA 處理）；自學習崩壞（用 ℓ≥2 驗證同 held-out 早停處理）；4 GB 訓練太慢（訓練期租 GPU 處理）。

---

## 13. 第四版：LFM 向量去向、「相搭度」、雲端成本

### 13.1 LFM 轉做圖之後，向量保唔保留？可唔可以傳畀新模型？

* **而家已經保留**：128 × 1024 向量存喺 `src/26_lfm_agent_embeddings.mbt`（int8，8-NN 一致率 99.7%），frame1 張圖係決定性、可重現。
* **但推理時仲要 LFM**：新題目、新答案要即時計向量先可以建圖。
* **可以傳畀新模型（向量頭蒸餾）**：
  1. 學生加一個 1536→1024 投影頭（約 157 萬參數），用 last-token pooling。
  2. 損失分三部分：(a) 直接對齊 LFM 向量（cosine），保持同一個 1024 維空間，舊記憶、舊圖照用；(b) 關係蒸餾，令 batch 內相似度矩陣 ≈ LFM，**保留 k-NN 圖結構**，呢樣先係 frame1 真正需要嘅；(c) InfoNCE 保留檢索能力。
  3. 中文：LFM 唔識中文，用中英平行句，以英文句嘅 LFM 向量做中文句嘅目標（或者改用 Qwen3-Embedding 做中文老師，見 13.2）。
  4. **骨幹凍結**，只訓練投影頭（加少量 LoRA），生成能力完全唔受影響。
  5. 驗收：學生向量建出嚟嘅圖同 LFM 圖嘅 8-NN 一致率 ≥ 0.9、λ₂ 差距 ≤ 10%（沿用 `embed_agents.py` 嘅一致率指標）。
* 完成之後就可以移除 LFM：**1 個權重檔同時做生成同埋提供 frame1 建圖嘅向量**。

### 13.2 係咪越「相搭」，留到畀新模型嘅功效就越多？

**大致係，但要分三層，而且太相似反而冇用。**

| 相搭層次 | 可以傳幾多 | 例子（已實測） |
|---|---|---|
| 詞表 id 一樣 | **完整機率分佈**（token 級 KD，連「暗知識」都傳到） | Limite、Qwen2.5-Coder-1.5B、Qwen3-Embedding-0.6B、Qwen3-0.6B 同 R1：151,643 個 id **完全一樣** |
| 詞表唔同 | 只傳到最後輸出嘅文字（序列級 KD） | K2（64K）、Qwen3.5（248K） |
| 形狀一樣 | 再加權重層操作（初始化、合併、移植層） | Qwen2.5-Coder-1.5B 同 R1：1536 × 28、8960、12/2 heads 完全一樣，**但** RoPE base（1e6 vs 1e4）同 embedding 綁定唔同 → 權重合併可以試，但風險高 |
| 同一祖先 + 細改動 | TIES / DARE 合併效果最好 | 同一學生上面訓練嘅專家 LoRA（§12.2 ③′） |

**反效果**：如果能力同錯誤都太相似（例如全部都係 Qwen 預訓練數據），「疊項」就會多，冇新嘢可以傳。最理想係**介面相搭（詞表、形狀）+ 能力互補（錯誤唔相關）**。仲有一個天花板：學生只有大約 1.5B 容量。

**可選嘅「Qwen 詞表家族」配置（只係建議，未改你嘅選擇）：**

| 角色 | 而家 | 可換成 | 好處 |
|---|---|---|---|
| 學生 | R1-1.5B 架構 | — | — |
| 數學 | Limite | 保留 | token KD ✓ |
| 程式 / 工具 | K2（64K 詞表） | Qwen2.5-Coder-1.5B-Instruct | 由序列級變 token 級 KD；形狀同學生一樣 |
| 向量 / 建圖 | LFM（無中文，512 ctx） | Qwen3-Embedding-0.6B | 同詞表、同樣 1024 維、100+ 語言（含中文）、32K ctx |
| 中文 / 多模態 | Qwen3.5-0.8B | 保留（序列級） | 多模態傳唔到純文字學生，除非學生加視覺塔 |

### 13.3 Hugging Face 成本：CPU $0.03 vs T4 $0.40

官方價（Spaces）：CPU Upgrade 8 vCPU / 32 GB RAM **US$0.03/h**；Nvidia T4 small 4 vCPU / 15 GB RAM / **16 GB VRAM** **US$0.40/h**；T4 medium 8 vCPU / 30 GB RAM **US$0.60/h**。GPU Space 唔會自動停，要手動暫停（T4 small 開足 30 日 = US$288）。

| 工作 | 8 核 CPU（$0.03） | T4 16 GB（$0.40） |
|---|---|---|
| frame1 128×640 | ✅ 數秒 | 唔需要 |
| LFM / Qwen3-Embedding 向量 | ✅（sandbox 2 vCPU 實測 128 條約 21 秒） | 唔需要 |
| 編排器、小規模測試、評測 | ✅ | — |
| 教師大量生成數據（百萬 token 級） | ❌ 太慢（估計單串流每秒十幾 token） | ✅ |
| Limite⊕R1 token 融合 | ⚠ 只可以做功能測試 | ✅ |
| 學生 QLoRA 訓練 / 自學習 | ❌ 實際上唔可行 | ✅（T4 冇 bf16，要用 fp16；冇 FlashAttention-2） |
| 部署最終 1.5B 模型（低並發 demo） | ✅（GGUF Q8，慢） | ✅ |

**建議混合用法**：CPU Space 長開（每日約 US$0.72），負責 frame1、編排同向量；T4 只喺生成數據同訓練時先開，用完即刻暫停。第一步先開 T4 跑 1 小時試點（約 US$0.4），實測 tok/s 之後再估總成本。

## 14. 第五版：Hugging Face Space 實機示範（2026-11-24）

**網址**：<https://huggingface.co/spaces/Transynergy/frame1-hexad>。原始碼放喺 [`space/`](../space/)，同 Space repo 同步。

### 14.1 架構

| 項目 | 選擇 | 原因 |
|---|---|---|
| 硬件 | ZeroGPU（RTX PRO 6000 Blackwell，半卡 48 GB） | creator 免費；PRO 帳戶最多可以開 10 個 |
| 精度 | 五個模型都用 bf16（約 9.5 GB） | VRAM 充裕，唔使量化。R1 8bit 只係 4 GB 方案先需要 |
| 依賴 | `transformers==5.15.0`、Python 3.12、torch 由 runtime 提供（2.11） | K2 驗證過 5.15；Limite remote code 嘅 API 靜態核對通過 |
| 版本 | 五個 repo 全部釘死 commit SHA | 有兩個 repo 帶自訂模型碼（K2、Limite） |

### 14.2 四個分頁

1. **路由 + 共識**
   - LFM 將 `query:` 向量同 frame1 128 張角色卡比較，取 top-16，用 softmax（T = 0.1）投票，得出家族權重。Auto 模式揀頭兩名。
   - 共識分兩步：先睇 `\boxed{}` 有冇多數票（可驗證領域，對應 §12）。冇多數就喺「完整答案」入面揀 LFM cosine 嘅 medoid。俾 token 上限截斷咗嘅答案唔可以做 winner。
2. **Limite ⊕ R1 token 級融合**
   - 兩個模型喺共享嘅 151 643 個常規 id 上，用 product（`α·log p_L + (1−α)·log p_R`）或者 mixture 混合。
   - 特殊 token 收窄成兩個 pseudo-token：END 同 `</think>`。`</think>` 喺 R1 係 151649，喺 Limite 係普通文字序列。
   - 有 budget forcing：用到 75% 預算就強制輸出 `</think>`。
3. **耦合圖**：module 25（A39）用 numpy 重寫，喺儲存咗嘅 int8 向量上逐位元重現 1275 條邊、489 條跨家族、λ₂ = 0.53975、CFL 0.0978。
4. **診斷**
   - LFM 官方 remote patch 係為 transformers 4.56 寫嘅。5.15 嘅 `Lfm2ShortConv` 多咗 `seq_idx` 參數，又多咗 `create_recurrent_attention_mask`，所以喺 Space 入面重寫咗雙向 attention 同非因果 short-conv。
   - 自檢結果：同 llama.cpp Q8_0 向量比，cos 平均 **0.9996**（最低 0.9995），相似度矩陣 Pearson r = **0.9999**。

### 14.3 實測（gradio_client，2026-11-24）

| 測試 | 結果 |
|---|---|
| 生成速度 | K2 57–60 · R1 67 · Qwen3.5 39 · Limite 32 tok/s |
| 「2024 有幾多個正因數」，All four，512 tok | K2 `\boxed{16}` ✓。R1 未諗完就俾截斷。Limite 截斷 |
| 同上，Auto，768 tok | 路由揀 Limite + K2，兩個都係 `\boxed{16}`，多數票通過，35 s |
| 廣東話問題，Auto | 路由去 Qwen3.5，答案流暢但概念唔準（0.8B 嘅上限） |
| 融合 7¹⁰⁰ mod 13，α = 0.5，product | 推理得出 9 ✓。top-1 一致率 86–88%，Limite 48.6% / R1 51.4% 主導，平均 JS 0.08 nats |
| 融合速度 | 約 10 tok/s：每步兩次 forward，Limite 48 層係瓶頸 |

### 14.4 成本

- ZeroGPU：hosting 免費，PRO 每日有 40 分鐘 GPU（2400 s）。部署當日測試用咗約 380 s，額外費用 **$0**。
- 超出額度先扣 credit，每 10 分鐘 $1。$10 預算可以買多約 100 分鐘。
- 未登入訪客每日只有 2 分鐘：一次 640-token 融合要預留約 92 s。

### 14.5 下一步

- 融合加速：兩個模型分開 CUDA stream 跑；或者用 AoTI 編譯 Limite decode。
- 4 GB 真機對照：T4 small 開 1 小時（約 $0.40），量 K2 Q4_K_M + R1 Q8_0 嘅實際速度。
- 學生 embedding head（§13.1）訓練完之後放上 Space，同 LFM 並排比較 8-NN 一致率。

---

## 15. 第六版設計文檔：逐模型對 sparse 權重合併（2026-11-26，設計、未執行）

§12.2 判咗「cross-teacher 權重合併 ❌」（四個 teacher 嘅形狀
1280×48 / 1536×28 / 1536×28 / 1024×24，Limite / R1 / K2 / Qwen3.5
逐對都唔啱位）、但「同基座 expert LoRA ✅」。本節把剩低嘅問題
——**邊個模型對、用邊層方法**——收成一張可行性矩陣 + 方法分層。
呢份係設計文檔（按第八輪後嘅決策），唔含可執行定理；
對齊後點融合嘅定律喺 [`DEFERRED.md`](DEFERRED.md) 第 3 節（T92–T94 / A73）。

### 15.1 逐模型對可行性矩陣

| 模型對 | 權重直併 | 可行方法 | 依據 |
|---|---|---|---|
| R1 ⊕ Limite | ❌（且詞表 151,643 **id 全同**，唔需要對齊） | output-KD / token-KD：兩個 teacher 對同一學生（R1-1.5B，§12.3）做軟目標 | 詞表逐位相同，融合只係分佈層面 |
| R1 ⊕ Qwen2.5-Coder | ❌ **W′ 高危**：rope base 1e6 vs 1e4、tie-breaking 規則不同 | 只做 output-KD；任何「搬 embedding/attention 欄位」嘅做法都唔接受 | 結構超參唔相容，搬咗會靜默壞 |
| K2 ⊕ R1 | ❌ 跨詞表 | **byte-align（§11.1 A40）**：交集 K2∩R1 = 42,986 token | `HEXAD_VOCAB_OVERLAP.txt` 實測 |
| K2 ⊕ Qwen3.5 | ❌ 跨詞表 | byte-align：K2∩Qwen3.5 = 50,932 | 同上 |
| R1 ⊕ Qwen3.5 | ❌ 跨詞表 | byte-align：R1∩Qwen3.5 = 131,586（最大交集） | 同上 |
| LFM ⊕ 任何 | ❌ 佢唔係 decoder 權重 | 只有 **圖 / 向量** 層面：k-NN 圖拉普拉斯（A39）、向量蒸餾（§13.1） | LFM 係 embedding 模型，冇 transformer 欄位可合 |

「邊個對都係 ❌ 直併」唔等於「乜都做唔到」：**權重層嘅正路係
同基座 expert 合併**——所有 teacher 先對齊到同一學生基座
（R1-1.5B，vocab 151,643）嘅 LoRA expert，再喺 expert 層做
TIES（trim 符號衝突、elect、merge）或 DARE（drop 0.5、
rescale `1/(1−p)`）。呢條路 §12.2 已判 ✅，本節只補「邊對 teacher
值得做」嘅排序：交集大嘅對（R1⊕Qwen3.5 131,586）byte-align 收益最高；
id 全同嘅對（R1⊕Limite）連對齊都唔使，直接 output-KD。

### 15.2 方法分層（由平貴到貴）

1. **output-KD / token-KD**（零結構風險）：只動軟目標，
   適用所有矩陣行；成本係每個 teacher 一次前向。
2. **byte-align token 級融合（A40，§11.1）**：跨詞表對先做，
   交集 token 上 logit 凸組合；對齊律嘅融合層見 DEFERRED.md T92–T94。
3. **同基座 expert LoRA + TIES/DARE**：唯一可以「減疊項」嘅權重層做法
   （§12.2 解讀二）；必須同一學生基座，唔同基座一律唔做。
4. **LFM 圖 / 向量**：獨立通道，唔與上面三層混。

### 15.3 誠實邊界

* 沙盒**冇 GPU、冇權重**：本節係設計文檔，一個 merge 都冇執行過；
  所有交集數字來自 [`HEXAD_VOCAB_OVERLAP.txt`](HEXAD_VOCAB_OVERLAP.txt) 嘅靜態量度，
  形狀 / rope / tie 診斷來自各 model card。
* A40 只喺 §11.1 實測過對齊度（K2⊕R1）；其餘對嘅 byte-align 係
  同一演算法、未跑過端到端。
* TIES/DARE 參數（drop 0.5、rescale `1/(1−p)`）係文獻預設，
  未喺本專案嘅任何模型上調過。

---

## 16. 第七版設計文檔：L3 蒸餾執行計劃（2026-11-26，設計、未執行）

§5.4 判咗 L3「可選、最後做」、§8 釘咗里程碑 M7、§12.3 選咗學生基座、
§13.1 給咗向量頭三段損失、§15 給咗逐模型對合併——蒸餾嘅材料散喺五處。
本節將佢哋收攏成一份**可執行嘅蒸餾計劃**：輸入、三條通道、資料過濾器、
驗收憑證、依賴成本。同 §15 一樣係設計文檔：
沙盒冇 GPU、冇權重，**一個訓練步都冇跑過**；對齊層嘅可執行定律見
[`DEFERRED.md`](DEFERRED.md) 第 3 節（T92–T94 / A73）。

### 16.1 學生基座與輸入（收斂後嘅決定）

* **學生**：R1-Distill-Qwen-1.5B（§12.3 選型；vocab 151,643 = Limite id 全同
  ⇒ 同詞表 token-KD 免對齊）。§5.4 嘅 L3-α（Qwen3.5-0.8B）/ L3-β（K2-0.9B）
  記錄為備選——**換基座 = §15 全張矩陣重判**（vocab / rope / tie 都會變），唔係免費切換。
* **教師信號三種**（按通道取）：
  1. L2 融合分佈：同詞表 α-混合；跨詞表先 A40 byte-align（§11.1），
     對齊後融合層定律 = DEFERRED T92–T94。
  2. 各 teacher 嘅 output logits（token-KD）：R1⊕Limite 直接混合；
     K2 / Qwen3.5 對走 byte-align 交集（42,986 / 50,932 / 131,586，§15）。
  3. LFM2.5-Embedding-350M 向量（1024-d）：int8 快照已喺
     `space/data/agent_vecs_int8.npy`，唔使再推理。
* **資料**：L1 題庫軌跡，只留通過 ℓ≥2 分級反射驗證嘅（§5.4）；
  frame1 圖拉普拉斯共識（A39）對教師衝突樣本降權——
  「蒸餾數據過濾器 + 自學裁判」係 §12.6 講嘅組合新意，呢度係佢嘅執行位。

### 16.2 三條通道（按 §12.2 分層，全部同基座）

| 通道 | 內容 | 損失 / 方法 | 出處 |
|---|---|---|---|
| (a) KD 主力 | 序列級 SFT 為主；同詞表加 token 級 KD；跨詞表 universal logit distillation 作實驗項 | `τ²·KL(softmax(z_T/τ) ‖ softmax(z_S/τ))`，τ ∈ {1,2,4} 掃描 | §5.4、Hinton 2015 形式 |
| (b) 向量頭 | 學生加 1536→1024 投影頭（≈1.57M 參數）、last-token pooling、**骨幹凍結** | 三段：cosine 對齊 LFM ＋ batch 相似度矩陣關係蒸餾（**保 k-NN 圖結構**）＋ InfoNCE；中文用中英平行句 | §13.1 原文 |
| (c) 專家合併 | 異詞表教師分別蒸成 R1-1.5B 嘅 expert LoRA | TIES / DARE（drop 0.5、rescale `1/(1−p)`）——只 §15 方法分層第 3 層 | §12.2、§15 |

三條通道互相獨立、可分別消融（M6 嘅對照組設計）；**唔做**跨基座欄位搬運、
**唔做**未過濾分佈直蒸（§12.6 兩個否決）。

### 16.3 驗收憑證（執行前釘死，唔事後凑）

| 憑證 | 閾值 | 狀態 |
|---|---|---|
| 學生向量圖 vs LFM 圖 8-NN 一致率 | **≥ 0.90** | §13.1 已釘 |
| 同上 λ₂ 差距 | **≤ 10%** | §13.1 已釘 |
| KD 收斂 | 驗證集 KL 隨 step 單調下降（τ 各一條曲線） | 本節設計值 |
| top-1 一致率 | ≥ 最好單教師；對照 M6 全量基線 | §8 M7「對比 M6」 |
| 合併回落 | TIES/DARE 合併後 KL 唔差過合併前 5% | 本節設計值（§12.2 解讀二） |
| 資源 | 4 GB VRAM 內完成；tok/s 對照 §14.3 | §12.5 |

### 16.4 依賴與成本

* **GPU 必須**：KD 訓練唔係推理——ZeroGPU PRO 每日 40 min（§14.4）或
  T4 $0.40/h（§13.3）；§13.3 嘅 $0.03 CPU 係推理成本，唔適用於訓練。
* **前置已齊**：M6 基線（§14.3 實測）、A40 對齊度（K2⊕R1 已實測）、
  LFM 向量快照、frame1 wasm 驗證器（幾百 KB，保留作推理閘門，§5.4）。
* **產出**（M7）：一個學生 checkpoint + model card，對比 M6。

### 16.5 誠實邊界

* 設計文檔：零訓練步驟。閾值分兩類——「已釘」（0.90 / 10%，§13.1 有出處）
  同「設計值」（KL 單調、5% 回落，執行時先知啱唔啱），冇混充成實測。
* 引用數字全部有出處：151,643 / 1.9 GB（§12.3）、1.57M 參數 / 99.7%（§13.1）、
  1,275 + 489 條邊（A39）、交集數字（§15）、成本（§13.3 / §14.4）。
* `τ²·KL` 形式係 KD 文獻標準，未喺本專案任何模型上調過；
  §13.1 嘅 8-NN ≥ 0.9 係計劃閾值，對照組係已實測嘅 LFM 自身 99.7%。
