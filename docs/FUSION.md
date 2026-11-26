# 同詞表 token 融合（T81–T83 / A65–A67）

> 對象：Hugging Face Space [`Transynergy/frame1-hexad`](https://huggingface.co/spaces/Transynergy/frame1-hexad)
> 第 4 個 tab「🔁 同詞表融合 × 3 循環」，以及 MoonBit 認證模組 `src/39_shared_vocab_fusion.mbt`。

## 0. 一句話

`Limite-1B-Violetto`、`Qwen2.5-Coder-1.5B`、`Qwen3-0.6B` 與同家族的 `DeepSeek-R1-Distill-Qwen-1.5B`、
`Qwen3.5-0.8B` **共用同一個 151 643 個 regular id 的 Qwen2 BPE 詞表**，所以可以在
**同一個標籤空間**上逐步投票：每一步只產生**一個共同 token**，跑完一輪把結果重餵給所有成員，
**循環 3 次**。不是文本層 ensemble——三個模型從頭到尾在看同一串 id。

## 1. 三個對象（`space/fusion.py`）

| 對象 | 內容 |
|---|---|
| `SharedVocab` | 共同座標系：擁有者 tokenizer（Qwen2 家族任一）、`V = 151 643`、兩個偽 token：`END = V`、`THINK_CLOSE = V+1`、寬度 `V+2`、alias 合併表 |
| `FusionMember` | 一個成員：`tokenizer / model / prompt_suffix / end_ids / think_close_ids / weight / same_vocab`，執行期帶 KV cache 與最後一格 logits |
| `run_cycles(cycles=3)` | 三循環主控：第 *n* 輪用「第 *n−1* 輪的共同 id 序列」重錨（`bootstrap="own"`），`seed + n` 做種，答案長度預算累加 |

## 2. 三條融合律（T81）

```
product  (log-linear)   p ∝ Π_i p_i^{w_i}          Σ w_i = 1
mixture  (linear)       p = Σ_i w_i p_i
assisted (speculative)  p ∝ [p_s − p_w]_+ + p_w
```

三條都寫成 `log p` 的向量運算，`−∞` 代表硬零。權重一律 `fusion_weights(n, weight)`：
均分，或 `w₀ = weight`、其餘 `(1 − weight)/(n − 1)`。

**T81（乘積律的等冪與硬否決）**
1. 若所有成員給同一個分佈 `p`，則 `P(p, …, p) = p`（因為 `Σ w_i = 1`）。
2. 若某成員在 id `x` 給 0，則 `P` 在 `x` 給 0——**乘積律是支持集的交集**。

> ⇒ 這就是為什麼融合必須在**共同詞表**上做：任何一個 id 只要不在某成員的詞表裡，
> 在乘積律下就變成永久禁用的 token（`fuse_vocab_coverage` 量度「至少一個成員可發」與
> 「所有成員都可發」的 id 數）。操作上的解法是先把 alias id 併進同一個共同槽（`SharedVocab.aliases`）。

**T82（JS 散度對遮罩安全）** `JS = (1/k)Σ D(p_i‖p̄)`，`p = 0` 的格位貢獻**恰為 0**；
即使某成員在整條 support 上為 0，JS 仍然有限。實作以 `torch.where(p > 0, …)` 遮罩，
而不是靠 IEEE 的 `0·log 0 = nan`。

**T83（循環可結合）** 把每輪的 token 流串接就得到多輪流；預算可加；取樣種子
`(seed + cycle) mod m` 只依賴 (seed, 循環數)，與 wall-clock 無關 ⇒ 3 循環結果可重現。

## 3. 兩份實作、同一組恆等式

| 實作 | 檔案 | 驗證方式 | 現況 |
|---|---|---|---|
| Python（GPU 端） | `space/fusion.py` | `python3 space/test_fusion.py` | **13/13 通過**（CPU，假 tokenizer + 假模型，無需下載權重） |
| Python 冒煙 | `space/demo_fusion_cpu.py` | 三個合成分佈走完 `run_cycles(3)`，輸出 `space/fusion_demo.csv` | 通過（**數字是合成的**，只證明管線／診斷正確） |
| MoonBit（認證端） | `src/39_shared_vocab_fusion.mbt` | `moon test --target wasm|js` | **146/146 通過**（兩 target，含 T81–T83），`moon fmt --check` 乾淨 |

MoonBit 端刻意用 `Array[Double]` + 手寫 `log_softmax`，不碰 GPU：**數值律與模型伺服器解耦**，
所以 GPU 側的融合只要有疑問，跑一句 `moon test` 就知道是律錯還是接線錯。

### 本輪真正抓到的三個 bug（都已修）

1. **乘積律的權重是指數，不是加項**。誤寫成 `log p_i + log w_i` 時，等權（`Σ w = 1`）會退化成
   「每個成員指數 1」，分佈被過度尖銳化。T81 的等冪性（三個相同成員必須原樣返回）一測即中；
   Python 側同步修正並加回歸測試 `test_product_weights_are_exponents`。
2. **硬否決被「跳過」吃掉**。把 `−∞` 格位單純 `continue` 掉，其他成員仍會貢獻質量 ⇒ 不再是交集。
   改為獨立 `veto[]` 遮罩，最後統一套用。
3. **`mixture` 的 0 權重殘留**。`+log(max(w, 1e-12))` 令 `w = 0` 的 `−∞` 格變成 `−27.6`；
   現在 `w = 0` 的成員在進入 `logsumexp` 前就被剔除。
4. （附帶）`fuse_cycle_seed` 在 `Int32` 下 `2147483646 + 2` 會溢出成負數 ⇒ 改用 `Int64` 做加法。

## 4. UI（app.py 第 4 個 tab）

| 控制項 | 作用 |
|---|---|
| 成員 checkbox | `Limite / Qwen2.5-Coder-1.5B / Qwen3-0.6B` 預設全選；另列 `DeepSeek-R1-Distill-Qwen-1.5B`、`Qwen3.5-0.8B`，**剔選任何組合**都走同一條管線 |
| 融合律 radio | `product / mixture / assisted` |
| 循環 slider | 1–5，**預設 3** |
| 權重 slider | 0 = 均分；> 0 = 第一個成員佔比 |
| 巡迴模式 checkbox | 把勾選名單的**每一個子集**（單體 / 全部配對 / 三體）逐個跑 3 循環，輸出可比較的摘要表 |
| 進階 | 每輪 token 預算（32–768）、temperature、top-p、seed、重餵策略 `own/none`、75 % 預算強制關閉 thinking |

輸出 5 個：答案 markdown（含 `\boxed{}`）、串接的 thinking、診斷 DataFrame、
每輪貢獻 DataFrame、**逐 token trace CSV**（`cycle, step, fused_id, fused_token, picks, top1_agree, js`；
巡迴模式再加 `combo` 欄）。

「其餘都係同詞 token 融合」由 `fusion.plan(roster, max_size)` 兌現：3 個成員 ⇒ 7 個組合
（3 單體 + 3 配對 + 1 三體）；4 個成員且 `max_size=2` ⇒ 10 個組合（6 個配對全覆蓋）。
所有組合共用同一 prompt / seed / 融合律，所以摘要表可直接比較。

## 5. 誠實邊界

- **沙盒冇 gradio / transformers / GPU**：`app.py` 只做過 `ast.parse` 靜態檢查 + 靜態接線測試
  （`test_static_app_wiring`），從未真機跑過。真機首次運行若失敗，最可能是模型 id／revision 或
  `past_key_values` 介面差異，而不是融合數學。
- 三成員 bf16 約 4 GB VRAM；巡迴模式是**順序**載入，14 個組合的 VRAM 不會同時佔用，但耗時是
  單次運行的十幾倍（ZeroGPU 有時間配額）。
- `assisted` 律在成員分佈完全相同時等於 `p_s`；不相同時峰值跟隨 strong 成員，
  但這**不是**「強模型必然更好」的主張，只是拒絕取樣的代數結果。
- T81 的「交集」結論只在 `weights > 0` 對所有成員成立時有效；`weight=1` 退化為單成員。
