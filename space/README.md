---
title: frame1-Hexad
emoji: 🧬
colorFrom: purple
colorTo: blue
sdk: gradio
sdk_version: 6.32.0
app_file: app.py
python_version: "3.12"
startup_duration_timeout: 1h
short_description: 5 small LLMs as one population - route, fuse, couple
license: apache-2.0
models:
  - IFM/K2-Horizon-0.9B
  - paradigma-inc/limite-1b-violetto
  - LiquidAI/LFM2.5-Embedding-350M
  - Qwen/Qwen3.5-0.8B
  - Qwen/Qwen2.5-Coder-1.5B-Instruct
  - Qwen/Qwen3-0.6B
  - deepseek-ai/DeepSeek-R1-Distill-Qwen-1.5B
tags:
  - multi-agent
  - ensemble
  - logit-fusion
  - embeddings
  - zerogpu
---

# frame1-Hexad

Five small open models run as one population, on top of the
[frame1](https://github.com/vv0101v01/frame1) emergence framework (MoonBit, commit `1545a42`).

| Model | Role in the Hexad |
|---|---|
| K2-Horizon-0.9B | general reasoning, planning, code and STEM |
| Limite-1B-Violetto | competition maths |
| DeepSeek-R1-Distill-Qwen-1.5B | long chain-of-thought |
| Qwen3.5-0.8B | fast multilingual chat |
| LFM2.5-Embedding-350M | router, consensus judge, agent coupling graph |

**Tabs**

1. **Router + consensus.** LFM embeds the query and finds the nearest of frame1's 128 agent role cards. Their votes decide which generators answer. Consensus is decided first by a majority on `\boxed{}` answers. Otherwise the answer most similar to the others (LFM cosine) wins.
2. **Limite ⊕ R1 fusion.** The two models share 151,643 Qwen2 token ids, so they are decoded together at the token level. Each step samples from the product (log-linear) or the mixture of the two distributions, with α weighting Limite.
3. **frame1 coupling graph.** Port of frame1 module 25: symmetric k-NN on the LFM vectors, heat-kernel weights, cross-family bridges and a Kruskal repair pass. Compared with frame1's ring, λ₂ rises from 0.0024 to 0.540.
4. **🔁 同詞表融合 × 3 循環 (shared-vocabulary fusion, 3 cycles).** 三個（或更多）同屬 Qwen2 BPE 家族的模型——
   `Limite-1B-Violetto` + `Qwen2.5-Coder-1.5B` + `Qwen3-0.6B`——在**同一個 151,643-id 標籤空間**上逐步投票：
   每步由 `product` / `mixture` / `assisted` 三種融合律之一合成一個分佈，**只採樣一個共同 token**；
   跑完一輪後把融合結果重餵給所有成員再跑下一輪，**預設循環 3 次**（UI 可調 1–5）。
   勾選「巡迴模式」則把名單的**每一個**子集（單體 / 全部配對 / 三體）都用同一條管線跑 3 循環，
   結果可直接比較（同一 prompt / seed / 融合律）。輸出：答案、串接的 thinking、診斷表、
   每輪貢獻表、以及**逐 token trace CSV**（`cycle, step, fused_id, fused_token, picks, top1_agree, js`）。
   數值核心 `fusion.py` 完全 GPU-free，可用 `python3 test_fusion.py`（13 個測試）與
   `python3 demo_fusion_cpu.py`（假模型冒煙測試）在 CPU 上驗證；MoonBit 端 `src/39_shared_vocab_fusion.mbt`
   把三條融合律（T81–T83）在 wasm/js 兩個 target 上一併釘死，兩邊共用同一組恆等式。
5. **Diagnostics.** Shows the load status of each model. It also checks that the HF-transformers LFM encoder reproduces the llama.cpp Q8_0 vectors frame1 was built from.

`data/` holds the 128 role cards, their int8 LFM vectors (frame1 module 26) and the frame1 MoonBit trace for 128 agents × 640 steps.

Runs on ZeroGPU. All five models are loaded in bf16, about 10 GB of VRAM in total.
The fusion tab lazily loads the three Qwen-family members on first use (bf16, ~4 GB for the trio);
each extra member adds VRAM, and the sweep mode runs the subsets **sequentially**, so a full
sweep of five models is 14 model-loads' worth of work rather than one.

**Fusion, in one paragraph.** `p ∝ Π p_i^{w_i}` (product, hard veto: an id that any member puts at
probability 0 stays at 0) · `p = Σ w_i p_i` (mixture) · `p ∝ [p_s − p_w]_+ + p_w` (assisted /
speculative). One cycle = one shared decode of up to *N* tokens; cycles chain through the shared id
stream, so a 3-cycle run is literally the concatenation of three streams — certified by T83 in
`src/39_shared_vocab_fusion.mbt`, and enforced by `test_fusion.py` on CPU.
