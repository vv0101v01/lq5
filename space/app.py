"""frame1-Hexad: five small models, one population.

K2-Horizon-0.9B · Limite-1B-Violetto · Qwen3.5-0.8B · DeepSeek-R1-Distill-Qwen-1.5B generate;
LFM2.5-Embedding-350M routes, scores consensus and builds the agent coupling graph (frame1 module 25).
"""
import spaces  # must be imported before torch / anything that touches CUDA

import json
import math
import os
import re
import tempfile
import time
import traceback

import gradio as gr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
import torch.nn.functional as F  # noqa: E402
import transformers  # noqa: E402
from transformers import AutoModelForCausalLM, AutoTokenizer  # noqa: E402

import graph  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DTYPE = torch.bfloat16

# --------------------------------------------------------------------------------------
# Model registry (revisions pinned: two repos ship custom modelling code)
# --------------------------------------------------------------------------------------
GEN_SPECS = {
    "K2-Horizon-0.9B": dict(
        repo="IFM/K2-Horizon-0.9B",
        rev="6b4f28e81b8a77c8b44465a93a452eea198d4ceb",
        trust=True,
        tmpl={"reasoning_effort": "low"},
        think_prefilled=True,  # template ends with <ifm|think_faster>
        temperature=1.0,
        top_p=0.95,
        role="general · planning · code · STEM",
    ),
    "Limite-1B-Violetto": dict(
        repo="paradigma-inc/limite-1b-violetto",
        rev="190ee8d9c6777b7c5f8e7592cc00caf495a90ca7",
        trust=True,
        tmpl={},
        think_prefilled=False,
        temperature=0.6,
        top_p=0.95,
        role="competition maths (fixed \\boxed{} system prompt)",
    ),
    "DeepSeek-R1-Distill-Qwen-1.5B": dict(
        repo="deepseek-ai/DeepSeek-R1-Distill-Qwen-1.5B",
        rev="ad9f0ae0864d7fbcd1cd905e3c6c5b069cc8b562",
        trust=False,
        tmpl={},
        think_prefilled=True,  # template ends with <think>\n
        temperature=0.6,
        top_p=0.95,
        role="long chain-of-thought reasoning",
    ),
    "Qwen3.5-0.8B": dict(
        repo="Qwen/Qwen3.5-0.8B",
        rev="2fc06364715b967f1860aea9cf38778875588b17",
        trust=False,
        tmpl={"enable_thinking": False},
        think_prefilled=False,
        temperature=0.7,
        top_p=0.8,
        role="fast multilingual chat (non-thinking)",
    ),
    # --- 同一 Qwen2 BPE 家族：常規 id 與 Limite 完全重合，故可與它做「同詞表 token 融合」
    "Qwen2.5-Coder-1.5B": dict(
        repo="Qwen/Qwen2.5-Coder-1.5B-Instruct",
        rev="main",  # 官方 instruct 版沒有 pinned hash 時用 main（同族 tokenizer 不受影響）
        trust=False,
        tmpl={},
        think_prefilled=False,
        temperature=0.3,
        top_p=0.95,
        role="code (Qwen2.5-Coder, shared 151 643-token vocabulary)",
    ),
    "Qwen3-0.6B": dict(
        repo="Qwen/Qwen3-0.6B",
        rev="main",
        trust=False,
        tmpl={"enable_thinking": True},
        think_prefilled=True,  # 模板以 <think>\n 收尾
        temperature=0.6,
        top_p=0.95,
        role="reasoning (Qwen3, shared 151 643-token vocabulary)",
    ),
}
GEN_NAMES = list(GEN_SPECS)
SHORT = {
    "K2-Horizon-0.9B": "K2",
    "Limite-1B-Violetto": "Limite",
    "DeepSeek-R1-Distill-Qwen-1.5B": "R1",
    "Qwen3.5-0.8B": "Qwen3.5",
    "Qwen2.5-Coder-1.5B": "Coder",
    "Qwen3-0.6B": "Qwen3",
}
LFM_REPO, LFM_REV = "LiquidAI/LFM2.5-Embedding-350M", "f35ae2c91d687658dbf1f2b449382f0b019b9808"

# frame1 family index (module 26): 0 K2, 1 Limite, 2 R1, 3 Qwen3.5, 4 LFM
FAMILY_TO_GEN = {0: "K2-Horizon-0.9B", 1: "Limite-1B-Violetto", 2: "DeepSeek-R1-Distill-Qwen-1.5B", 3: "Qwen3.5-0.8B"}

# 同詞表融合成員（Round 7 · 通道十一）：全部落在 Qwen2 家族的 151 643 個常規 id 上
FUSION_NAMES = ["Limite-1B-Violetto", "Qwen2.5-Coder-1.5B", "Qwen3-0.6B"]
FUSION_DEFAULT = FUSION_NAMES[:]          # 預設三個一起融合
FUSION_EXTRA = ["DeepSeek-R1-Distill-Qwen-1.5B", "Qwen3.5-0.8B"]  # 其餘同家族模型亦可加入

LOAD_LOG: list[str] = []
GENS: dict[str, tuple] = {}
LIMIT_TOKENIZER = None  # 共同詞表的擁有者（Limite/Qwen2 家族），載入後設定


def _stop_ids(tok) -> list[int]:
    ids = set()
    if tok.eos_token_id is not None:
        ids.add(int(tok.eos_token_id))
    vocab = tok.get_vocab()
    for t in ("<|im_end|>", "<|ifm|im_end|>", "<|endoftext|>", "<｜end▁of▁sentence｜>"):
        if t in vocab:
            ids.add(int(vocab[t]))
    return sorted(ids)


for _name, _s in GEN_SPECS.items():
    _t0 = time.perf_counter()
    try:
        _tok = AutoTokenizer.from_pretrained(_s["repo"], revision=_s["rev"], trust_remote_code=_s["trust"])
        _kw = dict(revision=_s["rev"], trust_remote_code=_s["trust"], dtype=DTYPE)
        if _name == "Limite-1B-Violetto":
            _kw["attn_implementation"] = "sdpa"
        try:
            _model = AutoModelForCausalLM.from_pretrained(_s["repo"], **_kw)
        except Exception as _e1:  # e.g. a VLM checkpoint the text-only class cannot map
            from transformers import AutoModelForImageTextToText

            LOAD_LOG.append(f"⚠️ {_name}: CausalLM load failed ({type(_e1).__name__}), using ImageTextToText")
            _model = AutoModelForImageTextToText.from_pretrained(_s["repo"], **_kw)
        _model = _model.to("cuda").eval()
        GENS[_name] = (_tok, _model, _stop_ids(_tok))
        if _name == "Limite-1B-Violetto":
            LIMIT_TOKENIZER = _tok
        LOAD_LOG.append(f"✅ {_name}: {type(_model).__name__}, stop ids {_stop_ids(_tok)}, {time.perf_counter() - _t0:.0f}s")
    except Exception as _e:  # keep the Space up even if one architecture fails
        LOAD_LOG.append(f"❌ {_name}: {type(_e).__name__}: {_e}")
        traceback.print_exc()

# ---- 融合專用成員：Qwen2.5-Coder / Qwen3（延遲載入，第一次用到才下載，避免拖慢 Space 啟動）
FUSION_LAZY = {
    "Qwen2.5-Coder-1.5B": GEN_SPECS["Qwen2.5-Coder-1.5B"],
    "Qwen3-0.6B": GEN_SPECS["Qwen3-0.6B"],
}


def ensure_loaded(name: str):
    """回傳 (tok, model, stop_ids)；融合專用成員第一次呼叫時才載入。"""
    if name not in GENS:
        spec = FUSION_LAZY.get(name) or GEN_SPECS[name]
        t0 = time.perf_counter()
        tok = AutoTokenizer.from_pretrained(spec["repo"], revision=spec["rev"], trust_remote_code=spec["trust"])
        kw = dict(revision=spec["rev"], trust_remote_code=spec["trust"], dtype=DTYPE)
        try:
            model = AutoModelForCausalLM.from_pretrained(spec["repo"], **kw)
        except Exception:
            from transformers import AutoModelForImageTextToText

            model = AutoModelForImageTextToText.from_pretrained(spec["repo"], **kw)
        model = model.to("cuda").eval()
        GENS[name] = (tok, model, _stop_ids(tok))
        LOAD_LOG.append(f"✅ (lazy) {name}: loaded on demand in {time.perf_counter() - t0:.0f}s")
    return GENS[name]


# ---- LFM2.5-Embedding: native Lfm2Model + encoder patches (bidirectional attention, non-causal short conv).
# The repo's own remote patch targets transformers 4.56; this re-implements it against the current API.
LFM_TOK = LFM = None
try:
    from transformers.models.lfm2 import modeling_lfm2 as _lfm2_mod

    def _noncausal_conv_forward(self, hidden_states, *args, **kwargs):
        bcx = self.in_proj(hidden_states).transpose(-1, -2)
        b, c, x = bcx.chunk(3, dim=-2)
        bx = b * x
        k = self.conv.weight.shape[-1]
        y = F.conv1d(bx, self.conv.weight, self.conv.bias, padding=k // 2, groups=bx.shape[1])
        y = y[..., : bx.shape[-1]]
        return self.out_proj((c * y).transpose(-1, -2).contiguous())

    _lfm2_mod.create_causal_mask = lambda *a, **k: None  # only the lfm2 module is patched
    LFM_TOK = AutoTokenizer.from_pretrained(LFM_REPO, revision=LFM_REV)
    LFM = _lfm2_mod.Lfm2Model.from_pretrained(LFM_REPO, revision=LFM_REV, dtype=DTYPE, attn_implementation="sdpa")
    for _m in LFM.modules():
        if isinstance(_m, _lfm2_mod.Lfm2Attention):
            _m.is_causal = False
        if isinstance(_m, _lfm2_mod.Lfm2ShortConv):
            _m.forward = _noncausal_conv_forward.__get__(_m)
    LFM = LFM.to("cuda").eval()
    LOAD_LOG.append("✅ LFM2.5-Embedding-350M: Lfm2Model + bidirectional/non-causal patches, CLS pooling")
except Exception as _e:
    LOAD_LOG.append(f"❌ LFM2.5-Embedding-350M: {type(_e).__name__}: {_e}")
    traceback.print_exc()

# ---- frame1 agent population (128 role cards + real LFM vectors from llama.cpp Q8_0, int8)
AGENTS = json.load(open(os.path.join(HERE, "data", "agents.json"), encoding="utf-8"))
FAMILY_NAMES = AGENTS["family_names"]
AGENT_FAM = np.array([a["family"] for a in AGENTS["agents"]])
AGENT_TEXT = [a["text"] for a in AGENTS["agents"]]
AGENT_VEC = np.load(os.path.join(HERE, "data", "agent_vecs_int8.npy")).astype(np.float64)
AGENT_UNIT = AGENT_VEC / np.linalg.norm(AGENT_VEC, axis=1, keepdims=True)
TRACE = pd.read_csv(os.path.join(HERE, "data", "hexad_trace.csv"))

# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
CLOSE_RE = re.compile(r"</(?:ifm\|)?think(?:_fast|_faster)?>")
OPEN_RE = re.compile(r"<(?:ifm\|)?think(?:_fast|_faster)?>")
BOXED_RE = re.compile(r"\\boxed\{")


def last_boxed(text: str):
    """Content of the last \\boxed{...} (balanced braces)."""
    starts = [m.end() for m in BOXED_RE.finditer(text or "")]
    if not starts:
        return None
    i = starts[-1]
    depth, j = 1, i
    while j < len(text) and depth:
        depth += {"{": 1, "}": -1}.get(text[j], 0)
        j += 1
    return text[i : j - 1].strip() if depth == 0 else None


def norm_answer(a):
    if a is None:
        return None
    a = a.replace(" ", "").replace("\\!", "").replace("dfrac", "frac").replace("tfrac", "frac").rstrip(".")
    return a.strip("$")


def split_reasoning(raw: str, tok, prefilled: bool):
    for t in sorted(set(tok.all_special_tokens), key=len, reverse=True):
        if CLOSE_RE.fullmatch(t) or OPEN_RE.fullmatch(t):
            continue
        raw = raw.replace(t, "")
    m = None
    for m in CLOSE_RE.finditer(raw):
        pass
    if m is not None:
        reasoning, answer, truncated = raw[: m.start()], raw[m.end() :], False
    elif prefilled or OPEN_RE.search(raw[:40]):
        reasoning, answer, truncated = raw, "", True
    else:
        reasoning, answer, truncated = "", raw, False
    reasoning = OPEN_RE.sub("", reasoning).strip()
    return reasoning, answer.strip(), truncated


@torch.inference_mode()
def lfm_embed(texts):
    """LFM2.5-Embedding vectors (CLS pooling, L2-normalised), one text at a time (no padding)."""
    out = []
    for t in texts:
        ids = LFM_TOK(t, truncation=True, max_length=512, return_tensors="pt")["input_ids"]
        bos = LFM_TOK.bos_token_id
        if bos is not None and int(ids[0, 0]) != bos:
            ids = torch.cat([torch.tensor([[bos]]), ids], dim=1)
        h = LFM(input_ids=ids.to("cuda"), use_cache=False).last_hidden_state[0, 0].float()
        out.append(F.normalize(h, dim=0).cpu().numpy())
    return np.stack(out).astype(np.float64)


def route(qvec: np.ndarray, top_k: int = 16, temp: float = 0.1):
    s = AGENT_UNIT @ qvec
    top = np.argsort(-s)[:top_k]
    w = np.exp((s[top] - s[top].max()) / temp)
    w /= w.sum()
    fam = np.zeros(5)
    for i, wi in zip(top, w):
        fam[AGENT_FAM[i]] += wi
    gen = np.array([fam[f] for f in range(4)])
    gen = gen / gen.sum() if gen.sum() > 0 else np.full(4, 0.25)
    agents = pd.DataFrame(
        {
            "agent": [int(i) for i in top[:8]],
            "family": [FAMILY_NAMES[AGENT_FAM[i]] for i in top[:8]],
            "cosine": [round(float(s[i]), 4) for i in top[:8]],
            "role card": [AGENT_TEXT[i].removeprefix("document: ")[:110] for i in top[:8]],
        }
    )
    return {FAMILY_TO_GEN[f]: float(gen[f]) for f in range(4)}, agents


@torch.inference_mode()
def generate_one(name: str, prompt: str, max_new_tokens: int, seed: int):
    tok, model, stops = GENS[name]
    spec = GEN_SPECS[name]
    text = tok.apply_chat_template([{"role": "user", "content": prompt}], tokenize=False, add_generation_prompt=True, **spec["tmpl"])
    enc = tok(text, return_tensors="pt", add_special_tokens=False)
    ids = enc["input_ids"].to("cuda")
    torch.manual_seed(seed)
    t0 = time.perf_counter()
    out = model.generate(
        input_ids=ids,
        attention_mask=torch.ones_like(ids),
        max_new_tokens=int(max_new_tokens),
        do_sample=True,
        temperature=spec["temperature"],
        top_p=spec["top_p"],
        eos_token_id=stops,
        pad_token_id=stops[0],
    )
    new = out[0, ids.shape[1] :]
    dt = time.perf_counter() - t0
    raw = tok.decode(new, skip_special_tokens=False)
    reasoning, answer, truncated = split_reasoning(raw, tok, spec["think_prefilled"])
    cut = truncated or (int(new.numel()) >= int(max_new_tokens) and int(new[-1]) not in stops)
    return dict(name=name, reasoning=reasoning, answer=answer, truncated=truncated, cut=cut, tokens=int(new.numel()), seconds=dt)


def card_md(r: dict, weight=None) -> str:
    head = f"**{r['name']}** · {r['tokens']} tok · {r['seconds']:.1f}s · {r['tokens'] / max(r['seconds'], 1e-6):.0f} tok/s"
    if weight is not None:
        head += f" · router weight {weight:.2f}"
    box = last_boxed(r["answer"] or r["reasoning"])
    parts = [head]
    if box:
        parts.append(f"**boxed:** `{box}`")
    if r["answer"]:
        parts.append(r["answer"] + ("\n\n_(cut off at the token budget)_" if r.get("cut") else ""))
    elif r["truncated"]:
        parts.append("_(still thinking when the token budget ran out — raise max new tokens)_")
    if r["reasoning"]:
        rs = r["reasoning"]
        rs = rs if len(rs) < 6000 else rs[:2500] + "\n\n…\n\n" + rs[-2500:]
        parts.append(f"<details><summary>reasoning ({len(r['reasoning'])} chars)</summary>\n\n{rs}\n\n</details>")
    return "\n\n".join(parts)


def _router_duration(prompt, mode, manual, max_new_tokens, seed, *args, **kwargs):
    n = 4 if mode == "All four" else (max(1, len(manual or [])) if mode == "Manual" else 2)
    return int(min(170, 12 + n * (3 + int(max_new_tokens) * 0.045)))


@spaces.GPU(duration=_router_duration)
def hexad_router(prompt, mode="Auto (LFM router)", manual=None, max_new_tokens=768, seed=0):
    """Route a prompt through the frame1-Hexad population and return the consensus answer."""
    if not prompt or not prompt.strip():
        raise gr.Error("Enter a prompt.")
    if LFM is None:
        raise gr.Error("LFM2.5-Embedding failed to load — see Diagnostics.")
    t_all = time.perf_counter()
    weights, agents = route(lfm_embed(["query: " + prompt])[0])
    if mode == "All four":
        chosen = [n for n in GEN_NAMES if n in GENS]
    elif mode == "Manual":
        chosen = [n for n in (manual or []) if n in GENS]
    else:
        ranked = sorted(weights, key=weights.get, reverse=True)
        chosen = [n for n in ranked if n in GENS][:2]  # consensus needs at least two voices
    if not chosen:
        raise gr.Error("No generator selected / loaded.")
    wdf = pd.DataFrame({"model": [SHORT[n] for n in GEN_NAMES], "weight": [round(weights[n], 4) for n in GEN_NAMES]})
    cards = {n: ("_(not selected)_" if n not in chosen else "⏳ generating…") for n in GEN_NAMES}
    status = f"Router → **{', '.join(SHORT[n] for n in chosen)}** (LFM top-16 agent vote)"
    yield status, wdf, agents, *[cards[n] for n in GEN_NAMES], None

    results = []
    for n in chosen:
        r = generate_one(n, prompt, max_new_tokens, int(seed))
        results.append(r)
        cards[n] = card_md(r, weights[n])
        yield status + f"\n\n✅ {SHORT[n]} done ({r['seconds']:.1f}s)", wdf, agents, *[cards[x] for x in GEN_NAMES], None

    # ---- consensus: (1) majority on \boxed{} answers (verifiable), (2) LFM agreement medoid
    texts = [(r["answer"] or r["reasoning"][-1500:] or "(empty)") for r in results]
    vecs = lfm_embed(["document: " + t[:2000] for t in texts])
    sim = vecs @ vecs.T
    rw = np.array([weights[r["name"]] for r in results]) + 0.05
    if len(results) > 1:
        agree = np.array([sum(sim[i, j] * rw[j] for j in range(len(results)) if j != i) / sum(rw[j] for j in range(len(results)) if j != i) for i in range(len(results))])
    else:
        agree = np.ones(1)
    boxes = [norm_answer(last_boxed(r["answer"] or r["reasoning"])) for r in results]
    votes = {}
    for i, b in enumerate(boxes):
        if b:
            votes.setdefault(b, []).append(i)
    best_box = max(votes.items(), key=lambda kv: (len(kv[1]), sum(agree[i] for i in kv[1])), default=(None, []))
    complete = [i for i, r in enumerate(results) if r["answer"] and not r["cut"]]
    pool = complete or list(range(len(results)))
    if best_box[0] is not None and len(best_box[1]) >= 2:
        cand = [i for i in best_box[1] if i in pool] or best_box[1]
        pick = max(cand, key=lambda i: agree[i])
        why = f"majority on `\\boxed{{{best_box[0]}}}` ({len(best_box[1])}/{len(results)} models)"
    elif len(results) > 1:
        pick = max(pool, key=lambda i: agree[i] + 0.05 * rw[i])
        why = "LFM agreement medoid among complete answers (cosine to the other answers)"
    else:
        pick, why = 0, "single model"
    simdf = pd.DataFrame(np.round(sim, 3), columns=[SHORT[r["name"]] for r in results])
    simdf.insert(0, "", [SHORT[r["name"]] for r in results])
    simdf["agreement"] = np.round(agree, 3)
    win = results[pick]
    summary = (
        f"### 🧬 Hexad answer — {SHORT[win['name']]}\n"
        f"_picked by {why}; total {time.perf_counter() - t_all:.1f}s_\n\n"
        + (win["answer"] or "_(winner was cut off while thinking; see its card)_")
    )
    yield summary, wdf, agents, *[cards[x] for x in GEN_NAMES], simdf


# --------------------------------------------------------------------------------------
# Limite ⊕ R1 token-level logit fusion (shared 151 643-token Qwen2 vocabulary)
# --------------------------------------------------------------------------------------
SHARED_V = 151643  # regular BPE ids identical in both tokenizers
END, THINK_CLOSE = SHARED_V, SHARED_V + 1  # two pseudo-tokens appended to the fused vocabulary


def _fusion_duration(prompt, alpha, rule, max_new_tokens, *args, **kwargs):
    return int(min(175, 15 + int(max_new_tokens) * 0.12))


def _top_p_sample(logp: torch.Tensor, temperature: float, top_p: float, gen: torch.Generator) -> int:
    if temperature <= 0:
        return int(torch.argmax(logp))
    lv, li = torch.topk(logp, 512)  # sorted descending; mass outside the top-512 is negligible for top-p
    sp = torch.softmax(lv / temperature, dim=-1)
    sp = sp * ((torch.cumsum(sp, 0) - sp) < top_p)
    return int(li[torch.multinomial(sp / sp.sum(), 1, generator=gen)])


@spaces.GPU(duration=_fusion_duration)
@torch.inference_mode()
def fuse_limite_r1(prompt, alpha=0.5, rule="product (log-linear)", max_new_tokens=640, temperature=0.6, top_p=0.95, seed=0, force_answer=True):
    """Decode with p ∝ p_Limite^α · p_R1^(1-α) (or the α-mixture) over the shared Qwen2 vocabulary."""
    a, b = "Limite-1B-Violetto", "DeepSeek-R1-Distill-Qwen-1.5B"
    if a not in GENS or b not in GENS:
        raise gr.Error("Limite or R1 failed to load — see Diagnostics.")
    if not prompt or not prompt.strip():
        raise gr.Error("Enter a prompt.")
    (tl, ml, _), (tr, mr, _) = GENS[a], GENS[b]
    msgs = [{"role": "user", "content": prompt}]
    pl = tl.apply_chat_template(msgs, tokenize=False, add_generation_prompt=True) + "<think>\n"
    pr = tr.apply_chat_template(msgs, tokenize=False, add_generation_prompt=True)  # ends with <think>\n
    il = tl(pl, return_tensors="pt", add_special_tokens=False)["input_ids"].to("cuda")
    ir = tr(pr, return_tensors="pt", add_special_tokens=False)["input_ids"].to("cuda")
    l_end = [151643, 151645]  # <|endoftext|>, <|im_end|>
    r_end = [151643]  # <｜end▁of▁sentence｜>
    r_close = tr.convert_tokens_to_ids("</think>")  # 151649 (special in R1)
    l_close_seq = tl("</think>", add_special_tokens=False)["input_ids"]  # plain text in Limite
    gen = torch.Generator(device="cuda").manual_seed(int(seed))

    out_l = ml(input_ids=il, use_cache=True)
    out_r = mr(input_ids=ir, use_cache=True)
    cache_l, cache_r = out_l.past_key_values, out_r.past_key_values
    logit_l, logit_r = out_l.logits[0, -1].float(), out_r.logits[0, -1].float()
    pieces = []
    st_agree = torch.zeros((), device="cuda"); st_lead_l = torch.zeros((), device="cuda"); st_js = torch.zeros((), device="cuda")
    n_steps, closed, forced = 0, False, False
    t0 = time.perf_counter()
    for step in range(int(max_new_tokens)):
        lpl, lpr = torch.log_softmax(logit_l[: len(tl)], -1), torch.log_softmax(logit_r[: len(tr)], -1)
        vl = torch.cat([lpl[:SHARED_V], torch.logsumexp(lpl[l_end], 0)[None], lpl[l_close_seq[0]][None]])
        vr = torch.cat([lpr[:SHARED_V], torch.logsumexp(lpr[r_end], 0)[None], lpr[r_close][None]])
        if rule.startswith("product"):
            fused = alpha * vl + (1 - alpha) * vr
        else:
            fused = torch.logaddexp(math.log(max(alpha, 1e-9)) + vl, math.log(max(1 - alpha, 1e-9)) + vr)
        fused = torch.log_softmax(fused, -1)
        if force_answer and not closed and step >= int(0.75 * int(max_new_tokens)):
            tok_id, forced = THINK_CLOSE, True  # budget forcing: stop thinking, leave room for the answer
        else:
            tok_id = _top_p_sample(fused, float(temperature), float(top_p), gen)
        closed = closed or tok_id == THINK_CLOSE
        # diagnostics
        n_steps += 1
        st_agree += (torch.argmax(vl) == torch.argmax(vr)).float()
        st_lead_l += (vl[tok_id] >= vr[tok_id]).float()
        el, er = vl.exp(), vr.exp()
        lpm = (0.5 * (el + er)).clamp_min(1e-30).log()
        st_js += 0.5 * (el * (vl - lpm)).sum() + 0.5 * (er * (vr - lpm)).sum()
        if tok_id == END:
            break
        if tok_id == THINK_CLOSE:
            feed_l, feed_r, text = l_close_seq, [r_close], "</think>"
        else:
            feed_l, feed_r, text = [tok_id], [tok_id], None
        pieces.append(text if text is not None else tok_id)
        out_l = ml(input_ids=torch.tensor([feed_l], device="cuda"), past_key_values=cache_l, use_cache=True)
        out_r = mr(input_ids=torch.tensor([feed_r], device="cuda"), past_key_values=cache_r, use_cache=True)
        cache_l, cache_r = out_l.past_key_values, out_r.past_key_values
        logit_l, logit_r = out_l.logits[0, -1].float(), out_r.logits[0, -1].float()
    dt = time.perf_counter() - t0
    # decode
    text, buf = [], []
    for p in pieces:
        if isinstance(p, int):
            buf.append(p)
        else:
            text.append(tl.decode(buf)); buf = []; text.append(p)
    text.append(tl.decode(buf))
    raw = "".join(text)
    reasoning, answer, truncated = split_reasoning(raw, tl, True)
    n = max(n_steps, 1)
    stats = dict(n=n_steps, agree=float(st_agree), lead_l=float(st_lead_l), js=float(st_js))
    stats["lead_r"] = n_steps - stats["lead_l"]
    sdf = pd.DataFrame(
        {
            "metric": ["tokens", "seconds", "tok/s", "top-1 agreement Limite vs R1", "tokens led by Limite", "tokens led by R1", "mean JS divergence (nats)", "α", "rule"],
            "value": [stats["n"], round(dt, 1), round(stats["n"] / max(dt, 1e-6), 1), f"{100 * stats['agree'] / n:.1f}%", f"{100 * stats['lead_l'] / n:.1f}%", f"{100 * stats['lead_r'] / n:.1f}%", round(stats["js"] / n, 4), alpha, rule],
        }
    )
    box = last_boxed(answer or reasoning)
    md = ("_(thinking was closed by budget forcing at 75% of the budget)_\n\n" if forced else "") + (f"**boxed:** `{box}`\n\n" if box else "") + (answer or "_(still thinking when the budget ran out)_")
    return md, reasoning, sdf


# --------------------------------------------------------------------------------------
# 同詞表 token 融合：Limite · Qwen2.5-Coder · Qwen3（+ 同家族其餘模型），循環 N 輪
# --------------------------------------------------------------------------------------
import fusion as FZ  # noqa: E402  (space/fusion.py — 純張量，CPU 可測)

FUSION_RULE_CHOICES = ["product (log-linear)", "mixture (linear)", "assisted (speculative)"]


def _fusion_member(name: str, tok, model, shared: FZ.SharedVocab) -> "FZ.FusionMember":
    """把某個成員的 tokenizer / 模型包成 fusion.FusionMember（同詞表路徑）。"""
    stop = _stop_ids(tok)
    vocab = tok.get_vocab()
    end_ids = sorted({i for i in stop if i < shared.V}) or [int(tok.eos_token_id)]
    think_close = [int(vocab[t]) for t in FZ.THINK_CLOSE_TOKENS if t in vocab]
    if name == "Limite-1B-Violetto":
        # Limite 的官方模板不含 <think>，舊管線手動補上；保持一致
        return FZ.FusionMember(
            name=name,
            tokenizer=tok,
            model=model,
            prompt_suffix="<think>\n",
            end_ids=end_ids,
            think_close_ids=think_close
            or [int(tok("</think>", add_special_tokens=False)["input_ids"][0])],
            same_vocab=True,
        )
    return FZ.FusionMember(
        name=name,
        tokenizer=tok,
        model=model,
        end_ids=end_ids,
        think_close_ids=think_close,
        same_vocab=True,  # Qwen2 家族：id 逐個對齊
    )


def _fusion_duration_multi(prompt, members, rule, weight, cycles, max_new_tokens, *args, **kwargs):
    n = max(len(members or []), 1)
    return int(min(175, 18 + n * 6 + int(cycles or 1) * int(max_new_tokens or 64) * 0.05))


@spaces.GPU(duration=_fusion_duration_multi)
@torch.inference_mode()
def fuse_shared_vocab(
    prompt,
    members=None,
    rule="product (log-linear)",
    weight=0.0,
    cycles=3,
    max_new_tokens=96,
    temperature=0.6,
    top_p=0.95,
    seed=0,
    bootstrap="own",
    force_answer=True,
    sweep=False,
    max_size=3,
):
    """**同詞表 token 級融合 + 循環 N 次**（預設 3 輪）。

    每個成員都在同一個 Qwen2 151 643-token 標籤空間投票，每一步只生成一個共同 token；
    一輪跑完後把融合結果重餵給所有成員再跑下一輪，`cycles=3` 就是「循環重複進行 3 次」。
    """
    members = members or FUSION_DEFAULT
    if not prompt or not prompt.strip():
        raise gr.Error("Enter a prompt.")
    if not members:
        raise gr.Error("Pick at least one model.")
    if sweep:  # 「其餘都係同詞 token 融合」：把 roster 的**每一個**子集都跑 3 輪
        return _fusion_sweep(prompt, members, cycle_opts=dict(
            rule=rule, weight=weight, cycles=int(cycles), max_new_tokens=int(max_new_tokens),
            temperature=temperature, top_p=top_p, seed=seed, bootstrap=bootstrap,
            force_answer=force_answer, max_size=int(max_size),
        ))
    owner_tok = LIMIT_TOKENIZER
    if owner_tok is None:  # Limite 載入失敗時，借用同家族任何一個 tokenizer 當共同座標系
        for n in members:
            try:
                owner_tok = ensure_loaded(n)[0]
                break
            except Exception:
                continue
    if owner_tok is None:
        raise gr.Error("No shared-vocabulary tokenizer available — see Diagnostics.")
    shared = FZ.SharedVocab(owner_tok, size=FZ.QWEN_FAMILY_SHARED_V)
    fm, names = [], []
    for n in members:
        try:
            tok, model, _ = ensure_loaded(n)
        except Exception as e:  # 單一成員失敗不拖垮整個融合
            LOAD_LOG.append(f"❌ fusion member {n}: {type(e).__name__}: {e}")
            continue
        fm.append(_fusion_member(n, tok, model, shared))
        names.append(n)
    if not fm:
        raise gr.Error("No fusion member could be loaded — see Diagnostics.")

    t0 = time.perf_counter()
    res = FZ.run_cycles(
        fm,
        shared,
        prompt,
        cycles=int(cycles),
        weight=(None if float(weight) <= 0 else float(weight)),
        rule=rule,
        temperature=float(temperature),
        top_p=float(top_p),
        seed=int(seed),
        budget=int(max_new_tokens),
        bootstrap=bootstrap,
        force_answer=bool(force_answer),
        device="cuda",
    )
    dt = time.perf_counter() - t0

    raw = res.text
    reasoning, answer, truncated = split_reasoning(raw, owner_tok, True)
    box = last_boxed(answer or reasoning)
    n = max(res.stats.steps, 1)
    rows = [
        ("members", " · ".join(names)),
        ("vocabulary", f"shared Qwen2 regular ids = {shared.V} (+2 pseudo: END, </think>)"),
        ("rule", FZ.canonical_rule(rule)),
        ("cycles", int(cycles)),
        ("tokens (all cycles)", res.stats.steps),
        ("seconds", round(dt, 1)),
        ("tok/s", round(res.stats.steps / max(dt, 1e-6), 1)),
        ("top-1 agreement between members", f"{100 * res.stats.agreement / n:.1f}%"),
        ("mean JS divergence (nats)", round(res.stats.js / n, 4)),
        ("shared-vocab coverage", f"{100 * res.stats.covered / n:.1f}%"),
        ("per-cycle tokens", " / ".join(str(c.stats.steps) for c in res.cycles)),
        ("weights", " · ".join(f"{nm}={w:.2f}" for nm, w in zip(names, res.weights))),
    ]
    sdf = pd.DataFrame({"metric": [r[0] for r in rows], "value": [r[1] for r in rows]})
    cyc = pd.DataFrame(
        [
            dict(cycle=i + 1, tokens=c.stats.steps, chars=len(c.text),
                 agree=round(c.stats.agreement / max(c.stats.steps, 1), 3),
                 js=round(c.stats.js / max(c.stats.steps, 1), 4))
            for i, c in enumerate(res.cycles)
        ]
    )
    md = (
        f"**循環 {int(cycles)} 輪 · {' + '.join(names)} · 同詞表 token 融合 ({FZ.canonical_rule(rule)})**\n\n"
        + (f"**boxed:** `{box}`\n\n" if box else "")
        + (answer or "_(no answer text — the whole budget went to thinking)_")
    )
    # 逐 token trace：每一輪、每一步、每個成員抽到的 id / top-1 / 該步 JS
    trace_rows = []
    for ci, c in enumerate(res.cycles):
        for st in getattr(c, "trace", []) or []:
            trace_rows.append(
                dict(cycle=ci + 1, step=st["step"], fused_id=st["fused_id"],
                     fused_token=owner_tok.decode([st["fused_id"]], skip_special_tokens=False),
                     picks=" · ".join(str(x) for x in st["picks"]),
                     top1_agree=st["agree"], js=round(st["js"], 5))
            )
    trace_path = None
    if trace_rows:
        trace_path = os.path.join(
            tempfile.gettempdir(),
            f"fusion_trace_{int(time.time())}_{int(cycles)}x{int(max_new_tokens)}.csv",
        )
        pd.DataFrame(trace_rows).to_csv(trace_path, index=False)
    return md, reasoning, sdf, cyc, trace_path


def _fusion_sweep(
    prompt,
    roster,
    cycle_opts,
):
    """把 roster 的**所有子集**（1..max_size 個成員）逐一走同一條同詞表融合管線。

    這兌現「其餘都係同詞 token 融合」：不是只有 Limite/Qwen2.5-Coder/Qwen3 三體，
    而是 roster 內任何配對、任何單體，都用同一套 SharedVocab + 3 循環解碼，
    彼此可比較（同一 prompt、同一 seed、同一融合律）。輸出：每個組合的摘要 DataFrame
    + 逐 token trace CSV + 合併的 markdown 表。
    """
    roster = list(roster or FUSION_DEFAULT + FUSION_EXTRA)
    if not prompt or not prompt.strip():
        raise gr.Error("Enter a prompt.")
    plan = FZ.fusion_plan(roster, max_size=cycle_opts.pop("max_size", 3))
    t0 = time.perf_counter()
    summary, trace_rows = [], []
    for item in plan:
        try:
            md, _reason, _stats, cyc, tpath = fuse_shared_vocab(prompt, item["members"], **cycle_opts)
        except Exception as e:  # 單一組合失敗（例如 OOM）不拖垮整條巡迴
            summary.append(dict(combo=item["label"], size=item["size"], members=len(item["members"]),
                                cycles=cycle_opts.get("cycles", 3), tokens=None, agree=None,
                                js=None, chars=None, note=f"{type(e).__name__}: {e}"))
            continue
        c = cyc.iloc[-1] if len(cyc) else {}
        summary.append(
            dict(
                combo=item["label"],
                size=item["size"],
                members=len(item["members"]),
                cycles=int(cycle_opts.get("cycles", 3)),
                tokens=int(cyc["tokens"].sum()) if len(cyc) else 0,
                agree=round(float(cyc["agree"].mean()), 3) if len(cyc) else None,
                js=round(float(cyc["js"].mean()), 4) if len(cyc) else None,
                chars=int(cyc["chars"].sum()) if len(cyc) else 0,
                note=str(c.get("tokens", "")),
            )
        )
        if tpath and os.path.exists(tpath):
            df = pd.read_csv(tpath)
            df.insert(0, "combo", item["label"])
            trace_rows.append(df)
    dt = time.perf_counter() - t0
    sdf = pd.DataFrame(summary)
    traces = pd.concat(trace_rows, ignore_index=True) if trace_rows else pd.DataFrame()
    tpath = None
    if len(traces):
        tpath = os.path.join(tempfile.gettempdir(), f"fusion_sweep_{int(time.time())}.csv")
        traces.to_csv(tpath, index=False)
    ok = int(sdf["tokens"].notna().sum()) if len(sdf) else 0
    md = (
        f"**同詞表融合巡迴：{ok}/{len(plan)} 個組合 × {int(cycle_opts.get('cycles', 3))} 循環**（"
        f"{round(dt, 1)} s · 同一 prompt / seed / 融合律，故結果可直接比較）\n\n"
        + "```\n" + sdf.to_string(index=False) + "\n```"
    )
    return md, "", sdf, traces.head(200), tpath


# --------------------------------------------------------------------------------------
# frame1 graph tab (CPU) + custom population (GPU for LFM)
# --------------------------------------------------------------------------------------
FAM_COLORS = ["#e4572e", "#7b2cbf", "#2a9d8f", "#f4a261", "#3a86ff"]


def _graph_figure(emb, fam, w, names, title):
    xy = graph.pca2(emb)
    fig, ax = plt.subplots(figsize=(7, 5.6), dpi=110)
    iu = np.argwhere(np.triu(w, 1) > 0)
    for i, j in iu:
        cross = fam[i] != fam[j]
        ax.plot(*zip(xy[i], xy[j]), color="#999" if not cross else "#111", lw=0.35 if not cross else 0.4, alpha=0.35 if not cross else 0.18, zorder=1)
    for f in sorted(set(fam.tolist())):
        m = fam == f
        ax.scatter(xy[m, 0], xy[m, 1], s=26, color=FAM_COLORS[f % len(FAM_COLORS)], label=names[f], zorder=2, edgecolor="white", lw=0.4)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, loc="best", frameon=False)
    ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    return fig


def _dynamics(w_ring, w_knn, emb, steps=200):
    x0 = graph.pca2(emb)
    rows = []
    for name, w in (("ring", w_ring), ("LFM k-NN", w_knn)):
        d = graph.diffuse(w, x0, steps)
        r = graph.kuramoto(w, steps)
        rows += [dict(step=t, value=float(d[t] / d[0]), series=f"{name} · concept dispersion") for t in range(steps + 1)]
        rows += [dict(step=t, value=float(r[t]), series=f"{name} · Kuramoto R") for t in range(steps + 1)]
    return pd.DataFrame(rows)


def _stats_df(w_ring, w_knn, fam):
    a, b = graph.graph_stats(w_ring, fam), graph.graph_stats(w_knn, fam)
    keys = ["agents", "edges", "cross_family_edges", "components", "lambda2", "lambda_max", "heat_CFL_dt"]
    fmt = lambda v: f"{v:.5g}" if isinstance(v, float) else str(v)  # noqa: E731
    df = pd.DataFrame({"metric": keys, "ring": [fmt(a[k]) for k in keys], "LFM k-NN": [fmt(b[k]) for k in keys]})
    return df, b["lambda2"] / max(a["lambda2"], 1e-12)


def frame1_graph(k=8, tau=0.25):
    """frame1-Hexad 128-agent coupling graph from the stored LFM vectors (no GPU)."""
    w_knn = graph.knn_graph(AGENT_VEC, AGENT_FAM, int(k), float(tau))
    w_ring = graph.ring_graph(len(AGENT_FAM))
    df, ratio = _stats_df(w_ring, w_knn, AGENT_FAM)
    fig = _graph_figure(AGENT_VEC, AGENT_FAM, w_knn, FAMILY_NAMES, f"128 agents · LFM k-NN (k={int(k)}, τ={tau}) · λ₂ ×{ratio:.0f} vs ring")
    return df, fig, _dynamics(w_ring, w_knn, AGENT_VEC)


@spaces.GPU(duration=40)
def custom_graph(lines, k=3, tau=0.25):
    """Embed your own agents with LFM2.5 ("family | description" per line) and build the coupling graph."""
    if LFM is None:
        raise gr.Error("LFM2.5-Embedding failed to load — see Diagnostics.")
    rows = [l.strip() for l in (lines or "").splitlines() if l.strip()][:96]
    if len(rows) < 4:
        raise gr.Error("Give at least 4 lines.")
    fams, texts = [], []
    for r in rows:
        f, _, t = r.partition("|") if "|" in r else ("agents", "", r)
        fams.append(f.strip()); texts.append(t.strip())
    names = list(dict.fromkeys(fams))
    fam = np.array([names.index(f) for f in fams])
    emb = lfm_embed(["document: " + t for t in texts])
    w_knn = graph.knn_graph(emb, fam, int(k), float(tau))
    w_ring = graph.ring_graph(len(rows))
    df, ratio = _stats_df(w_ring, w_knn, fam)
    fig = _graph_figure(emb, fam, w_knn, names, f"{len(rows)} agents · LFM k-NN (k={int(k)}, τ={tau}) · λ₂ ×{ratio:.1f} vs ring")
    return df, fig, _dynamics(w_ring, w_knn, emb, steps=120)


@spaces.GPU(duration=30)
def lfm_selfcheck(n=16):
    """Compare HF-transformers LFM vectors with the llama.cpp Q8_0 vectors frame1 was built from."""
    if LFM is None:
        raise gr.Error("LFM failed to load.")
    idx = np.linspace(0, len(AGENT_TEXT) - 1, int(n)).astype(int)
    v = lfm_embed([AGENT_TEXT[i] for i in idx])
    cos = np.sum(v * AGENT_UNIT[idx], axis=1)
    s_hf, s_cpp = v @ v.T, AGENT_UNIT[idx] @ AGENT_UNIT[idx].T
    iu = np.triu_indices(len(idx), 1)
    corr = float(np.corrcoef(s_hf[iu], s_cpp[iu])[0, 1])
    return (
        f"**LFM self-check** on {len(idx)} frame1 role cards: cos(HF, llama.cpp Q8_0) mean **{cos.mean():.4f}**, min {cos.min():.4f}; "
        f"similarity-matrix Pearson r = **{corr:.4f}**  (≥0.98 ⇒ the encoder patches reproduce the GGUF embedder)"
    )


def diagnostics_md():
    return "\n".join(
        [f"- {l}" for l in LOAD_LOG]
        + [f"- torch {torch.__version__} · transformers {transformers.__version__} · gradio {gr.__version__}"]
        + [f"- 同詞表融合成員可用：{', '.join(n for n in FUSION_DEFAULT if n in GENS)}（其餘首次使用時延遲載入）"]
    )


# --------------------------------------------------------------------------------------
# UI
# --------------------------------------------------------------------------------------
CSS = """
#col-container { max-width: 1180px; margin: 0 auto; }
.dark .gradio-container { color: var(--body-text-color); }
"""
INTRO = """# 🧬 frame1-Hexad — many small models, one population
**K2-Horizon-0.9B · Limite-1B-Violetto · DeepSeek-R1-Distill-Qwen-1.5B · Qwen3.5-0.8B** generate,
**LFM2.5-Embedding-350M** routes, judges consensus and wires the 128-agent coupling graph of
[frame1](https://github.com/vv0101v01/frame1) (MoonBit)。五個細模型合成一個群體：路由 → 生成 → 共識 · LFM k-NN 耦合圖。

**🔁 同詞表 token 融合（Round 7 · 通道十一）**：Limite · Qwen2.5-Coder-1.5B · Qwen3
共用 151 643 個 Qwen2 常規 BPE id，因此可以逐一 token 投票解碼；一輪跑完把結果重餵給所有成員，
**循環重複 3 次**（次數可調 1–5）。R1 / Qwen3.5 亦屬同一家族，可隨時加入同一條融合管線。
"""

with gr.Blocks(title="frame1-Hexad") as demo:
    with gr.Column(elem_id="col-container"):
        gr.Markdown(INTRO)
        with gr.Tab("🧭 Router + consensus"):
            with gr.Row():
                prompt = gr.Textbox(show_label=False, placeholder="Ask anything — maths, code, planning, 中文…", container=False, scale=5, lines=2)
                run = gr.Button("Run Hexad", variant="primary", scale=1)
            with gr.Row():
                mode = gr.Radio(["Auto (LFM router)", "All four", "Manual"], value="Auto (LFM router)", label="Who answers")
                manual = gr.CheckboxGroup(GEN_NAMES, value=GEN_NAMES[:2], label="Manual selection")
            with gr.Accordion("Advanced", open=False):
                max_new = gr.Slider(64, 1536, value=768, step=32, label="max new tokens per model")
                seed = gr.Number(value=0, precision=0, label="seed")
            summary = gr.Markdown()
            with gr.Row():
                wplot = gr.BarPlot(x="model", y="weight", title="LFM router weights (top-16 agent vote)", height=220, scale=1)
                agents_df = gr.Dataframe(label="Nearest frame1 agents to the query", wrap=True, scale=2)
            cards = []
            with gr.Row():
                for n in GEN_NAMES[:2]:
                    with gr.Column():
                        cards.append(gr.Markdown(label=n, container=True))
            with gr.Row():
                for n in GEN_NAMES[2:]:
                    with gr.Column():
                        cards.append(gr.Markdown(label=n, container=True))
            sim_df = gr.Dataframe(label="Consensus: LFM cosine between answers")
            outs = [summary, wplot, agents_df, *cards, sim_df]
            run.click(hexad_router, [prompt, mode, manual, max_new, seed], outs)
            prompt.submit(hexad_router, [prompt, mode, manual, max_new, seed], outs)
            gr.Examples(
                examples=[
                    ["What is the sum of all positive integers n such that n^2 + 12n - 2007 is a perfect square? Put the final answer in \\boxed{}."],
                    ["Write a Python function that returns the longest palindromic substring of a string, and explain its complexity."],
                    ["Plan a 3-step experiment to test whether a k-NN coupling graph synchronises agents faster than a ring."],
                    ["用廣東話解釋乜嘢係代數連通度 (algebraic connectivity)。"],
                ],
                inputs=[prompt],
                outputs=outs,
                fn=hexad_router,
                cache_examples=True,
                cache_mode="lazy",
            )

        with gr.Tab("🔗 Limite ⊕ R1 fusion"):
            gr.Markdown(
                "Both models share the same 151 643 regular Qwen2 token ids, so they can be decoded **as one model**: "
                "each step samples from `p ∝ p_Limite^α · p_R1^(1−α)` (product) or `α·p_Limite + (1−α)·p_R1` (mixture). "
                "Special tokens are mapped to two pseudo-tokens (END, `</think>`)."
            )
            with gr.Row():
                fprompt = gr.Textbox(show_label=False, placeholder="A maths / reasoning question…", container=False, scale=5, lines=2)
                frun = gr.Button("Fuse", variant="primary", scale=1)
            with gr.Row():
                alpha = gr.Slider(0, 1, value=0.5, step=0.05, label="α (weight of Limite; 0 = pure R1, 1 = pure Limite)")
                rule = gr.Radio(["product (log-linear)", "mixture (linear)"], value="product (log-linear)", label="fusion rule")
            with gr.Accordion("Advanced", open=False):
                fmax = gr.Slider(64, 1280, value=640, step=32, label="max new tokens")
                ftemp = gr.Slider(0, 1.2, value=0.6, step=0.05, label="temperature (0 = greedy)")
                ftopp = gr.Slider(0.5, 1, value=0.95, step=0.01, label="top-p")
                fseed = gr.Number(value=0, precision=0, label="seed")
                fforce = gr.Checkbox(value=True, label="budget forcing: close </think> at 75% of the budget")
            fanswer = gr.Markdown()
            with gr.Row():
                freason = gr.Textbox(label="fused reasoning", lines=12, scale=3)
                fstats = gr.Dataframe(label="fusion diagnostics", scale=2)
            fouts = [fanswer, freason, fstats]
            finputs = [fprompt, alpha, rule, fmax, ftemp, ftopp, fseed, fforce]
            frun.click(fuse_limite_r1, finputs, fouts)
            gr.Examples(
                examples=[["If x + 3 = 8, what is x?"], ["How many positive divisors does 2024 have?"], ["Find the remainder when 7^100 is divided by 13."]],
                inputs=[fprompt],
                outputs=fouts,
                fn=fuse_limite_r1,
                cache_examples=True,
                cache_mode="lazy",
            )

        with gr.Tab("🔁 同詞表融合 × 3 循環"):
            gr.Markdown(
                "**Limite · Qwen2.5-Coder-1.5B · Qwen3** 三個模型共用同一套 151 643 個 Qwen2 常規 BPE id，"
                "所以可以在**同一個標籤空間**裡逐 token 投票：每步只生成一個共同 token，"
                "分佈為 `p ∝ Π p_i^{w_i}`（product）、`Σ w_i p_i`（mixture）、或拒絕採樣的 `assisted`。"
                "一輪跑完會把融合結果重餵給所有成員，**再跑下一輪**——預設循環 3 次。"
            )
            with gr.Row():
                uprompt = gr.Textbox(
                    show_label=False,
                    placeholder="A maths / code / reasoning question…（中英皆可）",
                    container=False,
                    scale=5,
                    lines=2,
                )
                urun = gr.Button("Fuse × 3", variant="primary", scale=1)
            with gr.Row():
                umembers = gr.CheckboxGroup(
                    FUSION_DEFAULT + FUSION_EXTRA,
                    value=FUSION_DEFAULT,
                    label="融合成員（同詞表；剔選任何一個都走同一條 token 融合管線）",
                )
                urule = gr.Radio(FUSION_RULE_CHOICES, value=FUSION_RULE_CHOICES[0], label="融合律")
            with gr.Row():
                ucycles = gr.Slider(1, 5, value=3, step=1, label="循環次數（預設 3）")
                uweight = gr.Slider(0, 1, value=0.0, step=0.05, label="權重（0 = 均分；>0 = 第一個成員佔比）")
            with gr.Accordion("Advanced", open=False):
                umax = gr.Slider(32, 768, value=96, step=32, label="每輪 max new tokens")
                utemp = gr.Slider(0, 1.2, value=0.6, step=0.05, label="temperature（0 = greedy）")
                utopp = gr.Slider(0.5, 1, value=0.95, step=0.01, label="top-p")
                useed = gr.Number(value=0, precision=0, label="seed")
                uboot = gr.Radio(
                    ["own", "none"],
                    value="own",
                    label="每輪重錨：own = 把前幾輪的融合結果重餵給所有成員；none = 每輪只看原始 prompt",
                )
                uforce = gr.Checkbox(value=True, label="預算強制：75% 預算時關閉 thinking")
                usweep = gr.Checkbox(
                    value=False,
                    label="巡迴模式：把勾選名單的**每一個**子集（單體 / 配對 / 三體）都用同詞表融合跑 3 循環",
                )
                umaxsize = gr.Slider(1, 4, value=3, step=1, label="巡迴最大成員數")
            uanswer = gr.Markdown()
            with gr.Row():
                ureason = gr.Textbox(label="融合後的 thinking（所有循環串接）", lines=14, scale=3)
                ustats = gr.Dataframe(label="融合診斷", scale=2)
            ucycle_df = gr.Dataframe(label="每輪貢獻（tokens / 字符 / 一致率 / JS）")
            utrace = gr.File(label="逐 token trace（CSV：每步 → 每輪 → 每成員的 top-1 / JS）")
            uouts = [uanswer, ureason, ustats, ucycle_df, utrace]
            uinputs = [
                uprompt, umembers, urule, uweight, ucycles, umax, utemp, utopp, useed, uboot,
                uforce, usweep, umaxsize,
            ]
            urun.click(fuse_shared_vocab, uinputs, uouts)
            uprompt.submit(fuse_shared_vocab, uinputs, uouts)
            gr.Examples(
                examples=[
                    ["Find the remainder when 7^100 is divided by 13. Put the answer in \\boxed{}."],
                    ["Write a Python function that returns the longest palindromic substring, and give its complexity."],
                    ["用廣東話解釋乜嘢係同詞表 token 融合 (shared-vocabulary token fusion)。"],
                ],
                inputs=[uprompt],
                outputs=uouts,
                fn=fuse_shared_vocab,
                cache_examples=False,
            )

        with gr.Tab("🕸️ frame1 coupling graph"):
            gr.Markdown(
                "frame1 Round 1 couples 128 agents on a **ring** (λ₂ = 0.0024). Hexad replaces it with a graph built from "
                "**real LFM2.5 vectors** of the agents' role cards: symmetric k-NN, heat-kernel weights `exp((cos−1)/τ)`, "
                "one bridge per agent to every other model family, Kruskal repair (frame1 module 25, ported to numpy — "
                "reproduces 1275 edges / λ₂ = 0.540 at k=8, τ=0.25)."
            )
            with gr.Row():
                gk = gr.Slider(1, 24, value=8, step=1, label="k")
                gtau = gr.Slider(0.05, 1.0, value=0.25, step=0.05, label="τ")
                gbtn = gr.Button("Rebuild 128-agent graph", variant="primary")
            with gr.Row():
                gstats = gr.Dataframe(label="ring vs LFM k-NN", scale=2)
                gfig = gr.Plot(label="agents (PCA of LFM vectors), black = cross-family edges", scale=3)
            gdyn = gr.LinePlot(x="step", y="value", color="series", title="numpy dynamics on both graphs (dispersion ↓, Kuramoto R ↑)", height=300)
            with gr.Accordion("frame1 MoonBit run: 128 agents × 640 steps (hexad_trace.csv)", open=False):
                tr = TRACE.melt(id_vars="step", var_name="series", value_name="value")
                gr.LinePlot(value=tr, x="step", y="value", color="series", height=320)
            with gr.Accordion("Your own population (LFM on GPU)", open=False):
                ctext = gr.Textbox(
                    lines=8,
                    label="one agent per line — `family | description`",
                    value="\n".join(
                        [
                            "math | proves inequalities and olympiad number theory",
                            "math | solves geometry with coordinates",
                            "math | checks algebraic manipulations step by step",
                            "code | writes Python data pipelines",
                            "code | debugs concurrency bugs in Rust",
                            "code | reviews TypeScript front-end code",
                            "chat | translates between Cantonese and English",
                            "chat | summarises long documents politely",
                            "chat | answers customer-support questions",
                        ]
                    ),
                )
                with gr.Row():
                    ck = gr.Slider(1, 12, value=2, step=1, label="k")
                    ctau = gr.Slider(0.05, 1.0, value=0.25, step=0.05, label="τ")
                    cbtn = gr.Button("Embed + build graph")
                with gr.Row():
                    cstats = gr.Dataframe(scale=2)
                    cfig = gr.Plot(scale=3)
                cdyn = gr.LinePlot(x="step", y="value", color="series", height=260)
            gouts = [gstats, gfig, gdyn]
            gbtn.click(frame1_graph, [gk, gtau], gouts)
            cbtn.click(custom_graph, [ctext, ck, ctau], [cstats, cfig, cdyn])

        with gr.Tab("🩺 Diagnostics"):
            diag = gr.Markdown(diagnostics_md())
            sc_btn = gr.Button("Run LFM self-check (HF vs llama.cpp vectors)")
            sc_out = gr.Markdown()
            sc_btn.click(lfm_selfcheck, [], sc_out)

    demo.load(frame1_graph, None, gouts)

if __name__ == "__main__":
    demo.queue(max_size=24)
    try:
        demo.launch(theme=gr.themes.Citrus(), css=CSS, mcp_server=True, ssr_mode=False)
    except (ImportError, ModuleNotFoundError, TypeError):
        demo.launch(theme=gr.themes.Citrus(), css=CSS, ssr_mode=False)
