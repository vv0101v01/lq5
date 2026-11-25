"""同詞表 token 級融合 (shared-vocabulary token fusion) — Limite · Qwen2.5-Coder-1.5B · Qwen3 …

本模組把 `space/app.py` 原本只有 Limite⊕R1 一條硬編碼管線的「同詞表 logit 融合」抽出來，
變成一條**任意條數**、**可循環 N 次**的通用管道：

  1. **同詞表** (shared vocabulary)。所有成員模型都在 *同一個* 標籤空間裡投票：
     以一個「詞表擁有者」的 tokenizer 為標籤座標系（本專案為 Limite／Qwen2 家族的
     151 643 個常規 BPE id），每個成員模型的 logits 都被投影到這個空間。

  2. **token 級融合**。每一步只生成 **一個** token：先取各成員在共同詞表上的
     log-softmax，再合成一個分佈，最後從這個分佈採樣。支援三種合成律：
       * `product`  (log-linear)  `p ∝ Π_i p_i^{w_i}`          —— 幾何平均／專家乘積
       * `mixture`  (linear)      `p ∝ Σ_i w_i · p_i`           —— 線性混合
       * `assisted` (speculative) 強模型提出、弱模型驗證：
             `p(x) ∝ max(0, p_strong(x) − p_weak(x))_+ + p_weak(x)`
          （拒絕採樣的分佈式寫法，強模型只需在其 top-k 上給出機率，弱模型給全詞表）

  3. **循環**。一輪 = 一次融合解碼。第 n 輪把第 n−1 輪的融合結果（共同 id 序列）
     重餵給每個成員（`bootstrap` 決定餵法），再往後解碼，直到 `max_new_tokens`
     或答案收斂。`run_cycles(..., cycles=3)` 就是「循環重複進行 3 次」。

  4. **可移植性**。Qwen2.5-Coder-1.5B／Qwen3 與 Limite 同屬 Qwen2 BPE 家族，
     共用 id；`tokenizer="shared"` 表示「直接吃共同 id」，`tokenizer="own"` 表示
     「用自己的 tokenizer 重新切詞」（id 不同時用，會有輕微損失，但融合仍然成立）。

設計原則：**沒有任何 GPU 邏輯藏在裡面**。`fuse_step` / `assisted_distribute` /
`fused_sequence` 都是純張量函數，可以在 CPU 上用微型模型與假 tokenizer 完整測試
（見 `space/test_fusion.py`）；`app.py` 只負責載入模型與畫 UI。
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence

import torch
import torch.nn.functional as F

# --------------------------------------------------------------------------------------
# 0. 常數與小工具
# --------------------------------------------------------------------------------------


THINK_CLOSE_TOKENS = ("</think>", "<｜end▁of▁thinking｜>", "<｜end▁of▁thinking｜>")
"""各家族 `</think>` 的寫法：Qwen3 是 `<｜end▁of▁thinking｜>` 這個特殊 token，Limite 是普通文字。"""


@dataclass
class FusionMember:
    """一個融合成員：模型 + tokenizer + 它在共同詞表上的可用範圍。

    Attributes
    ----------
    name        : 顯示名
    tokenizer   : 「own」= 與詞表擁有者共用 tokenizer（直接吃共同 id）；
                  「own」以外的情況請傳一個帶 `apply_chat_template` 的 tokenizer 物件，
                  並把 `tokenizer` 設成該物件本身。
    model       : 任何 `model(input_ids=…, past_key_values=…, use_cache=True)` 可呼叫物件
    prompt_suffix : 官方模板之後要硬接的字串（例如 Limite 的 `<think>\n`）
    end_ids     : 該模型的句子結束 id（會被映射到偽 token END）
    think_close_ids : 該模型的 `</think>` id 清單（會被映射到偽 token THINK_CLOSE）
    weight      : 融合權重（僅在 `product` / `mixture` 使用）
    same_vocab  : 該成員的 id 與共同詞表逐個對齊（Qwen2 家族為 True）：
                  此時「上一輪的共同 id」可以直接餵回去，不必解碼再切詞
    """

    name: str
    tokenizer: object
    model: object
    prompt_suffix: str = ""
    end_ids: Sequence[int] = field(default_factory=list)
    think_close_ids: Sequence[int] = field(default_factory=list)
    weight: float = 1.0
    same_vocab: bool = False
    # 執行期填充
    _cache: object = None
    _logits: Optional[torch.Tensor] = None
    _len: int = 0


# --------------------------------------------------------------------------------------
# 1. 共同詞表
# --------------------------------------------------------------------------------------


class SharedVocab:
    """共同詞表：以 `owner_tokenizer` 的常規 id 為標籤座標系。

    偽 token 追加在尾端：`END = V`、`THINK_CLOSE = V + 1`，故融合分佈長度為 `V + 2`。
    """

    def __init__(self, owner_tokenizer, size: Optional[int] = None):
        vocab = owner_tokenizer.get_vocab()
        specials = set(getattr(owner_tokenizer, "all_special_tokens", []) or [])
        self.owner = owner_tokenizer
        self.vocab: Dict[str, int] = dict(vocab)
        # 可用的「常規」id 上限：詞表大小或被明確指定的共同大小
        self.V = int(size if size is not None else _vocab_size(owner_tokenizer))
        self.id_to_token: Dict[int, str] = {}
        for tok, idx in self.vocab.items():
            if tok in specials:
                continue
            if 0 <= idx < self.V and idx not in self.id_to_token:
                self.id_to_token[idx] = tok
        self.END = self.V
        self.THINK_CLOSE = self.V + 1
        self.width = self.V + 2

    # -- 別名（同義 token 合併，讓不同模型的切詞差異不至於讓某個 id 恆為 −inf） -------------
    def add_aliases(self, aliases: Dict[str, Iterable[str]]) -> int:
        """`{canonical: [alias, …]}`：把別名的 id 指向 canonical 的 id。回傳別名數。"""
        added = 0
        for canonical, names in aliases.items():
            cid = self.vocab.get(canonical)
            if cid is None:
                continue
            for a in names:
                aid = self.vocab.get(a)
                if aid is not None and aid != cid:
                    self.id_to_token[aid] = canonical
                    added += 1
        return added

    def encode(self, text: str, **kw) -> List[int]:
        return [int(i) for i in self.owner(text, add_special_tokens=False, **kw)["input_ids"]]

    def decode(self, ids: Sequence[int]) -> str:
        clean = [int(i) for i in ids if 0 <= int(i) < self.V]
        if not clean:
            return ""
        return self.owner.decode(clean, skip_special_tokens=True)


def _vocab_size(tok) -> int:
    n = len(tok)
    try:
        return int(n)
    except TypeError:  # 極少數 tokenizer 不支援 len()
        v = tok.get_vocab()
        return int(max(v.values())) + 1 if v else 0


# --------------------------------------------------------------------------------------
# 2. 投影：任一成員的 logits → 共同空間的 log-prob（長度 width = V + 2）
# --------------------------------------------------------------------------------------


def project_logprobs(
    member: FusionMember,
    logits: torch.Tensor,
    shared: SharedVocab,
    weight: float = 1.0,
    usable: Optional[int] = None,
) -> torch.Tensor:
    """把一個成員的 `logits[V_model]` 投影成共同空間上的加權 log 機率 `[width]`。

    * 超出共同範圍的 id 直接丟棄（它們不是共同詞表的詞）；
    * `end_ids` 用 logsumexp 聚合成一個 END；
    * `think_close_ids` 同理聚合成 THINK_CLOSE；
    * 若該成員的 logits 整個落在共同範圍內（同詞表家族），這就是無損投影。
    """
    logits = logits.float().reshape(-1)
    if usable is None:
        usable = int(logits.numel())
    n = min(int(logits.numel()), int(usable), shared.V)
    lp = torch.log_softmax(logits[:n], dim=-1)

    out = torch.full((shared.width,), float("-inf"), dtype=torch.float32, device=logits.device)
    out[:n] = lp

    def _gather(ids: Iterable[int]) -> Optional[torch.Tensor]:
        keep = [int(i) for i in ids if 0 <= int(i) < n]
        if not keep:
            return None
        return torch.logsumexp(lp[torch.tensor(keep, device=lp.device, dtype=torch.long)], dim=-1)

    e = _gather(member.end_ids)
    if e is not None:
        out[shared.END] = e
        # 從常規 id 中移除（避免同一個結尾被算兩次）
        for i in member.end_ids:
            if 0 <= int(i) < n:
                out[int(i)] = float("-inf")
    t = _gather(member.think_close_ids)
    if t is not None:
        out[shared.THINK_CLOSE] = t
        for i in member.think_close_ids:
            if 0 <= int(i) < n:
                out[int(i)] = float("-inf")

    if weight == 1.0:
        return out
    return out + _log(weight)


def _log(x: float) -> float:
    return torch.log(torch.tensor(max(float(x), 1e-12))).item()


# --------------------------------------------------------------------------------------
# 3. 三種融合律
# --------------------------------------------------------------------------------------


def fusion_weights(n: int, weight: Optional[float] = None) -> List[float]:
    """第 i 個成員的權重。`weight` 給定時只有第 0 個成員為 `weight`，其餘均分 `1−weight`。"""
    if n <= 0:
        return []
    if weight is None:
        return [1.0 / n] * n
    w = min(max(float(weight), 0.0), 1.0)
    if n == 1:
        return [1.0]
    rest = (1.0 - w) / (n - 1)
    return [w] + [rest] * (n - 1)


RULE_ALIASES = {
    "product": "product",
    "product (log-linear)": "product",
    "log-linear": "product",
    "mixture": "mixture",
    "mixture (linear)": "mixture",
    "linear": "mixture",
    "assisted": "assisted",
    "assisted (speculative)": "assisted",
    "speculative": "assisted",
}


def canonical_rule(rule: str) -> str:
    """把 UI 上各式寫法收斂成三種：product / mixture / assisted。"""
    r = (rule or "").strip().lower()
    if r in RULE_ALIASES:
        return RULE_ALIASES[r]
    for k, v in RULE_ALIASES.items():
        if r.startswith(k[:4]):
            return v
    raise ValueError(f"unknown fusion rule: {rule}")


def fuse_step(
    logprob_list: Sequence[torch.Tensor],
    weights: Sequence[float],
    rule: str = "product",
    strong: int = 0,
    weak: int = 1,
) -> torch.Tensor:
    rule = canonical_rule(rule)
    """把各成員在共同空間上的 log 機率合成一個 log 機率。

    rule ∈ {"product", "mixture", "assisted"}；後者用 `strong` / `weak` 兩個索引。
    """
    if not logprob_list:
        raise ValueError("fuse_step needs at least one member")
    if rule == "product":
        # p ∝ Π p_i^{w_i} ⇒ log p = Σ w_i·log p_i  （權重是指數，不是加項！）
        # 這裡必須用乘法：等權 (Σw = 1) 時才是幾何平均；若誤寫成 +log w，
        # 等權的乘積律會悄悄退化成「每個成員指數 1」，把分佈過度尖銳化。
        total = None
        for lp, w in zip(logprob_list, weights):
            if float(w) <= 0.0:
                continue
            term = lp * float(w)
            total = term if total is None else total + term
        if total is None:
            raise ValueError("all fusion weights are zero")
        return F.log_softmax(total, dim=-1)
    if rule == "mixture":
        # 權重 0 的成員必須**完全**移除。用 +log(max(w,1e-12)) 近似會在 w=0 時
        # 留下一個 tiny 的質量（−∞ 格變成 −27.6 而不是 −∞），令 veto 語義失真。
        parts = [lp + _log(w) for lp, w in zip(logprob_list, weights) if float(w) > 0.0]
        if not parts:
            raise ValueError("all fusion weights are zero")
        return F.log_softmax(torch.logsumexp(torch.stack(parts), dim=0), dim=-1)
    if rule == "assisted":
        # 推測解碼：強模型提出、弱模型驗證（拒絕採樣分佈式寫法）
        return F.log_softmax(assisted_distribute(logprob_list[strong], logprob_list[weak]).log(), dim=-1)
    raise ValueError(f"unknown fusion rule: {rule}")


def assisted_distribute(p_strong: torch.Tensor, p_weak: torch.Tensor) -> torch.Tensor:
    """拒絕採樣／推測解碼的分佈式寫法：

        p(x) = max(0, p_strong(x) − p_weak(x))_+ + p_weak(x)

    性質（測試固定）：輸出總和 = 1；`p_strong == p_weak` 時輸出恰等於 `p_strong`；
    因此「弱模型接受全部、強模型只在其優勢處改寫」。
    """
    ps = p_strong.exp() if _is_log(p_strong) else p_strong
    pw = p_weak.exp() if _is_log(p_weak) else p_weak
    residual = (ps - pw).clamp_min(0.0)
    out = residual + pw
    return out / out.sum().clamp_min(1e-12)


def _is_log(x: torch.Tensor) -> bool:
    return bool((x <= 0).all().item())


# --------------------------------------------------------------------------------------
# 4. 診斷：每一步的一致性與 JS 散度
# --------------------------------------------------------------------------------------


@dataclass
class FusionStats:
    steps: int = 0
    agreement: float = 0.0  # top-1 一致率（所有成員兩兩 argmax 相同）
    lead: List[float] = field(default_factory=list)  # 每個成員「自己的 argmax 被選中」次數
    js: float = 0.0
    covered: float = 0.0  # 共同詞表裡至少一個成員給出有限機率的比例
    per_step: List[dict] = field(default_factory=list)

    def summary(self) -> Dict[str, float]:
        n = max(self.steps, 1)
        return dict(
            steps=self.steps,
            agreement=self.agreement / n,
            js=self.js / n,
            covered=self.covered / n,
        )


def js_divergence(logprob_list: Sequence[torch.Tensor]) -> torch.Tensor:
    """各成員分佈對平均分佈的對稱 KL 平均（Jensen–Shannon，單位 nat）。

    小心 `0 · log 0`：`p = 0` 的格位貢獻必須當成 0（用 `where` 遮罩），
    否則 `inf − inf` 會把整個診斷變成 NaN（本輪實測踩到的坑之一）。
    """
    probs = [lp.exp() for lp in logprob_list]  # −inf → 0
    stack = torch.stack(probs)
    mean = stack.mean(0).clamp_min(1e-30)
    logmean = mean.log()  # clamp 之後有限，不會有 −inf − (−inf)
    total = torch.zeros((), device=logmean.device)
    for p in probs:
        positive = p > 0
        safe = torch.where(positive, p * (p.clamp_min(1e-30).log() - logmean), torch.zeros_like(p))
        total = total + safe.sum()
    return total / max(len(logprob_list), 1)


def record_step(stats: FusionStats, logprob_list: Sequence[torch.Tensor], chosen: int) -> None:
    tops = [int(lp.argmax()) for lp in logprob_list]
    ok = all(t == tops[0] for t in tops)
    stats.steps += 1
    stats.agreement += 1.0 if ok else 0.0
    if not stats.lead:
        stats.lead = [0.0] * len(logprob_list)
    for i, lp in enumerate(logprob_list):
        if chosen == tops[i]:
            stats.lead[i] += 1.0
    stats.js += float(js_divergence(logprob_list))
    finite = torch.isfinite(torch.stack(list(logprob_list))).any(0).float().mean()
    stats.covered += float(finite)


# --------------------------------------------------------------------------------------
# 5. 解碼：每步只吐一個 token
# --------------------------------------------------------------------------------------


def _sample(logp: torch.Tensor, temperature: float, top_p: float, gen: torch.Generator) -> int:
    if temperature <= 0:
        return int(torch.argmax(logp))
    k = min(512, int(logp.numel()))
    lv, li = torch.topk(logp, k)
    sp = torch.softmax(lv / temperature, dim=-1)
    sp = sp * ((torch.cumsum(sp, 0) - sp) < top_p)
    if sp.sum() <= 0:
        return int(li[0])
    return int(li[torch.multinomial(sp / sp.sum(), 1, generator=gen)])


@dataclass
class CycleResult:
    ids: List[int]  # 本輪新增的共同 id
    text: str
    stats: FusionStats
    closed: bool  # 是否已出現 THINK_CLOSE
    forced: bool
    trace: List[dict] = field(default_factory=list)  # 逐步診斷（CSV 匯出用）


def prefill(member: FusionMember, input_ids: torch.Tensor) -> None:
    """把 `input_ids` 餵進成員模型並記住 cache / 最後一格 logits。"""
    out = member.model(input_ids=input_ids, use_cache=True)
    member._cache = out.past_key_values
    member._logits = out.logits[0, -1].float()


def step_member(member: FusionMember, ids: Sequence[int]) -> None:
    """吃進一個（或數個）token，更新 cache。"""
    t = torch.tensor([list(ids)], dtype=torch.long, device=member._logits.device)
    out = member.model(input_ids=t, past_key_values=member._cache, use_cache=True)
    member._cache = out.past_key_values
    member._logits = out.logits[0, -1].float()


def _bootstrap_ids(member: FusionMember, shared: SharedVocab, common_ids: Sequence[int]) -> List[int]:
    """把上一輪的共同 id 序列轉成該成員自己吃的 id。"""
    if member.same_vocab or member.tokenizer is shared.owner:
        return [int(i) for i in common_ids if 0 <= int(i) < shared.V]
    text = shared.decode(common_ids)
    if not text:
        return []
    return [int(i) for i in member.tokenizer(text, add_special_tokens=False)["input_ids"]]


def build_anchor(
    member: FusionMember,
    prompt_text: str,
    shared: SharedVocab,
    prev_common: Sequence[int] = (),
    bootstrap: str = "own",
    device: str = "cuda",
) -> torch.Tensor:
    """成員的錨點序列：官方聊天模板 + prompt（+ 上一輪結果）。"""
    tok = member.tokenizer
    prompt_part = tok.apply_chat_template(
        [{"role": "user", "content": prompt_text}], tokenize=False, add_generation_prompt=True
    ) + member.prompt_suffix
    ids = [int(i) for i in tok(prompt_part, add_special_tokens=False)["input_ids"]]
    if prev_common and bootstrap != "none":
        ids += _bootstrap_ids(member, shared, prev_common)
    return torch.tensor([ids], dtype=torch.long, device=device)


def fused_sequence(
    members: Sequence[FusionMember],
    shared: SharedVocab,
    prompt_text: str,
    weights: Optional[Sequence[float]] = None,
    rule: str = "product",
    temperature: float = 0.6,
    top_p: float = 0.95,
    seed: int = 0,
    budget: int = 128,
    prev_common: Sequence[int] = (),
    bootstrap: str = "own",
    strong: int = 0,
    weak: int = 1,
    force_answer: bool = True,
    device: str = "cuda",
    usable: Optional[Sequence[int]] = None,
) -> CycleResult:
    """一輪融合解碼：每一步只生成一個共同 token，直到 END／預算／`</think>`。"""
    if not members:
        raise ValueError("no members to fuse")
    w = list(weights) if weights is not None else fusion_weights(len(members))
    assert len(w) == len(members), "weights must match members"
    rule = canonical_rule(rule)
    gen = torch.Generator(device=device).manual_seed(int(seed))

    for m in members:
        prefill(m, build_anchor(m, prompt_text, shared, prev_common, bootstrap, device))

    stats = FusionStats()
    ids: List[int] = []
    chunks: List[str] = []
    trace: List[dict] = []
    closed, forced = False, False
    answer_at = int(0.75 * budget) if force_answer else budget + 1

    for step in range(int(budget)):
        lps = [
            project_logprobs(m, m._logits, shared, 1.0,
                             usable=usable[i] if usable is not None else None)
            for i, m in enumerate(members)
        ]
        fused = fuse_step(lps, w, rule=rule, strong=strong, weak=weak)
        if force_answer and not closed and step >= answer_at:
            tok_id, forced = shared.THINK_CLOSE, True
        else:
            tok_id = _sample(fused, float(temperature), float(top_p), gen)
        record_step(stats, lps, tok_id)
        trace.append(
            dict(
                step=int(step),
                fused_id=int(tok_id),
                picks=[int(torch.argmax(lp).item()) for lp in lps],
                agree=bool(all(int(torch.argmax(lp).item()) == int(torch.argmax(lps[0]).item()) for lp in lps)),
                js=float(js_divergence(lps).item()) if len(lps) > 1 else 0.0,
                forced=bool(forced and tok_id == shared.THINK_CLOSE),
            )
        )

        if tok_id == shared.END:
            break
        if tok_id == shared.THINK_CLOSE:
            closed = True
            # 每個成員吃自己的 </think>（同詞表成員直接吃共同 id）
            for m in members:
                ids_m = list(m.think_close_ids) or [tok_id]
                step_member(m, ids_m)
            chunks.append("</think>")
            continue
        for m in members:
            step_member(m, [tok_id])
        ids.append(int(tok_id))
        chunks.append(shared.decode([int(tok_id)]))

    text = "".join(chunks)
    return CycleResult(ids=ids, text=text, stats=stats, closed=closed, forced=forced, trace=trace)


@dataclass
class MultiCycleResult:
    cycles: List[CycleResult]
    text: str
    stats: FusionStats
    prompt: str
    rule: str
    weights: List[float]

    def delta(self) -> List[float]:
        """每一輪對最終答案的貢獻（字符數），用來看 3 輪是否遞減收斂。"""
        return [len(c.text) for c in self.cycles]


def run_cycles(
    members: Sequence[FusionMember],
    shared: SharedVocab,
    prompt_text: str,
    cycles: int = 3,
    weight: Optional[float] = None,
    rule: str = "product",
    temperature: float = 0.6,
    top_p: float = 0.95,
    seed: int = 0,
    budget: int = 128,
    bootstrap: str = "own",
    strong: int = 0,
    weak: int = 1,
    force_answer: bool = True,
    device: str = "cuda",
    usable: Optional[Sequence[int]] = None,
) -> MultiCycleResult:
    """**循環重複 N 次**：每一輪都從「prompt + 前幾輪的融合結果」重新出發。

    `cycles=3` 就是「循環重複進行 3 次」。第 n 輪的 `prev_common` 是前面所有輪的串接，
    所以三輪合起來仍然是一條連續的融合解碼，只是每輪重錨一次（bootstrap），
    這正是「多數模型在同一詞表上輪流覆核」的行為。
    """
    w = fusion_weights(len(members), weight)
    acc: List[int] = []
    outs: List[CycleResult] = []
    for c in range(max(int(cycles), 1)):
        r = fused_sequence(
            members,
            shared,
            prompt_text,
            weights=w,
            rule=rule,
            temperature=temperature,
            top_p=top_p,
            seed=int(seed) + c,
            budget=budget,
            prev_common=tuple(acc) if acc else (),
            bootstrap=bootstrap,
            strong=strong,
            weak=weak,
            force_answer=force_answer,
            device=device,
            usable=usable,
        )
        outs.append(r)
        acc.extend(r.ids)
        if r.closed and not r.ids:
            break
    total = FusionStats()
    for r in outs:
        total.steps += r.stats.steps
        total.agreement += r.stats.agreement
        total.js += r.stats.js
        total.covered += r.stats.covered
        total.per_step.extend(r.stats.per_step)
        if not total.lead:
            total.lead = list(r.stats.lead)
        else:
            for i, v in enumerate(r.stats.lead):
                if i < len(total.lead):
                    total.lead[i] += v
    return MultiCycleResult(
        cycles=outs,
        text=shared.decode(acc) if acc else "",
        stats=total,
        prompt=prompt_text,
        rule=rule,
        weights=w,
    )


# --------------------------------------------------------------------------------------
# 6. 便利函數：把 app.py 的 Limite ⊕ R1 變成「三模型 × 3 輪」
# --------------------------------------------------------------------------------------


QWEN_FAMILY_SHARED_V = 151643
"""Limite / Qwen2.5-Coder / Qwen3 共用的常規 BPE id 數（Qwen2 家族 151 643）。"""

def fusion_plan(roster: Sequence[str], max_size: int = 3) -> List[dict]:
    """**全部組合的輪次計劃**：把 `roster` 的每一個大小 1..max_size 的子集列出來。

    用途：用戶要求「其餘都係同詞 token 融合」——不只是 Limite + Qwen2.5-Coder +
    Qwen3 三體，任何單體、任何配對都要走同一條 `SharedVocab` + `run_cycles` 管線
    （而不是文本層 ensemble）。回傳按 (成員數, 名稱) 排序，元素是
    `{"label": ..., "members": [...], "size": n}`，可直接餵進 UI 的下拉選單。
    """
    names = [str(n) for n in roster]
    plan: List[dict] = []
    for size in range(1, max(1, int(max_size)) + 1):
        for combo in itertools.combinations(names, size):
            plan.append(
                dict(
                    size=size,
                    members=list(combo),
                    label=("單體 " if size == 1 else f"{size} 體 ") + " + ".join(combo),
                )
            )
    plan.sort(key=lambda d: (d["size"], d["members"]))
    return plan


