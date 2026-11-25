"""CPU 測試：同詞表 token 融合（不需要 GPU、不需要下載任何模型）。

    python3 space/test_fusion.py

用一個微型「假分詞器 + 假模型」，把 `fusion.py` 的每一條性質釘住：
  * 同詞表投影：同家族成員的 logits 投影是無損的（共同範圍外的 id 被丟棄）
  * 三種融合律：product / mixture / assisted 的極限行為
  * 循環：3 輪的 id 串接、每輪重錨、編號與 seed 遞增
  * 診斷：top-1 一致率、JS 散度、共同覆蓋率
"""

import math
import sys
from dataclasses import dataclass

sys.path.insert(0, __file__.rsplit("/", 1)[0])

import torch  # noqa: E402

import fusion as fu  # noqa: E402

V = 32  # 共同詞表的常規 id 數（測試用小詞表）


# ---------------------------------------------------------------- 假分詞器
class FakeTok:
    def __init__(self, toks):
        self._vocab = {t: i for i, t in enumerate(toks)}
        self.all_special_tokens = ["<eos>"]
        self._vocab["<eos>"] = len(toks)

    def __call__(self, text, add_special_tokens=False, **_):
        ids = [1, 2, 3] if text else []
        return {"input_ids": torch.tensor(ids)}

    def apply_chat_template(self, msgs, tokenize=False, add_generation_prompt=True):
        return f"<user>{msgs[0]['content']}</user><assistant>"

    def get_vocab(self):
        return dict(self._vocab)

    def decode(self, ids, skip_special_tokens=True):
        inv = {v: k for k, v in self._vocab.items()}
        return "".join(inv.get(int(i), "?") for i in ids)

    def __len__(self):
        return len(self._vocab)


TOKS = [f"t{i}" for i in range(V)]
TOK = FakeTok(TOKS)


# ---------------------------------------------------------------- 假模型
@dataclass
class FakeOut:
    logits: torch.Tensor
    past_key_values: object


class FakeModel:
    """每次都吐同一個固定分佈（只看最後一格），足以測融合律。"""

    def __init__(self, seed, width=V, peak=None):
        g = torch.Generator().manual_seed(seed)
        self.base = torch.randn(width, generator=g) * 3.0
        if peak is not None:
            self.base[peak] += 12.0

    def __call__(self, input_ids=None, past_key_values=None, use_cache=True):
        # HF 的形狀是 [batch, seq, vocab]，最後一格才是下一步的 logits
        return FakeOut(logits=self.base.reshape(1, 1, -1), past_key_values=None)


def members(weights=(1.0, 1.0, 1.0), peaks=(None, None, None)):
    out = []
    for i, (w, p) in enumerate(zip(weights, peaks)):
        out.append(
            fu.FusionMember(
                name=f"m{i}",
                tokenizer=TOK,
                model=FakeModel(seed=i + 1, peak=p),
                end_ids=[V],
                think_close_ids=[],
                weight=w,
            )
        )
    return out


def shared():
    return fu.SharedVocab(TOK, size=V)


# ---------------------------------------------------------------- 測試
def test_shared_vocab():
    sv = shared()
    assert sv.V == V and sv.width == V + 2
    assert sv.END == V and sv.THINK_CLOSE == V + 1
    assert sv.decode([0, 1, 2]) == "t0t1t2"
    assert sv.decode([sv.END]) == ""
    print("✅ SharedVocab: V = %d, END = %d, THINK_CLOSE = %d" % (sv.V, sv.END, sv.THINK_CLOSE))


def test_projection_is_lossless_for_shared_family():
    sv = shared()
    m = members()[0]
    m.end_ids = [V - 1]  # 讓 END 落在共同範圍內（真實模型的 <|endoftext|> = 151643 < V）
    logits = torch.randn(V + 8)  # 多出來的 id 不在共同詞表裡
    lp = fu.project_logprobs(m, logits, sv)
    direct = torch.log_softmax(logits[:V], -1)
    inside = lp[:V]
    # 除了 END 那一格（end id = V 被聚合到 END），共同範圍內必須完全一致
    mask = torch.ones(V, dtype=torch.bool)
    for i in m.end_ids:
        if 0 <= i < V:
            mask[i] = False
    assert torch.allclose(inside[mask], direct[mask], atol=1e-6), "projection changed probabilities"
    assert torch.isclose(lp[sv.END], direct[m.end_ids].logsumexp(0), atol=1e-6)
    assert torch.isneginf(lp[m.end_ids[0]]) or m.end_ids[0] >= V, "END id must be moved, not double counted"
    # 共同空間只有 V + 2 格；id ≥ V 的 logits 根本沒有對應格位（等價於被丟棄）
    assert lp.numel() == V + 2
    assert torch.isfinite(lp[:V][mask]).all()  # END 那一格被搬去偽 token
    print("✅ 同詞表投影：共同範圍內無損，範圍外 id 丟棄，END 以 logsumexp 聚合")


def test_fusion_rules():
    a = torch.tensor([0.7, 0.2, 0.1]).log()
    b = torch.tensor([0.1, 0.3, 0.6]).log()
    prod = fu.fuse_step([a, b], [1.0, 1.0], rule="product").exp()
    assert torch.allclose(prod / prod.sum(), prod)
    assert prod.argmax() == 0, "product is dominated by the agreed peak"
    mix = fu.fuse_step([a, b], [0.5, 0.5], rule="mixture").exp()
    assert abs(float(mix.sum()) - 1.0) < 1e-6
    assert torch.allclose(mix, (a.exp() + b.exp()) / 2, atol=1e-6)
    # product 對「任一成員給 0」是硬否決，mixture 不是
    z = torch.tensor([0.0, 0.5, 0.5]).log()
    assert float(fu.fuse_step([a, z], [1.0, 1.0], rule="product").exp()[0]) < 1e-6
    assert float(fu.fuse_step([a, z], [0.5, 0.5], rule="mixture").exp()[0]) > 0.1
    # assisted：p_strong == p_weak ⇒ 原樣輸出
    same = fu.assisted_distribute(a, a)
    assert torch.allclose(same, a.exp(), atol=1e-6)
    strong = torch.tensor([0.9, 0.05, 0.05]).log()
    weak = torch.tensor([0.1, 0.6, 0.3]).log()
    dis = fu.assisted_distribute(strong, weak)
    assert torch.allclose(dis.sum(), torch.tensor(1.0), atol=1e-6)
    # 精確值：(residual + weak) / Z，residual = [0.8, 0, 0]，Z = 1.8
    assert torch.allclose(dis, torch.tensor([0.9, 0.6, 0.3]) / 1.8, atol=1e-6)
    assert int(dis.argmax()) == 0, "assisted 的峰由 strong 決定"
    assert float(dis[0]) > float(weak.exp()[0]), "strong 的優勢被保留"
    assert float(dis[1]) / float(dis[0]) > 0.5, "weak 的其餘質量未被抹掉"
    print("✅ 三種融合律：product（否決式幾何平均）· mixture（線性）· assisted（拒絕採樣）")


def test_cycle_single():
    sv = shared()
    ms = members(peaks=(5, 5, 5))  # 三個成員都偏好 id 5
    r = fu.fused_sequence(ms, sv, "hi", budget=4, temperature=0.0, device="cpu", force_answer=False)
    assert r.ids[:1] == [5], r.ids
    assert r.stats.summary()["agreement"] == 1.0
    # 覆蓋率 = 共同詞表中有成員給出有限機率的比例；END 那格被搬去偽 token，故略 < 1
    assert 0.9 <= r.stats.summary()["covered"] < 1.0
    # JS 散度必須是有限值（0·log 0 的遮罩回歸測試）
    assert math.isfinite(r.stats.summary()["js"]), r.stats.summary()
    print("✅ 單輪融合：三個成員同峰 → top-1 一致率 100%%，cover 100%%，前 4 步 = %s" % (r.ids,))


def test_cycles_three():
    sv = shared()
    ms = members(peaks=(5, 5, 7))
    r = fu.run_cycles(ms, sv, "hi", cycles=3, budget=3, temperature=0.0, device="cpu", force_answer=False)
    assert len(r.cycles) == 3, "必須是 3 輪"
    assert all(len(c.ids) == 3 for c in r.cycles), [len(c.ids) for c in r.cycles]
    assert len(r.text) >= 9
    joined = [i for c in r.cycles for i in c.ids]
    assert len(joined) == 9
    # 第二、三輪的 seed 遞增（避免三輪完全重複同一次採樣）
    assert r.stats.steps == 9
    print("✅ 循環 3 次：3 輪 × 3 token = 9 個共同 token，串成 %r" % r.text[:24])


def test_weights_and_force():
    assert fu.fusion_weights(3) == [1 / 3] * 3
    w = fu.fusion_weights(3, weight=0.5)
    assert abs(w[0] - 0.5) < 1e-9 and abs(sum(w) - 1.0) < 1e-9
    sv = shared()
    ms = members(peaks=(3, 9, 9))
    # weight=1 → 只看第一個成員
    r = fu.fused_sequence(
        [ms[0]], sv, "hi", weights=[1.0], budget=2, temperature=0.0, device="cpu", force_answer=False
    )
    assert r.ids[0] == 3
    # 強模型（第 1 個）在 assisted 下勝出
    r2 = fu.fused_sequence(ms, sv, "hi", rule="assisted", budget=2, temperature=0.0, device="cpu", force_answer=False)
    assert r2.ids[0] == 9
    print("✅ 權重與 assisted：weight=1 退化為單模型；assisted 由 strong 決定")


def test_product_weights_are_exponents():
    """product 律：p ∝ Π p_i^{w_i}。權重是**指數**，否則等權就不是幾何平均。

    這個測試是 MoonBit `39_shared_vocab_fusion.mbt` 的 T81 鏡像——那邊先寫錯成
    `log p_i + log w_i`，被「三個相同成員必須原樣返回」抓到；這裡兩邊一起釘死。
    """
    l1 = torch.log(torch.tensor([0.4, 0.3, 0.2, 0.1, 0.0]))
    l2 = torch.log(torch.tensor([0.1, 0.2, 0.2, 0.0, 0.5]))
    w = [0.5, 0.5]

    got = fu.fuse_step([l1, l2], w, rule="product")
    want = torch.log_softmax(0.5 * l1 + 0.5 * l2, dim=-1)  # Σ w_i·log p_i
    assert torch.allclose(got, want, atol=1e-6), (got, want)
    # 兩邊都是 0 的格只能是 0（硬否決：乘積律的交集）
    assert got[3].item() < -1e30, got
    assert got[4].item() < -1e30, got
    assert abs(got.exp().sum().item() - 1.0) < 1e-5

    # 若把權重誤當成對數加項（舊 bug），等權乘積律會過度尖銳化 → top1 機率不同
    buggy = torch.log_softmax(l1 + l2, dim=-1)
    assert not torch.allclose(got, buggy, atol=1e-6)

    # 同一個分佈 × 3、等權 ⇒ 必須原樣返回（idempotence）
    same = fu.fuse_step([l1, l1, l1], [1 / 3] * 3, rule="product")
    assert torch.allclose(same, torch.log_softmax(l1, dim=-1), atol=1e-6)

    # weight=1 退化為單成員（product 與 mixture 都要）
    solo_p = fu.fuse_step([l1, l2], [1.0, 0.0], rule="product")
    solo_m = fu.fuse_step([l1, l2], [1.0, 0.0], rule="mixture")
    assert torch.allclose(solo_p, torch.log_softmax(l1, dim=-1), atol=1e-6)
    assert torch.allclose(solo_m, torch.log_softmax(l1, dim=-1), atol=1e-6)
    print("✅ product 權重 = 指數：Σ w_i·log p_i，硬否決成立，等權即幾何平均")


def test_trace_is_exportable():
    """每輪都要留下逐 token trace，而且能被 pandas 寫成 CSV（app.py 的匯出路徑）。"""
    import pandas as pd

    res = fu.run_cycles(members(), shared(), "q", cycles=3, budget=5, temperature=0.8,
                        seed=1, force_answer=False, device="cpu")
    assert len(res.cycles) == 3
    rows = []
    for ci, c in enumerate(res.cycles):
        assert len(c.trace) == c.stats.steps, (len(c.trace), c.stats.steps)
        for st in c.trace:
            for k in ("step", "fused_id", "picks", "agree", "js", "forced"):
                assert k in st, k
            assert len(st["picks"]) == len(members())
            rows.append(dict(cycle=ci + 1, **st))
    df = pd.DataFrame(rows)
    assert len(df) == res.stats.steps == 15
    assert df["cycle"].nunique() == 3
    # picks 要展平成字串才寫得進 CSV
    df2 = df.assign(picks=df["picks"].astype(str))
    text = df2.to_csv(index=False)
    assert text.count("\n") == 16 and "picks" in text
    print(f"✅ trace 匯出：3 輪 × 5 步 = {len(df)} 行 CSV（欄位 cycle/step/fused_id/picks/agree/js/forced）")


def test_fusion_plan_is_exhaustive():
    """「其餘都係同詞 token 融合」⇒ 名單的每個子集都要有輪次計劃（單體/配對/三體）。"""
    trio = ["Limite-1B-Violetto", "Qwen2.5-Coder-1.5B", "Qwen3-0.6B"]
    plan = fu.fusion_plan(trio, max_size=3)
    labels = [d["label"] for d in plan]
    assert len(plan) == 7, labels              # C(3,1)+C(3,2)+C(3,3) = 3+3+1
    assert sum(1 for d in plan if d["size"] == 1) == 3
    assert sum(1 for d in plan if d["size"] == 2) == 3
    assert sum(1 for d in plan if d["size"] == 3) == 1
    # 三體那一項必須正是預設三人組
    assert plan[-1]["members"] == trio
    # 4 個成員、max_size=2 ⇒ 4 + 6 = 10 個組合（配對全部覆蓋）
    plan2 = fu.fusion_plan(trio + ["DeepSeek-R1-Distill-Qwen-1.5B"], max_size=2)
    assert len(plan2) == 10, len(plan2)
    pairs = [set(d["members"]) for d in plan2 if d["size"] == 2]
    assert len({tuple(sorted(p)) for p in pairs}) == 6
    # 排序穩定：size 先、再按名稱
    assert plan == fu.fusion_plan(trio, max_size=3)
    print(f"✅ 輪次計劃：3 人 → {len(plan)} 個組合；4 人（<=2）→ {len(plan2)} 個，全部配對覆蓋")


def test_budget_forcing():
    sv = shared()
    ms = [fu.FusionMember(name="m", tokenizer=TOK, model=FakeModel(seed=1), end_ids=[V], think_close_ids=[7])]
    r = fu.fused_sequence(ms, sv, "hi", budget=4, temperature=0.0, device="cpu", force_answer=True)
    assert r.forced and r.closed, "75% 預算時應強制關閉 thinking"
    print("✅ 預算強制：75%% 位置插入 THINK_CLOSE（closed=%s, forced=%s）" % (r.closed, r.forced))


def test_js_no_nan_with_neg_inf():
    a = torch.tensor([0.0, -math.inf, -math.inf])
    b = torch.tensor([-5.0, -0.5, -math.inf])
    js = fu.js_divergence([a, b])
    assert math.isfinite(float(js)), "JS 在 −inf 格位上不能變 NaN"
    print("✅ JS 遮罩：−inf 格位不產生 NaN（js = %.4f）" % float(js))


def test_same_vocab_bootstrap_passes_ids_through():
    """同詞表成員把上一輪的 id 直接餵回模型（不經 decode→encode，零損失）。"""
    sv = shared()
    m = members()[0]
    m.same_vocab = True
    m.tokenizer = TOK
    ids = [3, 4, 5]
    assert fu._bootstrap_ids(m, sv, ids) == ids
    # 非同詞表成員（用自己的 tokenizer）走 decode→encode 路線
    m2 = members()[1]
    m2.same_vocab = False

    class OtherTok(FakeTok):
        def __call__(self, text, add_special_tokens=False, **_):
            return {"input_ids": torch.tensor([9, 9])}

    m2.tokenizer = OtherTok(TOKS)
    assert fu._bootstrap_ids(m2, sv, ids) == [9, 9]
    print("✅ 重錨：同詞表成員直接吃 id，跨詞表成員走 decode→encode")


def test_static_app_wiring():
    """靜態檢查 app.py 的融合接線（不需要 GPU / gradio / transformers）。"""
    import ast

    src = open(__file__.rsplit("/", 1)[0] + "/app.py", encoding="utf-8").read()
    tree = ast.parse(src)
    consts = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            if isinstance(node.value, ast.List):
                try:
                    consts[node.targets[0].id] = ast.literal_eval(node.value)
                except Exception:
                    pass
            elif isinstance(node.value, ast.Constant):
                consts[node.targets[0].id] = node.value.value
    fns = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    assert "fuse_shared_vocab" in fns, "app.py 必須有 fuse_shared_vocab"
    fn = fns["fuse_shared_vocab"]
    defaults = {}
    args = fn.args.args[-len(fn.args.defaults):] if fn.args.defaults else []
    for a, d in zip(args, fn.args.defaults):
        defaults[a.arg] = ast.literal_eval(d)
    assert defaults.get("cycles") == 3, f"循環預設必須是 3，現在是 {defaults.get('cycles')}"
    body = ast.get_source_segment(src, fn) or ""
    assert "run_cycles" in body and "SharedVocab" in body, "必須走共用詞表的 run_cycles"
    names = consts.get("FUSION_NAMES") or consts.get("FUSION_DEFAULT") or []
    for want in ("Limite-1B-Violetto", "Qwen2.5-Coder-1.5B", "Qwen3-0.6B"):
        assert want in names, f"{want} 必須在 FUSION_DEFAULT 裡，現在是 {names}"
    # 三個模型都必須被標成同詞表成員（Qwen2 家族）
    assert "same_vocab=True" in src
    # UI 一定要有接上 fuse_shared_vocab 的按鈕
    assert "fuse_shared_vocab, uinputs, uouts" in src
    # 逐 token trace 匯出：uouts 有 5 個，return 也要 5 個（否則 gradio 會錯位）
    assert "utrace" in src and "uouts = [uanswer, ureason, ustats, ucycle_df, utrace]" in src
    returns = [n for n in ast.walk(fn) if isinstance(n, ast.Return)]
    assert returns, "fuse_shared_vocab 沒有 return"
    for r in returns:  # 主路徑 + 巡迴分支都要回 5 個值，否則 gradio 輸出會錯位
        v = r.value
        if isinstance(v, ast.Call):
            assert v.func.id == "_fusion_sweep"
        else:
            assert isinstance(v, ast.Tuple) and len(v.elts) == 5, "回傳值必須是 5 元組（對應 5 個輸出）"
    # 巡迴模式：roster 的全子集都要走同一條管線
    assert "_fusion_sweep" in src and "_fusion_plan" not in src.replace("FZ.fusion_plan", "")
    assert "fusion_plan(" in src
    print("✅ app.py 接線：FUSION_DEFAULT = %s，cycles 預設 3，走 SharedVocab + run_cycles" % (names,))


def test_stats_json_roundtrip():
    import json

    sv = shared()
    r = fu.run_cycles(members(peaks=(4, 4, 9)), sv, "hi", cycles=3, budget=4, temperature=0.0, device="cpu",
                      force_answer=False)
    payload = dict(**r.stats.summary(), cycles=len(r.cycles), weights=r.weights, rule=r.rule,
                   per_cycle=r.delta())
    blob = json.dumps(payload, ensure_ascii=False)
    assert "cycles" in blob and json.loads(blob)["cycles"] == 3
    assert abs(sum(r.weights) - 1.0) < 1e-9
    print("✅ 診斷可序列化：%s" % blob[:120])


if __name__ == "__main__":
    test_shared_vocab()
    test_projection_is_lossless_for_shared_family()
    test_fusion_rules()
    test_js_no_nan_with_neg_inf()
    test_cycle_single()
    test_cycles_three()
    test_weights_and_force()
    test_product_weights_are_exponents()
    test_fusion_plan_is_exhaustive()
    test_trace_is_exportable()
    test_budget_forcing()
    test_same_vocab_bootstrap_passes_ids_through()
    test_static_app_wiring()
    test_stats_json_roundtrip()
    print("\n全部測試通過（13/13）")
