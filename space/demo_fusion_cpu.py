#!/usr/bin/env python3
"""CPU 冒煙測試 / 演示：同詞表 token 融合 × 3 循環（不需 GPU、不需下載模型）。

    python3 space/demo_fusion_cpu.py

用三個「假成員」（各有偏好的 token）走完 `fusion.run_cycles(cycles=3)`，
把每一輪的 tokens / 字符 / top-1 一致率 / JS 散度印成表格，並寫出
`space/fusion_demo.csv`。**數字是合成的**——它只證明管線接得通、診斷算得對，
不是任何真實模型的表現。
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd  # noqa: E402
import torch  # noqa: E402

import fusion as FZ  # noqa: E402

V = 64
TOKS = [f"t{i}" for i in range(V)]


class Tok:
    all_special_tokens = ["<eos>"]

    def __init__(self):
        self._v = {t: i for i, t in enumerate(TOKS)}
        self._v["<eos>"] = V

    def __call__(self, text, add_special_tokens=False, **_):
        return {"input_ids": torch.tensor([1, 2, 3])}

    def apply_chat_template(self, msgs, tokenize=False, add_generation_prompt=True):
        return f"<user>{msgs[0]['content']}</user><assistant>"

    def get_vocab(self):
        return dict(self._v)

    def decode(self, ids, skip_special_tokens=True):
        inv = {v: k for k, v in self._v.items()}
        return "".join(inv.get(int(i), "?") for i in ids)

    def __len__(self):
        return len(self._v)


class Model:
    """固定分佈；`peaks` 讓它偏好幾個 token（模擬「三個模型各有專長」）。"""

    def __init__(self, seed, peaks):
        g = torch.Generator().manual_seed(seed)
        self.base = torch.randn(V, generator=g) * 2.0
        for p in peaks:
            self.base[p] += 9.0

    def __call__(self, input_ids=None, past_key_values=None, use_cache=True):
        class O:
            logits = self.base.reshape(1, 1, -1)
            past_key_values = None

        return O()


def main():
    tok = Tok()
    shared = FZ.SharedVocab(tok, size=V)
    trio = [
        ("Limite-1B-Violetto", [7, 8, 9]),      # 數學式推理
        ("Qwen2.5-Coder-1.5B", [20, 21]),       # 程式碼 token
        ("Qwen3-0.6B", [7, 40]),                # 與 Limite 有重疊、也各有主張
    ]
    members = [
        FZ.FusionMember(name=n, tokenizer=tok, model=Model(seed=i + 1, peaks=p), end_ids=[V],
                        same_vocab=True)
        for i, (n, p) in enumerate(trio)
    ]

    rows = []
    for rule in ("product", "mixture", "assisted"):
        res = FZ.run_cycles(members, shared, "demo", cycles=3, rule=rule, budget=8,
                            temperature=0.8, seed=0, force_answer=False, device="cpu")
        for i, c in enumerate(res.cycles):
            rows.append(
                dict(
                    rule=rule,
                    cycle=i + 1,
                    tokens=c.stats.steps,
                    text=c.text,
                    agreement=round(c.stats.agreement / max(c.stats.steps, 1), 3),
                    js=round(c.stats.js / max(c.stats.steps, 1), 4),
                    covered=round(c.stats.covered / max(c.stats.steps, 1), 3),
                )
            )
        print(f"\n[{rule}] 3 輪合併 = {res.stats.steps} tokens，text = {res.text!r}")
        print("  weights =", [round(w, 3) for w in res.weights])

    df = pd.DataFrame(rows)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fusion_demo.csv")
    df.to_csv(out, index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nwrote {out}")
    print("注意：以上是假模型（合成分佈）的冒煙測試，不是真實模型輸出。")


if __name__ == "__main__":
    main()
