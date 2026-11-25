"""Cross-tokenizer overlap analysis for frame1-Hexad (real tokenizer.json files)."""
import json, itertools
from tokenizers import Tokenizer
TOKDIR = __import__("os").environ.get("HEXAD_TOKDIR", "hexad_work/tok")
M = {
 "K2":     "IFM_K2-Horizon-0.9B",
 "R1":     "deepseek-ai_DeepSeek-R1-Distill-Qwen-1.5B",
 "Limite": "paradigma-inc_limite-1b-violetto",
 "Qwen35": "Qwen_Qwen3.5-0.8B",
 "LFM":    "LiquidAI_LFM2.5-Embedding-350M",
}
V = {}
for k, d in M.items():
    with open(f"{TOKDIR}/{d}/tokenizer.json", encoding="utf-8") as fh:
        V[k] = json.load(fh)["model"]["vocab"]
T = {k: Tokenizer.from_file(f"{TOKDIR}/{d}/tokenizer.json") for k, d in M.items()}
print("== base vocab sizes ==")
for k, v in V.items(): print(f"{k:7s} {len(v)}")
print("\n== token-string overlap (byte-level strings) ==")
for a, b in itertools.combinations(M, 2):
    inter = set(V[a]) & set(V[b])
    print(f"{a:7s} ∩ {b:7s} = {len(inter):7d}   ({len(inter)/len(V[a]):.1%} of {a}, {len(inter)/len(V[b]):.1%} of {b})")
same_id = sum(1 for s, i in V["R1"].items() if V["Limite"].get(s) == i)
print(f"\nR1 vs Limite identical (string,id) pairs: {same_id} / {len(V['R1'])}")
samples = {
 "en":   "The derivative of x^3 + 2x is 3x^2 + 2, so the critical points satisfy 3x^2 + 2 = 0, which has no real roots.",
 "math": "Let ABC be a triangle with AB=11, AC=13. Compute BC given that the orthocenter of AGO lies on BC. \\frac{a^2+b^2}{2}=\\sqrt{3}",
 "code": "def knn_laplacian(emb, k):\n    n = len(emb)\n    for i in range(n):\n        sims = [cos(emb[i], emb[j]) for j in range(n)]\n",
 "zh":   "設 agent 數目為 128，步數為 640，用 frame1 作為骨骼，將六個模型合成一個系統。",
 "yue":  "我哋將 ring laplacian 換成由 LFM 向量建嘅 k-NN 圖拉普拉斯，睇下收斂快咗幾多。",
}
def bounds(tok, text):
    enc = tok.encode(text, add_special_tokens=False)
    return [o[1] for o in enc.offsets], len(enc.ids), enc.tokens
print("\n== per-sequence alignment K2 vs R1 (same text) ==")
print("lang  | K2 tok | R1 tok | shared boundaries | identical tokens (same span)")
for lang, s in samples.items():
    bk, nk, tk = bounds(T["K2"], s); br, nr, tr = bounds(T["R1"], s)
    ek = T["K2"].encode(s, add_special_tokens=False); er = T["R1"].encode(s, add_special_tokens=False)
    spk = set(map(tuple, ek.offsets)); spr = set(map(tuple, er.offsets))
    shared_b = len(set(bk) & set(br)) / max(len(set(bk) | set(br)), 1)
    same = len(spk & spr) / max(nr, 1)
    print(f"{lang:5s} | {nk:6d} | {nr:6d} | {shared_b:17.1%} | {same:.1%} of R1 tokens")
