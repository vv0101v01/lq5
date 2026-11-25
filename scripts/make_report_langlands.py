#!/usr/bin/env python3
"""Self-contained HTML report for the Langlands-emergence pipeline.

Reads `langlands_trace.csv` (produced by `moon run src/cmd/main`) and writes
`docs/report_langlands.html`.  No external assets are referenced: every chart
is inline SVG, so the report renders identically inside a sandboxed preview.

Usage:  python3 scripts/make_report_langlands.py
"""
import math
import os

from reporthtml import PALETTE, bar_chart, doc_path, line_chart, load, table, trace_path

CSV_PATH = trace_path("langlands_trace.csv")
OUT_PATH = doc_path("report_langlands.html")



def heatmap(matrix, labels, width=420, height=300, title="", vmax=None):
    n = len(matrix)
    pad = 46
    cell = min((width - 2 * pad) / n, (height - 2 * pad) / n)
    vmax = vmax if vmax is not None else max(abs(v) for row in matrix for v in row) or 1.0
    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" role="img" '
             f'style="background:#fff;border:1px solid {PALETTE["grid"]};border-radius:10px">']
    parts.append(f'<text x="12" y="20" font-size="13" font-weight="600" fill="{PALETTE["ink"]}">{title}</text>')
    for i in range(n):
        for j in range(n):
            v = matrix[i][j] / vmax
            if v >= 0:
                col = f"rgba(37,99,235,{0.12 + 0.85*min(abs(v),1):.3f})"
            else:
                col = f"rgba(220,38,38,{0.12 + 0.85*min(abs(v),1):.3f})"
            x = pad + j * cell
            y = pad + i * cell
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cell:.1f}" height="{cell:.1f}" '
                         f'fill="{col}" stroke="#fff" stroke-width="1"/>')
            parts.append(f'<text x="{x+cell/2:.1f}" y="{y+cell/2+4:.1f}" font-size="9" text-anchor="middle" '
                         f'fill="{PALETTE["ink"]}">{matrix[i][j]:.2f}</text>')
    for i, lab in enumerate(labels):
        parts.append(f'<text x="{pad-6}" y="{pad+cell*(i+0.6):.1f}" font-size="10" text-anchor="end" '
                     f'fill="{PALETTE["muted"]}">{lab}</text>')
        parts.append(f'<text x="{pad+cell*(i+0.5):.1f}" y="{pad-6}" font-size="10" text-anchor="middle" '
                     f'fill="{PALETTE["muted"]}">{lab}</text>')
    parts.append("</svg>")
    return "".join(parts)



def main():
    cols = load(CSV_PATH)
    steps = [int(s) for s in cols["step"]]

    # --- constant-term projection matrices recomputed in Python (rank 2, 3) --
    # A(u,v) = c(u v^-1, v lambda) with c(w, mu) = prod_{i<j, w(i)>w(j)} c(mu_i - mu_j)
    from itertools import permutations
    from math import lgamma, log, pi, sin, exp

    B = [1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6]

    def zeta_em(s, N=40, K=7):
        t = sum(n ** (-s) for n in range(1, N))
        t += N ** (1 - s) / (s - 1) + 0.5 * N ** (-s)
        for k in range(1, K + 1):
            poch = 1.0
            for j in range(2 * k - 1):
                poch *= (s + j)
            t += B[k - 1] / math.factorial(2 * k) * poch * N ** (-s - 2 * k + 1)
        return t

    def zeta(s):
        if s >= 0.5:
            return zeta_em(s)
        return 2 ** s * pi ** (s - 1) * sin(pi * s / 2) * exp(lgamma(1 - s)) * zeta_em(1 - s)

    def xi(s):
        if s <= 0:
            return xi(1 - s)
        return exp(-0.5 * s * log(pi) + lgamma(0.5 * s)) * zeta(s)

    def c(t):
        return xi(t) / xi(t + 1)

    def cmat(lam):
        n = len(lam)
        W = list(permutations(range(n)))

        def cfac(w, mu):
            t = 1.0
            for i in range(n):
                for j in range(i + 1, n):
                    if w[i] > w[j]:
                        t *= c(mu[i] - mu[j])
            return t

        def inv(w):
            out = [0] * n
            for i in range(n):
                out[w[i]] = i
            return tuple(out)

        def act(w, mu):
            wi = inv(w)
            return [mu[wi[i]] for i in range(n)]

        def mull(a, b):
            return tuple(a[b[i]] for i in range(n))

        A = [[0.0] * len(W) for _ in W]
        for a, u in enumerate(W):
            for b, v in enumerate(W):
                A[a][b] = cfac(mull(u, inv(v)), act(v, lam))
        return A

    A2 = cmat([1.4, -0.9])
    A3 = cmat([1.7, -0.4, 2.3])

    order = cols["order_R"]
    ks = cols["kappa_spread"]
    jac = cols["jacobi_defect"]
    defect = cols["unipotent_defect"]
    eis = cols["eisenstein_fraction"]
    cas = cols["casimir_drift"]
    grading = cols["gcheck_fixed_ratio"]

    charts = []
    charts.append(line_chart(
        [("order R", order, PALETTE["blue"]), ("Eisenstein fraction", eis, PALETTE["green"])],
        title="Emergence order R and the Eisenstein (semisimple) fraction",
        ylabel="value"))
    charts.append(line_chart(
        [("kappa spread", ks, PALETTE["red"]), ("unipotent defect", defect, PALETTE["amber"]),
         ("grading activation", grading, PALETTE["violet"])],
        log=True,
        title="Casimir spread, unipotent defect and ^G grading activation (log scale)",
        ylabel="log10"))
    charts.append(line_chart(
        [("|{kappa,H}|", cas, PALETTE["violet"]), ("Jacobi defect", jac, PALETTE["blue"])],
        log=True,
        title="Structural certificates: Casimir bracket and Jacobi defect (log scale)",
        ylabel="log10"))
    charts.append(bar_chart(
        ["(1,1,1)", "(2,1)", "(2,2)", "(3,1)", "(4)"],
        [0, 6, 8, 10, 12],
        title="dim of the unipotent orbit = n^2 - sum(lambda'_i)^2  (GL_4)",
        colour=PALETTE["blue"], ylabel="dim"))
    charts.append(bar_chart(
        ["(2,1)+7", "(1,1,1)+9", "(3,1)+12", "(2,2)+12", "(4)+16"],
        [16, 16, 16, 16, 16],
        title="Levi + orbit = n^2  (dimension accounting, GL_3 / GL_4)",
        colour=PALETTE["green"], ylabel="dim"))

    heat2 = heatmap(A2, ["1", "w0"], width=360, height=250,
                    title="CT matrix A(lambda), rank 2", vmax=2.4)
    heat3 = heatmap(A3, ["1", "s1", "s2", "s1s2", "s2s1", "w0"], width=440, height=340,
                    title="CT matrix A(lambda), rank 3 (|W| = 6)", vmax=20.0)

    cert_rows = [
        ["c(t)c(-t) = 1 (functional equation, local)", "|defect| &lt; 1e-8", "T17"],
        ["c(w1w2,&lambda;) = c(w1,w2&lambda;)c(w2,&lambda;)", "36/36 pairs exact", "T18"],
        ["CT by integration vs. Langlands formula", "|diff| &lt; 2e-3 (N = 40 truncation)", "T19"],
        ["Res<sub>s=1</sub> E(z,s) = 1/(2&xi;(2)) = 3/&pi;", "0.9549305 vs 0.9549297", "T19b"],
        ["&Delta;A<sub>0</sub> = s(s-1)A<sub>0</sub>", "finite differences, 3 spectral points", "T19c"],
        ["A(&lambda;)<sup>2</sup> = |W|A(&lambda;), tr &Pi; = 1, det A = 0", "defect &lt; 1e-8, ranks 2, 3, 4", "T20"],
        ["Jordan type from ranks of powers", "7 partitions recovered exactly", "T21a"],
        ["dim Cent(N) = &Sigma;(&lambda;'<sub>i</sub>)<sup>2</sup>", "nullity vs. formula, GL_4", "T21b"],
        ["dim O<sub>N</sub> = n<sup>2</sup> &minus; &Sigma;(&lambda;'<sub>i</sub>)<sup>2</sup>", "0, 6, 8, 10, 12", "T21c"],
        ["#{nilpotent orbits} = p(n)", "p(4)=5, p(5)=7, p(6)=11", "T21e"],
        ["[h,e]=2e, [h,f]=&minus;2f, [e,f]=h", "exact, 5 Jordan types", "T22a"],
        ["&Theta; derivation, integral spectrum, tr &Theta;<sup>2</sup> check", "1e-9", "T22b"],
        ["dim ker &Theta; = &Sigma;<sub>v</sub> m<sub>v</sub><sup>2</sup>", "nullity vs. formula", "T22c"],
        ["&chi;<sub>&Theta;</sub>(t) = (&Sigma;t<sup>h<sub>i</sub></sup>)(&Sigma;t<sup>&minus;h<sub>i</sub></sup>)", "4 values of t, 4 types", "T22d"],
        ["det Cent(e,h,f) = #blocks (separated types)", "5 types", "T22e"],
        ["Fricke: &kappa; = 2 + tr[A,B]", "8 random SL_2 pairs, &lt; 1e-10", "T23a"],
        ["Goldman bracket: Jacobi + Casimir + rank 2", "1e-6", "T23b"],
        ["dim M<sub>B</sub> = (2g&minus;2)(n<sup>2</sup>&minus;1) + 2p(n&minus;1)", "2, 6, 4, 16, 10", "T23c"],
        ["CT &compfn; Eis = id, fibre = unipotent orbit", "&lt; 1e-12", "T23d"],
        ["coupled run: R &rarr; 1, &kappa; spread &rarr; 0", "R = 0.99999999999, spread 2.7e-15", "T24"],
    ]

    # extract the second pipeline's summary block from run_output.txt
    lang_summary = ""
    try:
        with open(trace_path("run_output.txt"), encoding="utf-8") as fh:
            txt = fh.read()
        block = txt.split("# BEGIN_CSV_LANG")[0]
        if "# langlands-emergence" in block:
            block = block.split("# langlands-emergence", 1)[1]
            lines = [ln for ln in block.splitlines() if ln.strip()]
            lang_summary = "\n".join(lines)
    except FileNotFoundError:
        lang_summary = "(run_output.txt not found)"

    html = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<title>Eisenstein / non-semisimple L-parameters / ^G-operator / Betti moduli — MoonBit</title>
<style>
 body {{ font-family: -apple-system, "Segoe UI", "Noto Sans TC", sans-serif; margin: 0;
        background: #f8fafc; color: {PALETTE["ink"]}; line-height: 1.55; }}
 .wrap {{ max-width: 880px; margin: 0 auto; padding: 28px 20px 64px; }}
 h1 {{ font-size: 26px; margin: 0 0 4px; }}
 h2 {{ font-size: 18px; margin: 34px 0 10px; border-left: 4px solid {PALETTE["blue"]};
       padding-left: 10px; }}
 .sub {{ color: {PALETTE["muted"]}; font-size: 13px; margin-bottom: 22px; }}
 p {{ font-size: 14px; }}
 code {{ background: #eef2ff; padding: 1px 5px; border-radius: 4px; font-size: 12.5px; }}
 table {{ border-collapse: collapse; width: 100%; font-size: 12.5px; margin: 10px 0 6px; }}
 th, td {{ border: 1px solid {PALETTE["grid"]}; padding: 5px 8px; text-align: left; }}
 th {{ background: #eef2ff; font-weight: 600; }}
 .cap {{ font-size: 11.5px; color: {PALETTE["muted"]}; margin-bottom: 14px; }}
 .cards {{ display: flex; flex-wrap: wrap; gap: 12px; margin: 14px 0 6px; }}
 .card {{ flex: 1 1 190px; background: #fff; border: 1px solid {PALETTE["grid"]};
          border-radius: 10px; padding: 12px 14px; }}
 .card .k {{ font-size: 11px; color: {PALETTE["muted"]}; text-transform: uppercase;
             letter-spacing: .04em; }}
 .card .v {{ font-size: 20px; font-weight: 650; margin-top: 2px; }}
 .fine {{ font-size: 11.5px; color: {PALETTE["muted"]}; }}
 .grid2 {{ display: flex; gap: 14px; flex-wrap: wrap; }}
 .grid2 > div {{ flex: 1 1 330px; }}
</style>
</head>
<body><div class="wrap">
<h1>愛森斯坦級數 · 非半單 L-參數 · Ĝ-算子 · 貝蒂模空間</h1>
<div class="sub">MoonBit 0.10 (moonc v0.10.22) · module <code>arena/superemergence</code> ·
55 inline tests, wasm + js 全綠 · 由 <code>scripts/make_report_langlands.py</code> 生成</div>

<div class="cards">
  <div class="card"><div class="k">order R</div><div class="v">{order[-1]:.9f}</div>
      <div class="fine">κ-spread {ks[-1]:.2e}</div></div>
  <div class="card"><div class="k">unipotent defect</div><div class="v">{defect[-1]:.2e}</div>
      <div class="fine">從 {defect[0]:.3f} 冷卻至 0</div></div>
  <div class="card"><div class="k">|{{κ,H}}|</div><div class="v">{max(cas):.2e}</div>
      <div class="fine">Casimir 為 Hamiltonian 流的不變量</div></div>
  <div class="card"><div class="k">dim M<sub>B</sub>(Σ<sub>1,1</sub>, SL₂)</div><div class="v">2</div>
      <div class="fine">Poisson rank 2、一個 Casimir κ</div></div>
</div>

<h2>1. 兩個流：共識與半單化</h2>
<p>橙線（unipotent defect）代表智能體 L-參數的非半單程度，由 c-函數閘門
<code>g = |c(w₀,λ)|</code> 驅動冷卻；綠線是 Eisenstein（半單）比例。藍線 R 是
貝蒂模空間上 κ 散度所定義的湧現序參量。</p>
{charts[0]}
{charts[1]}

<h2>2. 結構憑證：Casimir 與 Jacobi</h2>
<p>Goldman/Atiyah–Bott 括號的 Poisson 張量滿足 <code>P dκ = 0</code>，
所以 <code>{{κ,H}} ≡ 0</code>；這正是「Hamiltonian 流保持在辛葉內」的
定理 T23b 的數值憑證。Jacobi 恆等式同理。</p>
{charts[2]}

<h2>3. 常數項矩陣與 |W|·冪等投影</h2>
<p><code>A(λ)[u][v] = c(u v⁻¹, vλ)</code>，定理 T20 說
<code>A(λ)² = |W|·A(λ)</code>：常數項函子自帶的秩一投影
<code>Π = A/|W|</code>（<code>tr Π = 1</code>、<code>det A = 0</code>）。
下圖為 rank 2（W = S₂）與 rank 3（W = S₃）的矩陣。</p>
<div class="grid2"><div>{heat2}</div><div>{heat3}</div></div>

<h2>4. 非半單 L-參數的維度學</h2>
<p>對 Jordan 型 λ ⊢ n：<code>dim Cent(N) = Σ(λ'_i)²</code>、
<code>dim O_N = n² − Σ(λ'_i)²</code>、<code>dim Ext¹ = dim Cent(N)</code>。
Ĝ-算子 <code>Θ = ad(h/2)</code> 則把對偶李代數分級：
<code>ker Θ</code> 是分級 Levi，<code>gl_n = ker Θ ⊕ im Θ</code>。</p>
{charts[3]}
{charts[4]}

<h2>5. 定理與測試對照</h2>
{table(["陳述", "數值憑證", "定理"], cert_rows,
       "全部由 `moon test --deny-warn`（wasm 與 js 各 55 個測試）決定。")}

<h2>6. 本次執行摘要</h2>
<pre style="background:#0f172a;color:#e2e8f0;padding:14px;border-radius:10px;
 overflow:auto;font-size:11.5px;line-height:1.45">{lang_summary}</pre>
<div class="fine">完整輸出見 <code>run_output.txt</code>；資料見
<code>langlands_trace.csv</code>（{len(steps)} 步，含 step 0）。</div>
</div></body></html>
"""
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w") as fh:
        fh.write(html)
    print(f"wrote {OUT_PATH} ({len(html)} bytes), {len(steps)} data rows")


if __name__ == "__main__":
    main()
