// ============================================================================
// 39_shared_vocab_fusion.mbt -- 通道十一: 同詞表 token 級融合
//                               (Round 7 / Node 5, T81-T83, A65-A67)
//
// ABSTRACTION.  Several language models that share *one* token vocabulary can
// be decoded as a *single* system: at every step each member contributes a
// distribution on the common label space, and one token is sampled from a
// fused distribution.  The classical facts are:
//
//   P (product / log-linear)   p(x) ∝ Π_i p_i(x)^{w_i},  Σ w_i = 1
//   M (mixture / linear)       p(x) = Σ_i w_i p_i(x)
//   A (assisted / speculative) p(x) ∝ [p_s(x) − p_w(x)]_+ + p_w(x)
//
// COMPUTER ABSTRACTION.  Only the *log*-distributions matter, so a member is a
// vector `ℓ ∈ (ℝ ∪ {−∞})^n`; `−∞` is a hard zero and every law is written so
// that `0·log 0 = 0` (the mask is the whole content of T82).  Cycling is
// concatenation of the per-cycle token streams (T83).
//
//   T81 (product idempotence and hard veto).  If every member gives the same
//       distribution p, then P(p, …, p) = p (weights sum to 1).  If one member
//       assigns zero mass to x, then P assigns zero mass to x.  Hence the
//       product rule is the intersection of the members' *supports* — this is
//       why the fused vocabulary is the shared vocabulary: any id missing from
//       one member's vocabulary becomes a permanently forbidden token.
//       (Operationally: those ids must be merged into the shared slot BEFORE
//       fusing — see `alias_merge`.)
//   T82 (Jensen–Shannon is mask-safe).  JS(p₁,…,p_k) = (1/k)Σ D(p_i ‖ p̄) is
//       finite on every member whose support is non-empty, including the
//       positions where some other member is 0; `p = 0` contributes exactly 0.
//       The implementation masks `0·log 0` instead of relying on IEEE `nan`.
//   T83 (cycling is associative).  Concatenating cycles reproduces the
//       per-cycle streams; the token budget is additive, and sampling seeds are
//       derived from the cycle index by `(seed + cycle) mod m`, so the stream
//       depends only on (seed, #cycles) and never on wall-clock time.
//
// The point of certifying this in MoonBit (Integer/Double, no GPU) is that the
// *numerics of the fusion law* are separated from the *model server*: the
// Python implementation in `space/fusion.py` is line-for-line the same three
// laws, so a MoonBit test failure here means the Python fusion is wrong — and
// the GPU code never has to be trusted blindly.
// ============================================================================

// ---------------------------------------------------------------------------
// 1. Masked log-distributions.
// ---------------------------------------------------------------------------

///|
/// `−∞` stands for a hard zero (a token outside a member's support).
pub fn fuse_neg_inf() -> Double {
  -1.0e308
}

///|
/// Softmax of a logit vector, returning logs.  Implemented as `x − max x − log Σ
/// exp(x − max x)`; entries equal to `−∞` (or ≤ `fuse_neg_inf()`) map to
/// `fuse_neg_inf()` *exactly*, so that downstream masks never see `nan`.
pub fn fuse_log_softmax(logits : Array[Double]) -> Array[Double] {
  let mut mx = fuse_neg_inf()
  for v in logits {
    if v > mx {
      mx = v
    }
  }
  let mut total = 0.0
  for v in logits {
    if v > fuse_neg_inf() {
      total += @math.exp(v - mx)
    }
  }
  let logz = @math.ln(total)
  Array::makei(logits.length(), fn(i) {
    if logits[i] <= fuse_neg_inf() {
      fuse_neg_inf()
    } else {
      logits[i] - mx - logz
    }
  })
}

///|
/// Is the entry a hard zero?
pub fn fuse_is_zero(lp : Double) -> Bool {
  lp <= fuse_neg_inf() + 1.0
}

///|
/// `w` as a log-weight, clamped away from 0.
pub fn fuse_log_weight(w : Double) -> Double {
  @math.ln(if w < 1.0e-12 { 1.0e-12 } else { w })
}

// ---------------------------------------------------------------------------
// 2. The three fusion laws (T81).
// ---------------------------------------------------------------------------

///|
/// `product`: `p ∝ Π_i p_i^{w_i}` — the log-sum, then normalise.
pub fn fuse_product(
  logprobs : Array[Array[Double]],
  weights : Array[Double],
) -> Array[Double] {
  let n = logprobs[0].length()
  let acc = Array::make(n, fuse_neg_inf())
  let veto = Array::make(n, false)
  for i in 0..<logprobs.length() {
    let w = weights[i]
    for k in 0..<n {
      let lp = logprobs[i][k]
      if fuse_is_zero(lp) {
        // 乘積律是**硬否決**：任何一個成員在 k 給 0，融合後就是 0。
        // （只把這一項「跳過」是錯的——其他成員仍會貢獻質量，等於把否決權吃掉。）
        veto[k] = true
      } else {
        // 權重是**指數**（w_i·log p_i），不是對數空間的加項；兩者只在 w_i = 1 時相同，
        // 所以 Σ w_i = 1 的等權情形下寫錯會立刻被 T81 的 idempotence 抓到。
        let term = lp * w
        if acc[k] <= fuse_neg_inf() {
          acc[k] = term
        } else {
          acc[k] += term
        }
      }
    }
  }
  for k in 0..<n {
    if veto[k] {
      acc[k] = fuse_neg_inf()
    }
  }
  fuse_log_softmax(acc)
}

///|
/// `mixture`: `p = Σ_i w_i p_i` — log-sum-exp over the members.
pub fn fuse_mixture(
  logprobs : Array[Array[Double]],
  weights : Array[Double],
) -> Array[Double] {
  let n = logprobs[0].length()
  let acc = Array::make(n, fuse_neg_inf())
  for i in 0..<logprobs.length() {
    let lw = fuse_log_weight(weights[i])
    for k in 0..<n {
      let lp = logprobs[i][k]
      if !fuse_is_zero(lp) {
        let term = lp + lw
        acc[k] = if acc[k] <= fuse_neg_inf() {
          term
        } else {
          let a = acc[k]
          let m = if a > term { a } else { term }
          m + @math.ln(@math.exp(a - m) + @math.exp(term - m))
        }
      }
    }
  }
  fuse_log_softmax(acc)
}

///|
/// `assisted` (speculative decoding as a distribution): `p ∝ [p_s − p_w]_+ +
/// p_w`.  Both inputs are log-probabilities; the output is a log-probability.
pub fn fuse_assisted(
  strong : Array[Double],
  weak : Array[Double],
) -> Array[Double] {
  let n = strong.length()
  let raw = Array::make(n, 0.0)
  let mut total = 0.0
  for k in 0..<n {
    let ps = if fuse_is_zero(strong[k]) { 0.0 } else { @math.exp(strong[k]) }
    let pw = if fuse_is_zero(weak[k]) { 0.0 } else { @math.exp(weak[k]) }
    let diff = ps - pw
    let v = (if diff > 0.0 { diff } else { 0.0 }) + pw
    raw[k] = v
    total += v
  }
  if total <= 0.0 {
    return fuse_log_softmax(raw)
  }
  Array::makei(n, fn(k) {
    if raw[k] <= 0.0 {
      fuse_neg_inf()
    } else {
      @math.ln(raw[k]) - @math.ln(total)
    }
  })
}

///|
/// The three-law dispatcher: `0 = product`, `1 = mixture`, `2 = assisted`.
pub fn fuse_apply(
  logprobs : Array[Array[Double]],
  weights : Array[Double],
  rule : Int32,
  strong~ : Int32 = 0,
  weak~ : Int32 = 1,
) -> Array[Double] {
  match rule {
    0 => fuse_product(logprobs, weights)
    1 => fuse_mixture(logprobs, weights)
    2 => fuse_assisted(logprobs[strong], logprobs[weak])
    _ => fuse_product(logprobs, weights)
  }
}

///|
/// `p(x)` from a log-distribution (0 for hard zeros).
pub fn fuse_prob(lp : Double) -> Double {
  if fuse_is_zero(lp) {
    0.0
  } else {
    @math.exp(lp)
  }
}

///|
/// Total mass of a log-distribution (should be 1 up to rounding).
pub fn fuse_total(lp : Array[Double]) -> Double {
  let mut s = 0.0
  for v in lp {
    s += fuse_prob(v)
  }
  s
}

///|
/// `argmax` with `−∞` entries ignored (`−1` if the support is empty).
pub fn fuse_argmax(lp : Array[Double]) -> Int32 {
  let mut best = -1
  let mut bv = fuse_neg_inf()
  for i in 0..<lp.length() {
    if lp[i] > bv {
      bv = lp[i]
      best = i
    }
  }
  best
}

// ---------------------------------------------------------------------------
// 3. T82: mask-safe Jensen–Shannon.
// ---------------------------------------------------------------------------

///|
/// `JS(p₁..p_k) = (1/k) Σ D(p_i ‖ p̄)` in nats; `p = 0` contributes exactly 0.
pub fn fuse_js(logprobs : Array[Array[Double]]) -> Double {
  let k = logprobs.length()
  if k == 0 {
    return 0.0
  }
  let n = logprobs[0].length()
  let mean = Array::make(n, 0.0)
  for i in 0..<k {
    for j in 0..<n {
      mean[j] += fuse_prob(logprobs[i][j]) / k.to_double()
    }
  }
  let mut total = 0.0
  for i in 0..<k {
    for j in 0..<n {
      let p = fuse_prob(logprobs[i][j])
      if p > 0.0 {
        let q = if mean[j] < 1.0e-300 { 1.0e-300 } else { mean[j] }
        total += p * (@math.ln(p) - @math.ln(q))
      }
    }
  }
  total / k.to_double()
}

///|
/// Top-1 agreement: are all members' argmax the same token?
pub fn fuse_agreement(logprobs : Array[Array[Double]]) -> Bool {
  if logprobs.length() == 0 {
    return true
  }
  let a = fuse_argmax(logprobs[0])
  let mut same = true
  for i in 1..<logprobs.length() {
    if fuse_argmax(logprobs[i]) != a {
      same = false
    }
  }
  same
}

// ---------------------------------------------------------------------------
// 4. T83: cycling.
// ---------------------------------------------------------------------------

///|
/// A per-cycle record: how many common tokens it produced, and its seed.
pub(all) struct FuseCycle {
  index : Int32
  tokens : Int32
  seed : Int32
}

///|
/// The sampling seed of cycle `c`: `(seed + c) mod m`, with `m` the modulus
/// (2^31 − 1 in the Python implementation).  Purely a function of its inputs,
/// so a 3-cycle run is reproducible.
pub fn fuse_cycle_seed(seed : Int32, cycle : Int32, m : Int32) -> Int32 {
  // 用 Int64 做加法：seed 與 cycle 都 < m 時，Int32 的 s + c 會溢出（2147483646 + 2）
  let mm = m.to_int64()
  let s = (seed.to_int64() % mm + mm) % mm
  let c = (cycle.to_int64() % mm + mm) % mm
  ((s + c) % mm).to_int()
}

///|
/// T83: concatenate the per-cycle token streams.  Returns the joined stream and
/// its length (the budget is additive).
pub fn fuse_concat_cycles(
  cycles : Array[Array[Int32]],
) -> (Array[Int32], Int32) {
  let out : Array[Int32] = []
  for c in cycles {
    for t in c {
      out.push(t)
    }
  }
  (out, out.length())
}

///|
/// A whole multi-cycle run: the per-cycle records plus the joined stream.
pub(all) struct FuseRun {
  cycles : Array[FuseCycle]
  stream : Array[Int32]
  seed : Int32
  modulus : Int32
}

///|
/// Build a run from already-computed per-cycle streams (the GPU side computes
/// the streams; this is the bookkeeping that T83 certifies).
pub fn fuse_build_run(
  streams : Array[Array[Int32]],
  seed~ : Int32,
  modulus~ : Int32 = 2147483647,
) -> FuseRun {
  let cycles : Array[FuseCycle] = []
  for i in 0..<streams.length() {
    cycles.push(FuseCycle::{
      index: i,
      tokens: streams[i].length(),
      seed: fuse_cycle_seed(seed, i, modulus),
    })
  }
  let (stream, _) = fuse_concat_cycles(streams)
  FuseRun::{ cycles, stream, seed, modulus, }
}

// ---------------------------------------------------------------------------
// 5. Weights and the shared vocabulary (T81's operational corollary).
// ---------------------------------------------------------------------------

///|
/// `[1/k, …, 1/k]`, or `[w, (1−w)/(k−1), …]` when a first-member weight is given
/// (`w ≤ 0` means "equal weights").
pub fn fuse_weights(k : Int32, w~ : Double = 0.0) -> Array[Double] {
  if k <= 0 {
    return []
  }
  if w <= 0.0 || k == 1 {
    return Array::make(k, 1.0 / k.to_double())
  }
  let first = if w > 1.0 { 1.0 } else { w }
  let rest = (1.0 - first) / (k - 1).to_double()
  let out = Array::make(k, rest)
  out[0] = first
  out
}

///|
/// T81 corollary: an id that no member can emit is dead weight; an id that only
/// *some* members can emit becomes forbidden under the product rule.  The fix
/// is to merge alias ids into one shared slot before fusing.
/// Returns the number of ids that survive (are emitted by at least one member)
/// and the number that survive *all* members.
pub fn fuse_vocab_coverage(supports : Array[Array[Bool]]) -> (Int32, Int32) {
  if supports.length() == 0 {
    return (0, 0)
  }
  let n = supports[0].length()
  let any = Array::make(n, false)
  let all = Array::make(n, true)
  for s in supports {
    for i in 0..<n {
      if s[i] {
        any[i] = true
      } else {
        all[i] = false
      }
    }
  }
  let mut a = 0
  let mut b = 0
  for i in 0..<n {
    if any[i] {
      a += 1
    }
    if all[i] {
      b += 1
    }
  }
  (a, b)
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T81: product is idempotent on identical members and is a hard veto" {
  let p : Array[Double] = fuse_log_softmax([1.0, 2.0, 0.5, -1.0])
  let w = fuse_weights(3)
  // three identical members ⇒ the same distribution back
  let same = fuse_product([p, p, p], w)
  for i in 0..<p.length() {
    assert_true(
      approx(same[i], p[i], tol=1.0e-12),
      msg=fuse_prob(same[i]).to_string(),
    )
  }
  // one member with a hard zero ⇒ the fused distribution has a hard zero
  let vetoed = p.copy()
  vetoed[1] = fuse_neg_inf()
  let fused = fuse_product([p, vetoed], [0.5, 0.5])
  assert_true(
    fuse_prob(fused[1]) == 0.0,
    msg="veto cell = \{fuse_prob(fused[1])}",
  )
  assert_true(fuse_prob(fused[0]) > 0.0)
  assert_true(approx(fuse_total(fused), 1.0, tol=1.0e-12))
  // the mixture rule does *not* veto
  let mixed = fuse_mixture([p, vetoed], [0.5, 0.5])
  assert_true(fuse_prob(mixed[1]) > 0.1)
}

///|
test "T81: assisted keeps the strong model's peak and is a distribution" {
  let strong = fuse_log_softmax([9.0, 1.0, 0.0])
  let weak = fuse_log_softmax([0.0, 6.0, 1.0])
  let out = fuse_assisted(strong, weak)
  assert_true(approx(fuse_total(out), 1.0, tol=1.0e-12))
  assert_eq(fuse_argmax(out), 0) // the strong peak wins
  // identical members ⇒ identity
  let same = fuse_assisted(strong, strong)
  for i in 0..<strong.length() {
    assert_true(approx(same[i], strong[i], tol=1.0e-12))
  }
  // hand value: residual = [0.8, 0, 0], Z = 1.8 ⇒ p = [0.9, 0.6, 0.3]/1.8
  let s2 = fuse_log_softmax([2.302585, 0.0, 0.0]) // [0.9, 0.05, 0.05] up to 1e-3
  let w2 = fuse_log_softmax([0.0, 1.609438, 1.098612])
  let o2 = fuse_assisted(s2, w2)
  assert_true(approx(fuse_prob(o2[0]) / fuse_prob(o2[1]), 1.5, tol=1.0e-3))
}

///|
test "T81: the three-law dispatcher agrees with the direct laws" {
  let a : Array[Double] = fuse_log_softmax([2.0, 1.0, 0.0])
  let b : Array[Double] = fuse_log_softmax([0.0, 1.5, 0.5])
  let w = fuse_weights(2)
  let members = [a, b]
  let by_dispatch = [
    fuse_apply(members, w, 0),
    fuse_apply(members, w, 1),
    fuse_apply(members, w, 2),
  ]
  let direct = [
    fuse_product(members, w),
    fuse_mixture(members, w),
    fuse_assisted(a, b),
  ]
  for r in 0..<3 {
    for k in 0..<3 {
      assert_true(approx(by_dispatch[r][k], direct[r][k], tol=1.0e-12))
    }
    assert_true(approx(fuse_total(by_dispatch[r]), 1.0, tol=1.0e-12))
  }
  // unknown rule falls back to product
  let fallback = fuse_apply(members, w, 99)
  assert_true(approx(fallback[0], direct[0][0], tol=1.0e-12))
}

///|
test "T82: JS is finite (mask-safe) and zero on identical members" {
  let a : Array[Double] = fuse_log_softmax([1.0, 1.0, 1.0])
  let b = a.copy()
  b[2] = fuse_neg_inf()
  let js = fuse_js([a, b])
  assert_true(js > 0.0 && js < 1.0, msg=js.to_string())
  // identical ⇒ exactly 0
  assert_true(fuse_js([a, a, a]) < 1.0e-12)
  // the fully-zero member is the extreme case: still finite
  let zero = Array::make(3, fuse_neg_inf())
  let js2 = fuse_js([a, zero])
  assert_true(js2 < 10.0, msg=js2.to_string())
  // a member that is 0 everywhere must not poison JS of the others
  assert_true(fuse_js([a, a, zero]) < fuse_js([a, zero]) + 1.0)
  // agreement diagnostics
  assert_true(fuse_agreement([a, a.copy()]))
  assert_true(
    !fuse_agreement([fuse_log_softmax([1.0, 0.0]), fuse_log_softmax([0.0, 1.0])]),
  )
}

///|
test "T83: cycles concatenate and seeds are reproducible" {
  let run = fuse_build_run([[5, 6, 7], [8, 9], [10]], seed=0)
  assert_eq(run.cycles.length(), 3)
  assert_eq(run.stream.length(), 6)
  assert_eq(run.cycles[0].tokens, 3)
  assert_eq(run.cycles[1].tokens, 2)
  assert_eq(run.cycles[2].tokens, 1)
  // the 3-cycle stream is the concatenation, in order
  let expect = [5, 6, 7, 8, 9, 10]
  for i in 0..<6 {
    assert_eq(run.stream[i], expect[i])
  }
  // seeds: (seed + cycle) mod m, deterministic
  assert_eq(run.cycles[0].seed, 0)
  assert_eq(run.cycles[1].seed, 1)
  assert_eq(run.cycles[2].seed, 2)
  let r2 = fuse_build_run([[5, 6, 7], [8, 9], [10]], seed=0)
  assert_eq(r2.cycles[2].seed, run.cycles[2].seed)
  assert_eq(fuse_cycle_seed(2147483646, 2, 2147483647), 1)
}

///|
test "T81 corollary: coverage of the shared vocabulary" {
  // id 0 and 1 in both members, id 2 only in the first, id 3 in neither
  let (any, all) = fuse_vocab_coverage([
    [true, true, true, false],
    [true, true, false, false],
  ])
  assert_eq(any, 3)
  assert_eq(all, 2)
  let w = fuse_weights(3, w=0.5)
  assert_true(approx(w[0], 0.5, tol=1.0e-12))
  assert_true(approx(w[1], 0.25, tol=1.0e-12))
  assert_true(approx(w[0] + w[1] + w[2], 1.0, tol=1.0e-12))
}
