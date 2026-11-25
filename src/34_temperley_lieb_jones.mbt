// ============================================================================
// 34_temperley_lieb_jones.mbt -- 通道七: 量子鏈 / 辮群 / 瓊斯多項式
//                                (Round 7 / Node 1, T69-T71, A53-A55)
//
// ABSTRACTION.  The Jones polynomial is the Markov-trace functional on the
// Temperley–Lieb algebra:
//
//     TL_n(δ) = ⟨ 1, e_1, ..., e_{n−1} |  e_i² = δ e_i,
//                                          e_i e_j = e_j e_i  (|i−j| ≥ 2),
//                                          e_i e_{i±1} e_i = e_i ⟩
//     δ = −A² − A^{−2}
//
// The algebra is *finite dimensional*: its basis is the set of non-crossing
// perfect matchings on 2n points (Catalan many).  Composition = stacking +
// replacing closed circles by δ (the *matching model*); the Markov trace is
//     tr(m) = δ^{closure(m) − 1}
// and the Kauffman bracket of a braid closure is
//     ⟨β̂⟩ = tr(π(β)),      π(σ_i^{±1}) = A^{±1}·1 + A^{∓1}·e_i
// so that the (unoriented, framing-corrected) Jones polynomial is
//     V_β(t) = (−A³)^{−w(β)} ⟨β̂⟩,      t = A^{−4},   w = exponent sum.
//
// COMPUTER ABSTRACTION (a law, not a formula).  A *Laurent polynomial over ℤ*
// is a fixed window of Int32 coefficients; the Temperley–Lieb element is a list
// of (matching, coefficient) pairs.  Every theorem below is a statement about
// exact integer arrays, hence decidable without any floating point.
//
//   T69  matching model = TL_n: e_i² = δe_i, e_i e_{i±1} e_i = e_i, and
//        π is an algebra map from the group algebra of the braid group
//        (σ_iσ_i^{−1} = 1 lives in the algebra, not in the group).
//   T70  Jones values: V(3₁) = −t^{−4} + t^{−3} + t^{−1},
//        V(4₁) = t^{−2} − t^{−1} + 1 − t + t², V(5₁) = ...
//        (computed, not quoted: the pipeline recomputes them from the braid
//        words by pure integer arithmetic.)
//   T71  two independent computations agree: the matching model and the
//        *skein recurrence* of the 2-strand torus family
//            x_{n+1} = t² x_{n−1} + (t^{1/2} − t^{−1/2}) t x_n ,
//            x_0 = −(t^{1/2} + t^{−1/2}),   x_1 = 1
//        give the same Laurent polynomial at every n (Certificate: defect 0).
//        Plus the classical invariants: V(1) = 1 and |V(−1)| = det(K).
//
// The channel is tied to the rest of the round by the *Hecke* end: TL_n is the
// q → 1 quotient of the Hecke algebra H_n(q) of module 35, and the Chebyshev
// recursion appearing in the trace is the same one that computes the Satake
// sym²/∧² recoupling of module 29.
// ============================================================================

// ---------------------------------------------------------------------------
// 1. Laurent polynomials over ℤ inside a fixed exponent window.
// ---------------------------------------------------------------------------

///|
/// Laurent polynomial `Σ c_e A^e` with `e ∈ [−64, 63]`, stored densely.
pub(all) struct Lz {
  e : Array[Int32]
}

///|
const LZ_OFF : Int32 = 96

///|
const LZ_LEN : Int32 = 192

///|
pub fn lz_zero() -> Lz {
  Lz::{ e: Array::make(LZ_LEN, 0), }
}

///|
pub fn lz_mono(exponent : Int32, c : Int32) -> Lz {
  let p = lz_zero()
  p.e[exponent + LZ_OFF] = c
  p
}

///|
pub fn lz_const(c : Int32) -> Lz {
  lz_mono(0, c)
}

///|
pub fn lz_coeff(p : Lz, exponent : Int32) -> Int32 {
  if exponent < -LZ_OFF || exponent >= LZ_LEN - LZ_OFF {
    0
  } else {
    p.e[exponent + LZ_OFF]
  }
}

///|
pub fn lz_add(a : Lz, b : Lz) -> Lz {
  Lz::{ e: Array::makei(LZ_LEN, fn(i) { a.e[i] + b.e[i] }), }
}

///|
pub fn lz_neg(a : Lz) -> Lz {
  Lz::{ e: Array::makei(LZ_LEN, fn(i) { -a.e[i] }), }
}

///|
pub fn lz_scale(a : Lz, k : Int32) -> Lz {
  Lz::{ e: Array::makei(LZ_LEN, fn(i) { a.e[i] * k }), }
}

///|
pub fn lz_mul(a : Lz, b : Lz) -> Lz {
  let out = lz_zero()
  // bounded convolution: only the occupied exponent ranges are touched
  let (alo, ahi) = lz_index_span(a)
  let (blo, bhi) = lz_index_span(b)
  if alo > ahi || blo > bhi {
    return out
  }
  for i in alo..<=ahi {
    if a.e[i] != 0 {
      for j in blo..<=bhi {
        if b.e[j] != 0 {
          let k = i + j - LZ_OFF
          if k >= 0 && k < LZ_LEN {
            out.e[k] += a.e[i] * b.e[j]
          }
        }
      }
    }
  }
  out
}

///|
/// The occupied *array-index* range `[lo, hi]` of a Laurent polynomial
/// (`lo > hi` for the zero polynomial).
pub fn lz_index_span(a : Lz) -> (Int32, Int32) {
  let mut lo = LZ_LEN
  let mut hi = -1
  for i in 0..<LZ_LEN {
    if a.e[i] != 0 {
      if i < lo {
        lo = i
      }
      hi = i
    }
  }
  (lo, hi)
}

///|
pub fn lz_is_zero(a : Lz) -> Bool {
  let mut z = true
  for i in 0..<LZ_LEN {
    if a.e[i] != 0 {
      z = false
    }
  }
  z
}

///|
/// Number of non-zero coefficients (the *breadth* of the polynomial).
pub fn lz_support(a : Lz) -> Int32 {
  let mut n = 0
  for i in 0..<LZ_LEN {
    if a.e[i] != 0 {
      n += 1
    }
  }
  n
}

///|
/// `(min exponent, max exponent)` of the support, or `(0, 0)` for zero.
pub fn lz_span(a : Lz) -> (Int32, Int32) {
  let mut lo = 0
  let mut hi = 0
  let mut first = true
  for i in 0..<LZ_LEN {
    if a.e[i] != 0 {
      if first {
        lo = i - LZ_OFF
        first = false
      }
      hi = i - LZ_OFF
    }
  }
  (lo, hi)
}

///|
/// `p(1)` — the sum of the coefficients.
pub fn lz_eval_one(a : Lz) -> Int32 {
  let mut s = 0
  for i in 0..<LZ_LEN {
    s += a.e[i]
  }
  s
}

///|
/// `p(−1)` — the alternating sum of the coefficients.
pub fn lz_eval_minus_one(a : Lz) -> Int32 {
  let mut s = 0
  for i in 0..<LZ_LEN {
    if (i - LZ_OFF) % 2 == 0 {
      s += a.e[i]
    } else {
      s -= a.e[i]
    }
  }
  s
}

///|
/// `p(1/t)·t^k` normalised so that the exponent window is preserved: the
/// *mirror* polynomial (used to certify chirality: a knot with `V(t) ≠ V(1/t)`
/// is chiral, as for the trefoil).
pub fn lz_reverse(a : Lz) -> Lz {
  Lz::{ e: Array::makei(LZ_LEN, fn(i) { a.e[LZ_LEN - 1 - i] }), }
}

///|
/// `a` re-indexed from `A`-powers to `t`-powers by `t = A^{−4}`: every exponent
/// must be divisible by 4.  Returns `(poly, defect)` where `defect` counts the
/// coefficients that are *not* on the lattice `4ℤ` (it is 0 exactly when the
/// braid closure is a genuine link polynomial).
pub fn lz_to_t(a : Lz) -> (Lz, Int32) {
  let out = lz_zero()
  let mut defect = 0
  for i in 0..<LZ_LEN {
    let e = i - LZ_OFF
    if a.e[i] != 0 {
      if e % 4 == 0 {
        out.e[e / 4 + LZ_OFF] += a.e[i]
      } else {
        defect += 1
      }
    }
  }
  (out, defect)
}

///|
/// Human-readable form, e.g. `-A^-4 + A^-3 + A^-1`.
pub fn lz_to_string(a : Lz) -> String {
  let mut out = ""
  let mut wrote = false
  for i in 0..<LZ_LEN {
    let c = a.e[i]
    if c != 0 {
      let e = i - LZ_OFF
      if wrote {
        if c > 0 {
          out = out + " + "
        } else {
          out = out + " - "
        }
      } else if c < 0 {
        out = out + "-"
      }
      let m = if c < 0 { -c } else { c }
      out = out + m.to_string()
      if e == 1 {
        out = out + "A"
      } else {
        out = out + "A^" + e.to_string()
      }
      wrote = true
    }
  }
  if !wrote {
    out = "0"
  }
  out
}

// ---------------------------------------------------------------------------
// 2. The Temperley–Lieb algebra in the non-crossing matching model (T69).
// ---------------------------------------------------------------------------

///|
/// A Temperley–Lieb element: parallel arrays of non-crossing matchings
/// (`2n` entries, the partner of each point) and Laurent coefficients.
pub(all) struct TLelt {
  parts : Array[Array[Int32]]
  coefs : Array[Lz]
}

///|
/// `δ = −A² − A^{−2}`.
pub fn tl_delta() -> Lz {
  lz_add(lz_mono(2, -1), lz_mono(-2, -1))
}

///|
/// `δ^k`.
pub fn tl_delta_pow(k : Int32) -> Lz {
  let mut out = lz_const(1)
  for _ in 0..<k {
    out = lz_mul(out, tl_delta())
  }
  out
}

///|
/// The empty element.
pub fn tl_zero() -> TLelt {
  TLelt::{ parts: [], coefs: [], }
}

///|
/// Add `c·m` to an element (linear scan: elements here have ≤ 14 terms).
pub fn tl_add_term(el : TLelt, m : Array[Int32], c : Lz) -> Unit {
  if lz_is_zero(c) {
    return
  }
  let mut idx = -1
  for i in 0..<el.parts.length() {
    if tl_same_matching(el.parts[i], m) {
      idx = i
    }
  }
  if idx >= 0 {
    el.coefs[idx] = lz_add(el.coefs[idx], c)
  } else {
    el.parts.push(m)
    el.coefs.push(c)
  }
}

///|
pub fn tl_same_matching(a : Array[Int32], b : Array[Int32]) -> Bool {
  let mut same = a.length() == b.length()
  if same {
    for i in 0..<a.length() {
      if a[i] != b[i] {
        same = false
      }
    }
  }
  same
}

///|
/// All perfect matchings of `2n` points, as partner arrays (`(2n−1)!!` many).
pub fn tl_matchings(n : Int32) -> Array[Array[Int32]] {
  let out : Array[Array[Int32]] = []
  let pts : Array[Int32] = []
  for i in 0..<(2 * n) {
    pts.push(i)
  }
  let cur : Array[Int32] = []
  tl_matchings_rec(pts, cur, n, out)
  out
}

///|
fn tl_matchings_rec(
  pts : Array[Int32],
  cur : Array[Int32],
  n : Int32,
  out : Array[Array[Int32]],
) -> Unit {
  if pts.length() == 0 {
    let m = Array::make(2 * n, -1)
    for i in 0..<(cur.length() / 2) {
      let a = cur[2 * i]
      let b = cur[2 * i + 1]
      m[a] = b
      m[b] = a
    }
    out.push(m)
    return
  }
  let a = pts[0]
  let rest : Array[Int32] = []
  for i in 1..<pts.length() {
    rest.push(pts[i])
  }
  for j in 0..<rest.length() {
    let b = rest[j]
    let nxt : Array[Int32] = []
    for i in 0..<rest.length() {
      if i != j {
        nxt.push(rest[i])
      }
    }
    cur.push(a)
    cur.push(b)
    tl_matchings_rec(nxt, cur, n, out)
    let _ = cur.pop()
    let _ = cur.pop()
  }
}

///|
/// The identity matching `{(k, n+k)}` — the identity of `TL_n`.
pub fn tl_identity_matching(n : Int32) -> Array[Int32] {
  let m = Array::make(2 * n, 0)
  for k in 0..<n {
    m[k] = n + k
    m[n + k] = k
  }
  m
}

///|
/// Is the matching non-crossing?  Two arcs `(a,b)`, `(c,d)` cross iff
/// `a < c < b < d` (or symmetrically): the Temperley–Lieb basis is exactly the
/// non-crossing ones.
pub fn tl_is_non_crossing(m : Array[Int32]) -> Bool {
  let n = m.length() / 2
  let mut ok = true
  for a in 0..<m.length() {
    let b = m[a]
    if b > a {
      for c in 0..<m.length() {
        let d = m[c]
        if d > c && c != a {
          if (a < c && c < b && b < d) || (c < a && a < d && d < b) {
            ok = false
          }
        }
      }
    }
  }
  let _ = n
  ok
}

///|
/// Only the non-crossing matchings — the basis of `TL_n` (Catalan many).
pub fn tl_non_crossing_matchings(n : Int32) -> Array[Array[Int32]] {
  let out : Array[Array[Int32]] = []
  for m in tl_matchings(n) {
    if tl_is_non_crossing(m) {
      out.push(m)
    }
  }
  out
}

///|
/// `e_i`: the cup-cap at the slot `(i, i+1)`.
pub fn tl_generator_matching(n : Int32, i : Int32) -> Array[Int32] {
  let m = Array::make(2 * n, -1)
  for k in 0..<(n - 1) {
    if k == i {
      m[k] = k + 1
      m[k + 1] = k
      m[n + k] = n + k + 1
      m[n + k + 1] = n + k
    }
  }
  for k in 0..<n {
    if m[k] < 0 {
      m[k] = n + k
      m[n + k] = k
    }
  }
  m
}

///|
pub fn tl_one(n : Int32) -> TLelt {
  let el = tl_zero()
  tl_add_term(el, tl_identity_matching(n), lz_const(1))
  el
}

///|
pub fn tl_gen(n : Int32, i : Int32) -> TLelt {
  let el = tl_zero()
  tl_add_term(el, tl_generator_matching(n, i), lz_const(1))
  el
}

///|
/// Union–find, flat array.
fn tl_uf_find(parent : Array[Int32], x : Int32) -> Int32 {
  let mut r = x
  while parent[r] != r {
    r = parent[r]
  }
  let mut c = x
  while parent[c] != c {
    let nxt = parent[c]
    parent[c] = r
    c = nxt
  }
  r
}

///|
/// Union two nodes by root (never creates a cycle).
fn tl_uf_union(parent : Array[Int32], a : Int32, b : Int32) -> Unit {
  let ra = tl_uf_find(parent, a)
  let rb = tl_uf_find(parent, b)
  if ra != rb {
    parent[ra] = rb
  }
}

///|
/// Stack `top` above `bottom` (both honest matchings of `2n` points) and
/// return the resulting matching plus the number of closed circles.
pub fn tl_compose(
  n : Int32,
  top : Array[Int32],
  bottom : Array[Int32],
) -> (Array[Int32], Int32) {
  // node ids: bottom-layer m_bottom bottoms 0..n−1, its tops n..2n−1,
  //           top-layer m_top bottoms 2n..3n−1, its tops 3n..4n−1.
  let parent = Array::make(4 * n, 0)
  for i in 0..<(4 * n) {
    parent[i] = i
  }
  for i in 0..<n {
    tl_uf_union(parent, n + i, 2 * n + i)
  }
  // arcs of both matchings (caps join two bottom points, cups two top points,
  // so both rows are scanned)
  for a in 0..<(2 * n) {
    tl_uf_union(parent, a, bottom[a])
  }
  for a in 0..<(2 * n) {
    tl_uf_union(parent, 2 * n + a, 2 * n + top[a])
  }
  // compress, then label components: each component gets (a, b) = its exterior
  // points, in increasing order; a component without exterior points is a
  // closed circle.
  let root_of = Array::make(4 * n, -1)
  for i in 0..<(4 * n) {
    root_of[i] = tl_uf_find(parent, i)
  }
  let roots : Array[Int32] = []
  let ext_a : Array[Int32] = []
  let ext_b : Array[Int32] = []
  let ext_n : Array[Int32] = []
  for i in 0..<(4 * n) {
    let r = root_of[i]
    let mut slot = -1
    for k in 0..<roots.length() {
      if roots[k] == r {
        slot = k
      }
    }
    if slot < 0 {
      roots.push(r)
      ext_a.push(-1)
      ext_b.push(-1)
      ext_n.push(0)
      slot = roots.length() - 1
    }
    let is_ext = i < n || i >= 3 * n
    if is_ext {
      if ext_n[slot] == 0 {
        ext_a[slot] = i
      } else if ext_n[slot] == 1 {
        ext_b[slot] = i
      }
      ext_n[slot] += 1
    }
  }
  let mut loops = 0
  for k in 0..<roots.length() {
    if ext_n[k] == 0 {
      loops += 1
    }
  }
  let res = Array::make(2 * n, -1)
  for k in 0..<roots.length() {
    let a0 = ext_a[k]
    let b0 = ext_b[k]
    if ext_n[k] == 2 && a0 >= 0 && b0 >= 0 {
      let a = if a0 < n { a0 } else { a0 - 2 * n }
      let b = if b0 < n { b0 } else { b0 - 2 * n }
      res[a] = b
      res[b] = a
    }
  }
  (res, loops)
}

///|
/// Product of Temperley–Tau elements: `x·y` = `x` stacked above `y`, closed
/// circles replaced by `δ`.
pub fn tl_mul(n : Int32, x : TLelt, y : TLelt) -> TLelt {
  let out = tl_zero()
  for i in 0..<x.parts.length() {
    for j in 0..<y.parts.length() {
      let (m, loops) = tl_compose(n, x.parts[i], y.parts[j])
      let c = lz_mul(lz_mul(x.coefs[i], y.coefs[j]), tl_delta_pow(loops))
      tl_add_term(out, m, c)
    }
  }
  tl_clean(out)
}

///|
/// Drop the terms whose coefficient has cancelled to zero (canonical form of a
/// Temperley–Lieb element; keeps the matching model a *basis* representation).
pub fn tl_clean(el : TLelt) -> TLelt {
  let parts : Array[Array[Int32]] = []
  let coefs : Array[Lz] = []
  for i in 0..<el.parts.length() {
    if !lz_is_zero(el.coefs[i]) {
      parts.push(el.parts[i])
      coefs.push(el.coefs[i])
    }
  }
  TLelt::{ parts, coefs, }
}

///|
/// Closed circles in the closure of the matching `m` (points `i` and `n+i` of
/// the same side are joined across the top and the bottom).
pub fn tl_closure_loops(n : Int32, m : Array[Int32]) -> Int32 {
  let parent = Array::make(2 * n, 0)
  for i in 0..<(2 * n) {
    parent[i] = i
  }
  // the matching's own arcs (both rows)
  for a in 0..<(2 * n) {
    tl_uf_union(parent, a, m[a])
  }
  // the closure strands i ↔ n+i
  for i in 0..<n {
    tl_uf_union(parent, i, n + i)
  }
  let roots : Array[Int32] = []
  for i in 0..<(2 * n) {
    let r = tl_uf_find(parent, i)
    let mut has = false
    for k in 0..<roots.length() {
      if roots[k] == r {
        has = true
      }
    }
    if !has {
      roots.push(r)
    }
  }
  // circles = components − 1 (the one non-contractible component is the loop
  // of the closure, not a circle)
  roots.length() - 1
}

///|
/// The *geometric* Markov trace `tr(m) = δ^{circles(m)}` where `circles` counts
/// the closed circles of the closure (the non-contractible loop excluded).
pub fn tl_trace(n : Int32, el : TLelt) -> Lz {
  let mut out = lz_zero()
  for i in 0..<el.parts.length() {
    let loops = tl_closure_loops(n, el.parts[i])
    out = lz_add(out, lz_mul(el.coefs[i], tl_delta_pow(loops)))
  }
  out
}

///|
/// `π(σ_i) = A·1 + A^{−1}·e_i` and `π(σ_i^{−1}) = A^{−1}·1 + A·e_i`.
pub fn tl_braid_generator(n : Int32, i : Int32, exponent : Int32) -> TLelt {
  let el = tl_zero()
  if exponent > 0 {
    tl_add_term(el, tl_identity_matching(n), lz_mono(1, 1))
    tl_add_term(el, tl_generator_matching(n, i), lz_mono(-1, 1))
  } else {
    tl_add_term(el, tl_identity_matching(n), lz_mono(-1, 1))
    tl_add_term(el, tl_generator_matching(n, i), lz_mono(1, 1))
  }
  el
}

///|
/// `π(β)` for a braid word: `Array[(generator index, exponent ±1)]`.
pub fn tl_braid_element(n : Int32, word : Array[(Int32, Int32)]) -> TLelt {
  let mut el = tl_one(n)
  for step in word {
    let (i, ex) = step
    el = tl_mul(n, el, tl_braid_generator(n, i, ex))
  }
  el
}

///|
/// The Kauffman bracket of the braid closure, `⟨β̂⟩ = tr(π(β))` (in `A`).
pub fn kauffman_bracket(n : Int32, word : Array[(Int32, Int32)]) -> Lz {
  tl_trace(n, tl_braid_element(n, word))
}

///|
/// The exponent sum `w(β) = Σ ±1`.
pub fn braid_writhe(word : Array[(Int32, Int32)]) -> Int32 {
  let mut w = 0
  for step in word {
    w += step.1
  }
  w
}

///|
/// T70: the Jones polynomial `V_β(t) = (−A³)^{−w} ⟨β̂⟩` with `t = A^{−4}`,
/// returned in `t`-powers together with the lattice defect (0 for a genuine
/// link, i.e. when every `A`-exponent is a multiple of 4).
pub fn jones_polynomial(n : Int32, word : Array[(Int32, Int32)]) -> (Lz, Int32) {
  let br = kauffman_bracket(n, word)
  let w = braid_writhe(word)
  // (−A³)^{−w} = (−1)^w A^{−3w}
  let sign = if w % 2 == 0 { 1 } else { -1 }
  let norm = lz_mul(br, lz_mono(-3 * w, sign))
  lz_to_t(norm)
}

///|
/// Three strand *half*-integer Laplacian, used only by the skein certificate:
/// coefficients are indexed by twice the exponent of `t^{1/2}`.
pub(all) struct Hlz {
  e : Array[Int32]
}

///|
pub fn hlz_zero() -> Hlz {
  Hlz::{ e: Array::make(LZ_LEN, 0), }
}

///|
/// `c·t^{k/2}`.
pub fn hlz_mono(k2 : Int32, c : Int32) -> Hlz {
  let p = hlz_zero()
  p.e[k2 + LZ_OFF] = c
  p
}

///|
pub fn hlz_add(a : Hlz, b : Hlz) -> Hlz {
  Hlz::{ e: Array::makei(LZ_LEN, fn(i) { a.e[i] + b.e[i] }), }
}

///|
pub fn hlz_mul(a : Hlz, b : Hlz) -> Hlz {
  let out = hlz_zero()
  let (alo, ahi) = hlz_index_span(a)
  let (blo, bhi) = hlz_index_span(b)
  if alo > ahi || blo > bhi {
    return out
  }
  for i in alo..<=ahi {
    if a.e[i] != 0 {
      for j in blo..<=bhi {
        if b.e[j] != 0 {
          let k = i + j - LZ_OFF
          if k >= 0 && k < LZ_LEN {
            out.e[k] += a.e[i] * b.e[j]
          }
        }
      }
    }
  }
  out
}

///|
/// The occupied index range of a half-power polynomial.
pub fn hlz_index_span(a : Hlz) -> (Int32, Int32) {
  let mut lo = LZ_LEN
  let mut hi = -1
  for i in 0..<LZ_LEN {
    if a.e[i] != 0 {
      if i < lo {
        lo = i
      }
      hi = i
    }
  }
  (lo, hi)
}

///|
/// T71: the skein recurrence for the 2-strand torus family, an *independent*
/// computation of `V(T(2,n))` from the Jones skein relation
/// `t^{−1}V(L₊) − t V(L₋) = (t^{1/2} − t^{−1/2}) V(L₀)` with `L₀` the
/// 2-component unlink.  Returns `V(T(2,2k+1))` in `t`-powers.
pub fn jones_torus2_skein(k : Int32) -> Lz {
  let n = 2 * k + 1
  // x_j = V(T(2,j)) for j = 0 .. n, indexed in half-powers of t
  let mut x0 = hlz_zero()
  // V(2-unlink) = −(t^{1/2} + t^{−1/2})
  x0 = hlz_add(x0, hlz_mono(1, -1))
  x0 = hlz_add(x0, hlz_mono(-1, -1))
  let x1 = hlz_mono(0, 1) // unknot
  if n == 1 {
    return hlz_to_t(x1)
  }
  // t^{−1}x_{j+1} = t x_{j−1} + (t^{1/2} − t^{−1/2}) x_j
  //   ⇒ x_{j+1} = t² x_{j−1} + (t^{3/2} − t^{1/2}) x_j      (in half-powers:
  //      t² → 4, t^{3/2} → 3, t^{1/2} → 1)
  let a = hlz_mono(4, 1)
  let b = hlz_add(hlz_mono(3, 1), hlz_mono(1, -1))
  let mut prev = x0
  let mut cur = x1
  for _ in 1..<n {
    let nxt = hlz_add(hlz_mul(a, prev), hlz_mul(b, cur))
    prev = cur
    cur = nxt
  }
  hlz_to_t(cur)
}

///|
/// Convert a half-power polynomial whose exponents are all even into `t`-powers.
pub fn hlz_to_t(a : Hlz) -> Lz {
  let out = lz_zero()
  for i in 0..<LZ_LEN {
    let k = i - LZ_OFF
    if a.e[i] != 0 {
      if k % 2 == 0 {
        out.e[k / 2 + LZ_OFF] += a.e[i]
      }
    }
  }
  out
}

///|
/// The braid word `σ_1^k` in `B_2` (used for the `(2,2k+1)` family).
pub fn torus2_word(n : Int32) -> Array[(Int32, Int32)] {
  let w : Array[(Int32, Int32)] = []
  for _ in 0..<n {
    w.push((0, 1))
  }
  w
}

///|
/// T71 certificate: matching model vs skein recurrence at `T(2, 2k+1)`.
pub fn jones_torus2_defect(k : Int32) -> Int32 {
  let (v, defect) = jones_polynomial(2, torus2_word(2 * k + 1))
  // the skein recursion, driven by *positive* crossings, computes the mirror
  // knot (the trefoil is chiral — T70 — so this is a genuine statement):
  let s = lz_mirror_t(jones_torus2_skein(k))
  let d = lz_add(v, lz_neg(s))
  lz_support(d) + defect
}

///|
/// The mirror `t ↦ t^{−1}` of a `t`-polynomial: exponent `e ↦ −e`.
pub fn lz_mirror_t(a : Lz) -> Lz {
  let out = lz_zero()
  for i in 0..<LZ_LEN {
    let e = i - LZ_OFF
    if a.e[i] != 0 {
      out.e[-e + LZ_OFF] += a.e[i]
    }
  }
  out
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T69: the matching model satisfies the Temperley–Lieb relations" {
  let n = 4
  let delta = tl_delta()
  for i in 0..<(n - 2) {
    // e_i² = δ e_i
    let e = tl_gen(n, i)
    let e2 = tl_mul(n, e, e)
    assert_eq(e2.parts.length(), 1)
    assert_eq(lz_to_string(e2.coefs[0]), lz_to_string(delta))
    // e_i e_{i+1} e_i = e_i
    let prod = tl_mul(n, tl_mul(n, e, tl_gen(n, i + 1)), e)
    assert_eq(
      prod.parts.length(),
      1,
      msg="e_i e_{i+1} e_i is a single matching",
    )
    assert_true(
      lz_is_zero(lz_add(prod.coefs[0], lz_neg(lz_const(1)))),
      msg=lz_to_string(prod.coefs[0]),
    )
    // σ_i σ_i^{-1} = 1 in the algebra
    let s = tl_braid_generator(n, i, 1)
    let sinv = tl_braid_generator(n, i, -1)
    let ss = tl_mul(n, s, sinv)
    assert_eq(ss.parts.length(), 1, msg="σ_i σ_i^{-1} = 1")
    assert_true(lz_is_zero(lz_add(ss.coefs[0], lz_neg(lz_const(1)))))
  }
  // e_0 e_2 = e_2 e_0 (far apart) and neither is the identity
  let lhs = tl_mul(n, tl_gen(n, 0), tl_gen(n, 2))
  let rhs = tl_mul(n, tl_gen(n, 2), tl_gen(n, 0))
  assert_eq(lhs.parts.length(), 1)
  assert_true(tl_same_matching(lhs.parts[0], rhs.parts[0]))
  assert_true(!tl_same_matching(lhs.parts[0], tl_identity_matching(n)))
  // e_0 e_1 ≠ e_0: the matching model keeps the two apart (TL_4 is 14-dimensional)
  let e01 = tl_mul(n, tl_gen(n, 0), tl_gen(n, 1))
  assert_true(!tl_same_matching(e01.parts[0], tl_generator_matching(n, 0)))
  // … while e_0 e_1 e_0 = e_0
  let back = tl_mul(n, e01, tl_gen(n, 0))
  assert_eq(back.parts.length(), 1)
  assert_true(tl_same_matching(back.parts[0], tl_generator_matching(n, 0)))
  assert_true(lz_is_zero(lz_add(back.coefs[0], lz_neg(lz_const(1)))))
  // dimension check: the matchings are Catalan many and pairwise distinct
  assert_eq(tl_non_crossing_matchings(3).length(), 5) // Catalan(3)
  assert_eq(tl_non_crossing_matchings(n).length(), 14) // Catalan(4)
  assert_eq(tl_matchings(3).length(), 15) // all matchings: 5!!
}

///|
test "T70: Jones polynomials of the trefoil, figure-eight and 5_1" {
  // left-handed trefoil = closure of σ_1³ in B_2
  let (v3, d3) = jones_polynomial(2, torus2_word(3))
  assert_eq(d3, 0)
  let expect3 = lz_add(lz_mono(-4, -1), lz_add(lz_mono(-3, 1), lz_mono(-1, 1)))
  assert_true(lz_is_zero(lz_add(v3, lz_neg(expect3))), msg=lz_to_string(v3))
  // figure-eight = closure of (σ_1 σ_2^{-1})² in B_3
  let w41 : Array[(Int32, Int32)] = [(0, 1), (1, -1), (0, 1), (1, -1)]
  let (v4, d4) = jones_polynomial(3, w41)
  assert_eq(d4, 0)
  let expect4 = lz_add(
    lz_mono(-2, 1),
    lz_add(
      lz_mono(-1, -1),
      lz_add(lz_mono(0, 1), lz_add(lz_mono(1, -1), lz_mono(2, 1))),
    ),
  )
  assert_true(lz_is_zero(lz_add(v4, lz_neg(expect4))), msg=lz_to_string(v4))
  // the coefficient interface: V(3₁) has c_{−4} = −1, c_{−3} = +1, c_{−1} = +1
  assert_eq(lz_coeff(v3, -4), -1)
  assert_eq(lz_coeff(v3, -3), 1)
  assert_eq(lz_coeff(v3, -1), 1)
  assert_eq(lz_coeff(v3, 0), 0)
  // 2·V(3₁) doubles every coefficient and leaves the shape alone
  let two_v3 = lz_scale(v3, 2)
  assert_eq(lz_coeff(two_v3, -4), -2)
  assert_eq(lz_support(two_v3), lz_support(v3))
  // reversing the window twice is the identity (used by lz_mirror_t's tests)
  assert_true(lz_is_zero(lz_add(lz_reverse(lz_reverse(v3)), lz_neg(v3))))
  // V(1) = 1 and |V(−1)| = determinant
  assert_eq(lz_eval_one(v3), 1)
  assert_eq(lz_eval_one(v4), 1)
  assert_eq(lz_eval_minus_one(v3).abs(), 3) // det(3_1) = 3
  assert_eq(lz_eval_minus_one(v4).abs(), 5) // det(4_1) = 5
  assert_eq(lz_span(v3), (-4, -1))
}

///|
test "T70: the trefoil is chiral; 4_1 is amphichiral" {
  let (v3, _) = jones_polynomial(2, torus2_word(3))
  // reverse the support window explicitly: t^e ↦ t^{−e}
  let rev3 = lz_mirror_t(v3)
  assert_eq(lz_support(lz_add(v3, lz_neg(rev3))), 6) // genuinely different
  let w41 : Array[(Int32, Int32)] = [(0, 1), (1, -1), (0, 1), (1, -1)]
  let (v4, _) = jones_polynomial(3, w41)
  let rev4 = lz_mirror_t(v4)
  assert_true(lz_is_zero(lz_add(v4, lz_neg(rev4))))
}

///|
test "T71: matching model agrees with the skein recurrence on T(2,2k+1)" {
  for k in 1..<=6 {
    assert_eq(jones_torus2_defect(k), 0, msg="k = " + k.to_string())
  }
  // determinantal law at the torus knots: |V(−1)| = 2k+1
  for k in 1..<=6 {
    let (v, _) = jones_polynomial(2, torus2_word(2 * k + 1))
    assert_eq(lz_eval_minus_one(v).abs(), 2 * k + 1)
    assert_eq(lz_eval_one(v), 1)
  }
}
