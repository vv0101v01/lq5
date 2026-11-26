// =====================================================================
// 41 — 視窗無關的 T(2,n) Jones 計算（|V(−1)| = n 行列式定律的
// 視窗伸縮律）。第九輪 deferred 項目 2（docs/OCTAVE.md §1.4 / :67
// 「|V(−1)| 視窗定律」），對應 A71–A72，證書 T88–T91。
//
// 模組 34 的 Lz 視窗固定 ±96：n ≥ 17 時 A^{-3n} 正規化把範式推出
// 視窗（真值 spanA = (−3n+2, n+2)，真值 spanN = (−6n+2, −2n+2)），
// lz_mul 靜默丟項 —— V(1)/V(−1) 離開行等式而 lattice defect 仍為 0
// （T90 釘死的壞值）。本模組：
//   * Lz2：視窗 off 由呼叫者給出（與 Lz 同構，指數 ∈ [−off, off)）；
//   * TL_2 乘法帶丟項計數 (c0,c1)·(d0,d1) = (c0d0, c0d1+c1d0+δc1d1)，
//     c1 係數 O(n)（Int32 安全；係數絕對值 ≤ n+1）；
//   * t2_window_bound(n) = 6n−2：spanN[0] = −6n+2 = −bound 恰落在
//     下界（idx 0），spanN[1] = −2n+2 ≤ off−1（n ≥ 2），故 bound 之下
//     零丟項 —— 此即視窗伸縮律。
// 範圍：n ≥ 1 的鏈環/結；偶數 n 給連結（lattice defect > 0，屬預期）。
// =====================================================================

///|
/// 視窗化的 A-多項式：指數 e ∈ [−off, off)，長度 2·off。
pub(all) struct Lz2 {
  off : Int32
  e : Array[Int32]
}

///|
/// 零多項式。
pub fn lz2_zero(off : Int32) -> Lz2 {
  Lz2::{ off, e: Array::make(2 * off, 0), }
}

///|
/// c·A^exp；exp 必須落在視窗內（否則 panic —— 調用者受 bound 約束）。
pub fn lz2_mono(off : Int32, exp : Int32, c : Int32) -> Lz2 {
  let p = lz2_zero(off)
  p.e[exp + off] = c
  p
}

///|
/// 佔用指標範圍 [lo, hi]（零多項式返回 (off, off−1)）。
fn lz2_index_span(a : Lz2) -> (Int32, Int32) {
  let mut lo = 2 * a.off
  let mut hi = 0 - 1
  for i in 0..<(2 * a.off) {
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
/// 逐項加（兩側需同 off）。
pub fn lz2_add(a : Lz2, b : Lz2) -> Lz2 {
  let out = lz2_zero(a.off)
  for i in 0..<(2 * a.off) {
    out.e[i] = a.e[i] + b.e[i]
  }
  out
}

///|
/// 逐項乘；返回 (乘積, 丟項數) —— 模組 34 靜默丟棄的這裡計數。
pub fn lz2_mul(a : Lz2, b : Lz2) -> (Lz2, Int32) {
  let off = a.off
  let out = lz2_zero(off)
  let mut drops = 0
  let (alo, ahi) = lz2_index_span(a)
  let (blo, bhi) = lz2_index_span(b)
  if alo > ahi || blo > bhi {
    return (out, 0)
  }
  for i in alo..<=ahi {
    if a.e[i] != 0 {
      for j in blo..<=bhi {
        if b.e[j] != 0 {
          let k = i + j - off
          if k >= 0 && k < 2 * off {
            out.e[k] += a.e[i] * b.e[j]
          } else {
            drops = drops + 1
          }
        }
      }
    }
  }
  (out, drops)
}

///|
/// 支撐集的指數跨度 (lo, hi)；零多項式返回 (0, 0)。
pub fn lz2_span(a : Lz2) -> (Int32, Int32) {
  let (lo, hi) = lz2_index_span(a)
  if lo > hi {
    return (0, 0)
  }
  (lo - a.off, hi - a.off)
}

///|
/// A^exp 的係數（視窗外為 0）。
pub fn lz2_coeff(a : Lz2, exp : Int32) -> Int32 {
  let idx = exp + a.off
  if idx < 0 || idx >= 2 * a.off {
    return 0
  }
  a.e[idx]
}

///|
/// V(1)：係數總和。
pub fn lz2_eval_one(a : Lz2) -> Int32 {
  let mut s = 0
  for i in 0..<(2 * a.off) {
    s += a.e[i]
  }
  s
}

///|
/// V(−1)：Σ c_e·(−1)^e。
pub fn lz2_eval_minus_one(a : Lz2) -> Int32 {
  let mut s = 0
  for i in 0..<(2 * a.off) {
    if a.e[i] != 0 {
      if (i - a.off) % 2 == 0 {
        s += a.e[i]
      } else {
        s -= a.e[i]
      }
    }
  }
  s
}

///|
/// A-指數 → t-指數（t = A^{−4}）：e ≡ 0 (mod 4) 才是合法 t-冪，
/// 其餘計入 lattice defect（返回 (t-多項式, defect)）。
pub fn lz2_to_t(a : Lz2) -> (Lz2, Int32) {
  let off = a.off
  let out = lz2_zero(off)
  let mut defect = 0
  for i in 0..<(2 * off) {
    let c = a.e[i]
    if c != 0 {
      let e = i - off
      if e % 4 == 0 {
        out.e[e / 4 + off] += c
      } else {
        defect += 1
      }
    }
  }
  (out, defect)
}

///|
/// δ = A^2 + A^{−2}（圈值）。
fn lz2_delta(off : Int32) -> Lz2 {
  lz2_add(lz2_mono(off, 2, 0 - 1), lz2_mono(off, 0 - 2, 0 - 1))
}

///|
/// π(σ) = A·1 + A^{−1}·e 在 TL_2 中的 n 次冪：(c0, c1)。
/// 乘法規則（tl_compose 驗證）：e·e → e（1 圈），故
/// (c0,c1)·(d0,d1) = (c0d0, c0d1 + c1d0 + δ·c1d1)；c1 係數 O(n)。
fn t2_braid_power(n : Int32, off : Int32) -> (Lz2, Lz2, Int32) {
  let mut c0 = lz2_mono(off, 0, 1)
  let mut c1 = lz2_zero(off)
  let mut drops = 0
  let delta = lz2_delta(off)
  let d0 = lz2_mono(off, 1, 1)
  let d1 = lz2_mono(off, 0 - 1, 1)
  for _ in 0..<n {
    let (p0, dr0) = lz2_mul(c0, d0)
    let (q1, dr1) = lz2_mul(c0, d1)
    let (r1, dr2) = lz2_mul(c1, d0)
    let (s1a, dr3) = lz2_mul(c1, d1)
    let (s1, dr4) = lz2_mul(delta, s1a)
    drops += dr0 + dr1 + dr2 + dr3 + dr4
    c0 = p0
    c1 = lz2_add(lz2_add(q1, r1), s1)
  }
  (c0, c1, drops)
}

///|
/// ⟨T(2,n)⟩ = tr(π(σ)^n) = δ·c0 + c1（A-多項式），連同丟項數。
pub fn t2_bracket(n : Int32, off : Int32) -> (Lz2, Int32) {
  let (c0, c1, drops) = t2_braid_power(n, off)
  let (dc0, dr) = lz2_mul(lz2_delta(off), c0)
  (lz2_add(dc0, c1), drops + dr)
}

///|
/// 範式 (−A^3)^{−n}·⟨T(2,n)⟩ = (−1)^n A^{−3n}·⟨…⟩，連同丟項數。
pub fn t2_norm(n : Int32, off : Int32) -> (Lz2, Int32) {
  let (br, drops) = t2_bracket(n, off)
  let sign = if n % 2 == 0 { 1 } else { 0 - 1 }
  let (nm, dr) = lz2_mul(br, lz2_mono(off, 0 - 3 * n, sign))
  (nm, drops + dr)
}

///|
/// A72: 視窗下界 —— spanN = (−6n+2, −2n+2) 完整落在 [−off, off)
/// 所需的最小 off（n ≥ 1：6n−2 ≥ 4 ≥ 2 且 −3n ≥ −(6n−2)）。
pub fn t2_window_bound(n : Int32) -> Int32 {
  6 * n - 2
}

///|
/// A71: 視窗無關的 Jones 計算：返回 (t-多項式, 丟項數, lattice defect)。
/// 在 off = t2_window_bound(n) 下零丟項，V 滿足 |V(−1)| = n（奇數 n）。
pub fn jones_t2_scaled(n : Int32, off : Int32) -> (Lz2, Int32, Int32) {
  let (nm, drops) = t2_norm(n, off)
  let (tp, defect) = lz2_to_t(nm)
  (tp, drops, defect)
}

// ---------------------------------------------------------------------
// T88: 視窗伸縮律 —— bound 之下零丟項，且 span 準確命中解析邊界
// ---------------------------------------------------------------------

///|
test "T88: window-bound span law with zero drops" {
  // n = 1 是唯一例外（δ·c0 的 −1 冪與 c1 抵消）：
  let (br1, dr1) = t2_bracket(1, t2_window_bound(1))
  let (lo1, hi1) = lz2_span(br1)
  assert_true(dr1 == 0 && lo1 == 3 && hi1 == 3, msg="n=1 bracket span (3,3)")
  let (nm1, dn1) = t2_norm(1, t2_window_bound(1))
  let (nlo1, nhi1) = lz2_span(nm1)
  assert_true(dn1 == 0 && nlo1 == 0 && nhi1 == 0, msg="n=1 norm span (0,0)")
  // n ≥ 2：spanA = (−3n+2, n+2)，spanN = (−6n+2, −2n+2)，drops = 0
  let mut checked = 0
  for n in 2..<=41 {
    let off = t2_window_bound(n)
    let (br, dr) = t2_bracket(n, off)
    let (lo, hi) = lz2_span(br)
    assert_true(dr == 0, msg="bracket: zero drops at bound")
    assert_true(lo == 0 - 3 * n + 2 && hi == n + 2, msg="spanA = (-3n+2, n+2)")
    let (nm, dn) = t2_norm(n, off)
    let (nlo, nhi) = lz2_span(nm)
    assert_true(dn == 0, msg="norm: zero drops at bound")
    assert_true(
      nlo == 0 - 6 * n + 2 && nhi == 0 - 2 * n + 2,
      msg="spanN = (-6n+2, -2n+2)",
    )
    checked = checked + 1
  }
  assert_true(checked == 40, msg="n = 2..41 all checked")
}

// ---------------------------------------------------------------------
// T89: 行等式定律（奇數 n ≤ 61）+ bound−1 恰差一項即破
// ---------------------------------------------------------------------

///|
test "T89: determinant law |V(-1)| = n and the tight threshold" {
  // (a) 奇數 n：V(1) = 1，|V(−1)| = n，lattice defect = 0，drops = 0
  let mut tested = 0
  for i in 0..<30 {
    let n = 2 * i + 3
    let off = t2_window_bound(n)
    let (v, drops, defect) = jones_t2_scaled(n, off)
    assert_true(
      drops == 0 && defect == 0,
      msg="scaled: zero drops and zero defect",
    )
    assert_true(lz2_eval_one(v) == 1, msg="V(1) = 1")
    assert_true(lz2_eval_minus_one(v).abs() == n, msg="|V(-1)| = n")
    tested = tested + 1
  }
  assert_true(tested == 30, msg="odd n = 3..61 all checked")
  // (b) off = bound−1：恰在下界丟一項 ⇒ V(1) = 2、|V(−1)| = n−1
  for n in [17, 19, 21] {
    let off = t2_window_bound(n) - 1
    let (v, drops, defect) = jones_t2_scaled(n, off)
    assert_true(drops == 1, msg="exactly one dropped term below the bound")
    assert_true(defect == 0, msg="defect still invisible after truncation")
    assert_true(lz2_eval_one(v) == 2, msg="V(1) jumps 1 -> 2")
    let vm1 = lz2_eval_minus_one(v)
    assert_true(vm1 == n - 1 || vm1 == 0 - (n - 1), msg="|V(-1)| = n-1")
  }
}

// ---------------------------------------------------------------------
// T90: 模組 34 固定 ±96 視窗的靜默壞值（lattice defect 看不見）
// ---------------------------------------------------------------------

///|
test "T90: fixed-window truncation breaks the law silently" {
  let ns = [17, 19, 21, 25, 31]
  let want1 = [2, 1, 2, 2, 1]
  let wantm = [16, 0 - 15, 14, 12, 0 - 9]
  let mut i = 0
  while i < 5 {
    let n = ns[i]
    let (v34, defect) = jones_polynomial(2, torus2_word(n))
    // 釘死實測壞值：V(1)/V(−1) 離開行等式，而 lattice defect 仍為 0 ——
    // 窗口截斷對 lattice 檢測完全隱形，只有行等式本身能發現。
    assert_true(defect == 0, msg="defect stays 0: truncation is invisible")
    assert_true(lz_eval_one(v34) == want1[i], msg="fixed-window V(1) pinned")
    assert_true(
      lz_eval_minus_one(v34) == wantm[i],
      msg="fixed-window V(-1) pinned",
    )
    i = i + 1
  }
}

// ---------------------------------------------------------------------
// T91: 三路交叉驗證 —— n ≤ 15 與模組 34 一致；n ≥ 17 與 skein 一致
// ---------------------------------------------------------------------

///|
test "T91: cross-algorithm agreement scaled vs module34 vs skein" {
  // (a) 視窗尚未綁緊（n ≤ 15）：scaled 與模組 34 逐項一致
  for n in [3, 5, 7, 9, 11, 13, 15] {
    let (vref, _) = jones_polynomial(2, torus2_word(n))
    let (v, drops, defect) = jones_t2_scaled(n, t2_window_bound(n))
    assert_true(drops == 0 && defect == 0, msg="scaled clean at n<=15")
    let (rlo, rhi) = lz_span(vref)
    let (vlo, vhi) = lz2_span(v)
    assert_true(rlo == vlo && rhi == vhi, msg="span equal to module 34")
    let mut e = rlo
    while e <= rhi {
      assert_true(
        lz_coeff(vref, e) == lz2_coeff(v, e),
        msg="coeff equal to module 34",
      )
      e = e + 1
    }
    assert_true(lz_eval_one(vref) == lz2_eval_one(v), msg="V(1) equal")
    assert_true(
      lz_eval_minus_one(vref) == lz2_eval_minus_one(v),
      msg="V(-1) equal",
    )
  }
  // (b) 越過破壞點（k = 8,9,10,12,15 → n = 17,19,21,25,31）：
  // skein 走 t-半冪視窗（±96 ⇒ t ∈ [−48,47]，n ≤ 31 的真值 t-跨度
  // (−46,−15) 以內）不截斷，故與 scaled 一致 —— 模組 34 的
  // jones_polynomial（A-視窗）在此已壞（T90），三路不再同真。
  for k in [8, 9, 10, 12, 15] {
    let n = 2 * k + 1
    let sk = lz_mirror_t(jones_torus2_skein(k))
    let (v, drops, defect) = jones_t2_scaled(n, t2_window_bound(n))
    assert_true(drops == 0 && defect == 0, msg="scaled clean beyond breakdown")
    let (slo, shi) = lz_span(sk)
    let (vlo, vhi) = lz2_span(v)
    assert_true(slo == vlo && shi == vhi, msg="span equal to skein")
    let mut e = slo
    while e <= shi {
      assert_true(
        lz_coeff(sk, e) == lz2_coeff(v, e),
        msg="coeff equal to skein",
      )
      e = e + 1
    }
    assert_true(
      lz_eval_one(sk) == 1 && lz_eval_minus_one(sk).abs() == n,
      msg="skein itself satisfies the law",
    )
  }
}
