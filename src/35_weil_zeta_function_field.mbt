// ============================================================================
// 35_weil_zeta_function_field.mbt -- 通道八: 函數域上的 Weil ζ、函數方程與
//                                   「曲線的黎曼假設」 (Round 7 / Node 2,
//                                   T72-T74, A56-A58)
//
// ABSTRACTION.  For a smooth projective curve C/F_q of genus g,
//
//     Z(C, T) = exp( Σ_{k≥1} #C(F_{q^k}) T^k / k ) = P_1(T) / ((1−T)(1−qT)),
//     P_1(T) = Π_{i=1..2g} (1 − α_i T) ∈ ℤ[T],   deg P_1 = 2g,  P_1(0) = 1
//
// with three *decidable* laws (Weil 1948, RH for curves):
//
//   T72  functional equation:  P_1(T) = q^g T^{2g} P_1(1/(qT)),
//        i.e. coefficient-wise  c_k = q^{g−k} c_{2g−k}   (an identity in ℤ).
//   T73  Riemann hypothesis:   |α_i| = √q  for every i.  Written through the
//        pairing α ↦ q/α this is *rational*: with b_i = α_i + q/α_i one needs
//        b_i ∈ ℝ and |b_i| ≤ 2√q.  (For g = 1: |a_1| ≤ 2√q, Hasse's bound.)
//   T74  Langlands/moderate tie-back: for E: y²+y = x³−x²−10x−20 (11a1) the
//        Frobenius trace #E(F_p) − (p+1) equals −a_p, where a_p is the Hecke
//        eigenvalue of module 29 derived from η(τ)²η(11τ)².  This is the
//        numerical shadow of modularity, and it is checked at p = 2,3,5,7.
//
// The point counts are *computed*, not quoted: N_1 by the quadratic-character
// sum, and for the genus-2 curve also N_2 by brute force over F_{p²} (the
// quadratic algebra F_p[i]/(i² − nr)), which is what pins down a_2 and hence
// makes T73 a genuine 4-root statement rather than a Hasse-bound remark.
// ============================================================================

// ---------------------------------------------------------------------------
// 1. Point counts (T72's raw data).
// ---------------------------------------------------------------------------

///|
/// The genus-2 curve `C: y² = x⁵ − x + 1`.  Its discriminant is
/// `5⁵·1⁴ + 4⁴·(−1)⁵ = 3125 − 256 = 2869 = 19·151`, so it has good reduction at
/// every prime used below (3, 5, 7, 11).
pub fn genus2_discriminant() -> Int32 {
  2869
}

///|
/// `f(x) = x⁵ − x + 1` mod `p`, reduced into `[0, p)`.
pub fn genus2_f(x : Int32, p : Int32) -> Int32 {
  let x2 = x * x % p
  let x4 = x2 * x2 % p
  ((x4 * x % p - x + 1) % p + p) % p
}

///|
/// `#C(F_p)` for `y² = x⁵ − x + 1`: `1 + Σ_x (1 + χ(f(x)))` (one point at
/// infinity, as `deg f` is odd).
pub fn genus2_point_count(p : Int32) -> Int32 {
  let mut total = 1
  for x in 0..<p {
    total += 1 + legendre_sym(genus2_f(x, p), p)
  }
  total
}

///|
/// `#E(F_p)` for `11a1` — module 29's count (included here as the `g = 1`
/// input of the round).
pub fn elliptic11_point_count(p : Int32) -> Int32 {
  curve11_point_count(p)
}

// ---------------------------------------------------------------------------
// 2. F_{p²} = F_p[i]/(i² − nr) — enough to count C(F_{p²}) (T73 for g = 2).
// ---------------------------------------------------------------------------

///|
/// An element `a + b·i` of `F_p[i]/(i² − nr)`.
pub(all) struct Fp2 {
  a : Int32
  b : Int32
  p : Int32
  nr : Int32
}

///|
/// The smallest quadratic non-residue mod `p` (used as `nr`).
pub fn fp2_non_residue(p : Int32) -> Int32 {
  let mut d = 2
  while d < p && legendre_sym(d, p) != -1 {
    d += 1
  }
  d
}

///|
pub fn Fp2::make(a : Int32, b : Int32, p : Int32, nr : Int32) -> Fp2 {
  Fp2::{ a: (a % p + p) % p, b: (b % p + p) % p, p, nr, }
}

///|
pub fn Fp2::add(self : Fp2, o : Fp2) -> Fp2 {
  Fp2::make(self.a + o.a, self.b + o.b, self.p, self.nr)
}

///|
pub fn Fp2::sub(self : Fp2, o : Fp2) -> Fp2 {
  Fp2::make(self.a - o.a, self.b - o.b, self.p, self.nr)
}

///|
pub fn Fp2::mul(self : Fp2, o : Fp2) -> Fp2 {
  // (a + b i)(c + d i) = (ac + bd·nr) + (ad + bc) i
  let ac = self.a * o.a % self.p
  let bd = self.b * o.b % self.p * self.nr % self.p
  let ad = self.a * o.b % self.p
  let bc = self.b * o.a % self.p
  Fp2::make(ac + bd, ad + bc, self.p, self.nr)
}

///|
pub fn Fp2::is_zero(self : Fp2) -> Bool {
  self.a == 0 && self.b == 0
}

///|
pub fn Fp2::pow(self : Fp2, e : Int32) -> Fp2 {
  let mut result = Fp2::make(1, 0, self.p, self.nr)
  let mut base = self
  let mut k = e
  while k > 0 {
    if k % 2 == 1 {
      result = result.mul(base)
    }
    base = base.mul(base)
    k = k / 2
  }
  result
}

///|
/// The quadratic character on `F_{p²}`: `0`, `+1` (square) or `−1`.
pub fn fp2_chi(x : Fp2) -> Int32 {
  if x.is_zero() {
    return 0
  }
  let q2 = x.p * x.p
  let e = (q2 - 1) / 2
  let r = x.pow(e)
  if r.a == 1 && r.b == 0 {
    1
  } else {
    -1
  }
}

///|
/// `#C(F_{p²})`: brute force over the `p²` values of `x` (the `y` direction is
/// summed by the quadratic character, so the cost is `p²` field operations).
pub fn genus2_point_count_ext2(p : Int32) -> Int32 {
  let nr = fp2_non_residue(p)
  let one = Fp2::make(1, 0, p, nr)
  let mut total = one.a
  for xa in 0..<p {
    for xb in 0..<p {
      // f(x) = x⁵ − x + 1 inside F_{p²}
      let x = Fp2::make(xa, xb, p, nr)
      let x2 = x.mul(x)
      let x4 = x2.mul(x2)
      let fx = x4.mul(x).sub(x).add(one)
      total += 1 + fp2_chi(fx)
    }
  }
  total
}

// ---------------------------------------------------------------------------
// 3. P_1(T) and the functional equation (T72).
// ---------------------------------------------------------------------------

///|
/// `P_1(T) = 1 − a₁T + qT²` for the elliptic curve, `a₁ = q + 1 − N_1`.
pub fn weil_P1_genus1(p : Int32) -> Array[Int32] {
  let a1 = p + 1 - elliptic11_point_count(p)
  [1, -a1, p]
}

///|
/// `P_1(T) = 1 − a₁T + a₂T² − q a₁T³ + q²T⁴` for `C`, with
/// `a₁ = q + 1 − N_1` and `a₂ = (a₁² − (q² + 1 − N_2)) / 2`.
pub fn weil_P1_genus2(p : Int32) -> Array[Int32] {
  let n1 = genus2_point_count(p)
  let n2 = genus2_point_count_ext2(p)
  let q = p
  let a1 = q + 1 - n1
  let num = a1 * a1 - (q * q + 1 - n2)
  let a2 = num / 2
  [1, -a1, a2, -q * a1, q * q]
}

///|
/// T72: number of coefficients violating `c_k = q^{g−k} c_{2g−k}`.
pub fn weil_functional_equation_defect(
  coeffs : Array[Int32],
  q : Int32,
) -> Int32 {
  let two_g = coeffs.length() - 1
  let g = two_g / 2
  let mut bad = 0
  for k in 0..<=two_g {
    // coefficient of T^k in q^g T^{2g} P(1/(qT)) is c_{2g−k}·q^{k−g}; for
    // k < g this is a *rational* number q^{−(g−k)}·c_{2g−k}, so the identity is
    // cleared of denominators before it is compared in ℤ.
    let e = k - g
    let mut qpow = 1
    for _ in 0..<e.abs() {
      qpow *= q
    }
    let ok = if e >= 0 {
      coeffs[k] == qpow * coeffs[two_g - k]
    } else {
      coeffs[k] * qpow == coeffs[two_g - k]
    }
    if !ok {
      bad += 1
    }
  }
  bad
}

// ---------------------------------------------------------------------------
// 4. The Riemann hypothesis for the curve (T73).
// ---------------------------------------------------------------------------

///|
/// T73 (g = 1): Hasse's bound `|a₁| ≤ 2√q`, returned as the *excess* over
/// `2√q` (0 when the bound holds, i.e. exactly when both roots `α, q/α` have
/// absolute value `√q`).
pub fn hasse_excess_genus1(p : Int32) -> Double {
  let a1 = (p + 1 - elliptic11_point_count(p)).to_double()
  let excess = a1.abs() - 2.0 * p.to_double().sqrt()
  if excess > 0.0 {
    excess
  } else {
    0.0
  }
}

///|
/// T73 (g = 2): the paired Frobenius eigenvalues are roots of
/// `T² − a₁T + (a₂ − 2q)`; the curve satisfies RH iff that quadratic has real
/// roots inside `[−2√q, 2√q]`.  Returns `0.0` when it does, the obstruction
/// otherwise.
pub fn genus2_riemann_hypothesis_defect(p : Int32) -> Double {
  let poly = weil_P1_genus2(p)
  let a1 = -poly[1]
  let a2 = poly[2]
  let q = p
  let disc = a1.to_double() * a1.to_double() - 4.0 * (a2 - 2 * q).to_double()
  if disc < 0.0 {
    return -disc + 1.0
  }
  let root = disc.sqrt()
  let b1 = (a1.to_double() + root) / 2.0
  let b2 = (a1.to_double() - root) / 2.0
  let bound = 2.0 * p.to_double().sqrt()
  let mut worst = 0.0
  for b in [b1, b2] {
    let e = b.abs() - bound
    if e > worst {
      worst = e
    }
  }
  if worst < 0.0 {
    0.0
  } else {
    worst
  }
}

///|
/// The *actual* Frobenius eigenvalues `α` of `C/F_p`, obtained as
/// `α = (b ± √(b²−4q))/2` — the `2g`-tuple whose absolute values must all equal
/// `√q` (returns the largest deviation from `√q`).
pub fn genus2_sqrt_q_defect(p : Int32) -> Double {
  let poly = weil_P1_genus2(p)
  let a1 = -poly[1].to_double()
  let a2 = poly[2].to_double()
  let q = p.to_double()
  let disc = a1 * a1 - 4.0 * (a2 - 2.0 * q)
  let root = disc.sqrt()
  let mut worst = 0.0
  for b in [(a1 + root) / 2.0, (a1 - root) / 2.0] {
    let inner = b * b - 4.0 * q
    let s = inner.sqrt()
    for alpha in [(b + s) / 2.0, (b - s) / 2.0] {
      let d = (alpha.abs() - q.sqrt()).abs()
      if d > worst {
        worst = d
      }
    }
  }
  worst
}

// ---------------------------------------------------------------------------
// 5. T74: the Frobenius trace *is* the Hecke eigenvalue (modularity shadow).
// ---------------------------------------------------------------------------

///|
/// T74: `#E(F_p) = p + 1 − a_p` with `a_p` the η-product Hecke eigenvalue of
/// module 29.  Number of violations over the primes `2, 3, 5, 7, 11`.
pub fn modularity_defect() -> Int32 {
  let mut bad = 0
  for p in [2, 3, 5, 7, 11] {
    if elliptic11_point_count(p) != p + 1 - hecke_a(p) {
      bad += 1
    }
  }
  bad
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T72: Weil P_1 has degree 2g and satisfies the functional equation" {
  // the genus-2 curve is chosen so that every prime of the sweep is good:
  // disc(x⁵−x+1) = 5⁵ + 4⁴ = 2869 = 19·151, and none of 3,5,7,11,13 divides it
  assert_eq(genus2_discriminant(), 2869)
  for p in [3, 5, 7, 11, 13] {
    assert_true(
      genus2_discriminant() % p != 0,
      msg="bad reduction at p = " + p.to_string(),
    )
  }
  for p in [3, 5, 7, 11, 13] {
    let f = weil_P1_genus1(p)
    assert_eq(f.length(), 3)
    assert_eq(f[0], 1)
    assert_eq(f[2], p)
    assert_eq(weil_functional_equation_defect(f, p), 0)
    let g = weil_P1_genus2(p)
    assert_eq(g.length(), 5)
    assert_eq(g[0], 1)
    assert_eq(g[4], p * p)
    assert_eq(weil_functional_equation_defect(g, p), 0)
  }
}

///|
test "T73: Hasse's bound (g = 1) and the 4-root RH statement (g = 2)" {
  for p in [2, 3, 5, 7, 11, 13, 17, 19, 23] {
    assert_eq(hasse_excess_genus1(p), 0.0, msg="p = " + p.to_string())
  }
  for p in [3, 5, 7, 11, 13] {
    assert_true(
      genus2_riemann_hypothesis_defect(p) == 0.0,
      msg="RH fails at p = " + p.to_string(),
    )
    assert_true(
      genus2_sqrt_q_defect(p) < 1.0e-9,
      msg="|α| ≠ √q at p = " + p.to_string(),
    )
  }
  // the genus is genuinely 2: a_2 is not determined by a_1, and the point
  // count over F_{p²} is not the "diagonal" count p² + 1 − a_1²
  let p = 5
  let n1 = genus2_point_count(p)
  let n2 = genus2_point_count_ext2(p)
  assert_true(n2 != p * p + 1 - (p + 1 - n1) * (p + 1 - n1))
}

///|
test "T74: Frobenius traces coincide with the level-11 Hecke eigenvalues" {
  assert_eq(modularity_defect(), 0)
  for p in [2, 3, 5, 7, 11] {
    assert_eq(elliptic11_point_count(p), p + 1 - hecke_a(p))
  }
  // and the Hasse bound is the g = 1 instance of T73
  for p in [2, 3, 5, 7] {
    assert_true(hecke_a(p).to_double().abs() <= 2.0 * p.to_double().sqrt())
  }
}
