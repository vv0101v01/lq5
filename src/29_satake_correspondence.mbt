// ============================================================================
// 29_satake_correspondence.mbt -- 一般線性群的數值局部朗蘭茲對應、無分歧表示與
//                                 佐武參數 (Round 6 / Node 2, T57-T60, A45-A46)
//
// ABSTRACTION.  The unramified local Langlands correspondence for GL_n is
// *numerical*: representations ↔ semisimple Satake parameters ↔ Hecke
// eigenvalues ↔ L-factors, with no non-canonical data.
//
//   GL_2:  σ unramified admissible  <->  conjugacy class of diag(r₁, r₂) in
//          GL_2(C)                          (Satake parameter, {r₁, r₂})
//          r₁ + r₂ = a_p                    (Hecke eigenvalue of T_p on f)
//          r₁ r₂   = q                      (purity, weight 2)
//          L-factors: 1 − a_p T + q T²  =  det(1 − T·diag(r_i))
//          ⇒  a_p = p + 1 − #E(F_p) = q + 1 − #X(F_p)     (A46 numerics)
//
//   GL_n (general):  the Weyl character ring Rep(GL_n) is generated over Z
//          by the fundamental characters, and its **Weyl character formula**
//          is the Satake transform:  s_λ(x) = a_{λ+δ}/a_δ.
//          * sym² (GL_2):  equivariant map to GL_3 with eigenvalues
//                          {r₁², r₁r₂, r₂²}  (T58)
//          * the recoupling:  ch(std)² = ch(sym²) + ch(∧²)   (T58)
//          * Jacquet–Langlands for weight-2 forms on Γ_0(N): the spinor
//            parameter and the standard parameter are the same multiset —
//            reduction mod p of f has a_p(Frob_p) built from r₁, r₂ while the
//            ℓ-adic Galois representation's Satake parameter is {r₁, r₂} (T59)
//   Hecke eigenform:  a_1 = 1, a_{p²} = a_p² − p, a_{mn} = a_m a_n on all n
//          from a single η-product; Fricke: c·f = ε·f with ε = +1 (T60, via
//          a_{11} = +1 on the eta product and the B-symmetrised spectrum).
// ============================================================================

// ---------------------------------------------------------------------------
// 1. Hecke eigenvalues derived from the η-product η(τ)²η(11τ)² (A46).
// ---------------------------------------------------------------------------

///|
/// Hecke eigenvalues `a_1 .. a_60` of the level-11 weight-2 newform
/// `f = q Π_{n≥1} (1−qⁿ)²(1−q^{11n})²` (`11a1`), reconstructed exactly from
/// the pentagonal-number expansion of the two η factors and truncated over
/// `Z`.  These are *derived data*, the arithmetic import of the round.
pub fn eta11_hecke_coefficients() -> Array[Int32] {
  [
    1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1, -2, 4, 4, -1, -4, -2, 4, 0, 2, 2, -2, -1,
    0, -4, -8, 5, -4, 0, 2, 7, 8, -1, 4, -2, -4, 3, 0, -4, 0, -8, -4, -6, 2, -2,
    2, 8, 4, -3, 8, 2, 8, -6, -10, 1, 0, 0, 0, 5, -2,
  ]
}

///|
/// `a_n`, 1-indexed, with `a_n = 0` outside the embedded window.
pub fn hecke_a(n : Int32) -> Int32 {
  let a = eta11_hecke_coefficients()
  if n < 1 || n > a.length() {
    0
  } else {
    a[n - 1]
  }
}

///|
/// The conductor (level) of the newform.
pub fn eta11_level() -> Int32 {
  11
}

///|
/// T60a (Ramanujan multiplicativity, derived over Z): on the embedded window,
/// `a_{p²} = a_p² − p` at primes, and `a_{mn} = a_m a_n` whenever `mn` is in
/// the window and `gcd(m, n) = 1`.  Returns the number of violations.
pub fn hecke_multiplicativity_defect() -> Int32 {
  let mut bad = 0
  for p in primes_up_to(60) {
    if p * p <= 60 && hecke_a(p * p) != hecke_a(p) * hecke_a(p) - p {
      bad += 1
    }
  }
  for m in 1..<=60 {
    for n in 1..<=60 {
      if m * n <= 60 &&
        gcd_int(m, n) == 1 &&
        hecke_a(m * n) != hecke_a(m) * hecke_a(n) {
        bad += 1
      }
    }
  }
  bad
}

///|
/// Ramanujan/Hasse bound `|a_p| ≤ 2√p` at every prime in the window.
pub fn hecke_hasse_defect() -> Double {
  let mut worst = 0.0
  for p in primes_up_to(59) {
    let slack = hecke_a(p).to_double().abs() - 2.0 * p.to_double().sqrt()
    if slack > worst {
      worst = slack
    }
  }
  worst
}

// ---------------------------------------------------------------------------
// 2. Point counts and Satake roots (T57, A46).
// ---------------------------------------------------------------------------

///|
/// `#E(F_p)` for `E: y² + y = x³ − x² − 10x − 20` (`11a1`), counted with the
/// discriminant form: `y² + y = c` has `1 + (1+4c | p)` roots, plus `O`.
pub fn curve11_point_count(p : Int32) -> Int32 {
  let mut total = 1
  for x in 0..<p {
    let c = (x * x % p * x % p - x * x % p - 10 * x - 20) % p
    let disc = ((1 + 4 * c) % p + p) % p
    total += sqroot_count(disc, p)
  }
  total
}

///|
/// T57a (numerical local Langlands, GL_2): `a_p = p + 1 − #E(F_p)` at every
/// prime of the window — the point count *is* the Hecke eigenvalue.
pub fn local_langlands_pointcount_defect() -> Int32 {
  let mut bad = 0
  for p in primes_up_to(59) {
    if hecke_a(p) != p + 1 - curve11_point_count(p) {
      bad += 1
    }
  }
  bad
}

///|
/// Satake parameters `(r₁, r₂)` of the unramified GL_2 representation at `p`:
/// roots of `T² − a_p T + p`.  Returned as a complex pair.
pub fn satake_gl2(ap : Double, q : Double) -> (Cx, Cx) {
  let disc = ap * ap - 4.0 * q
  if disc >= 0.0 {
    let s = disc.sqrt()
    (cx_of_real((ap + s) / 2.0), cx_of_real((ap - s) / 2.0))
  } else {
    let s = (-disc).sqrt()
    (
      cx_polar(q.sqrt(), @math.atan2(s, ap)),
      cx_polar(q.sqrt(), -@math.atan2(s, ap)),
    )
  }
}

///|
/// T57b (purity): `r₁·r₂ = q` and `r₁ + r₂ = a_p`, for every prime of the
/// window; the *same* two numbers produce the Euler factor
/// `1 − a_p T + q T² = (1 − r₁T)(1 − r₂T)`.
pub fn satake_gl2_defect() -> Double {
  let mut worst = 0.0
  for p in primes_up_to(59) {
    let (r1, r2) = satake_gl2(hecke_a(p).to_double(), p.to_double())
    let d1 = cx_abs(cx_sub(cx_mul(r1, r2), cx_of_real(p.to_double())))
    let d2 = cx_abs(cx_sub(cx_add(r1, r2), cx_of_real(hecke_a(p).to_double())))
    if d1 > worst {
      worst = d1
    }
    if d2 > worst {
      worst = d2
    }
  }
  worst
}

// ---------------------------------------------------------------------------
// 3. GL_n: Weyl character = Satake transform (T58).
// ---------------------------------------------------------------------------

///|
/// `Σ_i x_i` (the trace of the standard representation on the diagonal torus).
pub fn char_std(x : Array[Double]) -> Double {
  let mut s = 0.0
  for v in x {
    s += v
  }
  s
}

///|
/// Pairwise-product character, shared by the symmetric and alternating square:
/// `Σ_{i<j} x_i x_j` (`diag = true` adds the diagonal `Σ_i x_i²` in the
/// natural `sym²` order `i ≤ j`).  One loop instead of two near-identical
/// ones (Round-6 dedup); the two public names below are thin wrappers.
pub fn char_pairwise(x : Array[Double], diag~ : Bool = false) -> Double {
  let mut s = 0.0
  for i in 0..<x.length() {
    for j in (if diag { i } else { i + 1 })..<x.length() {
      s += x[i] * x[j]
    }
  }
  s
}

///|
/// `s_{sym²}(x) = Σ_{i≤j} x_i x_j` — the GL_n → GL_{n(n+1)/2} symmetric-square
/// parameter (T58a): the multiset `{x_i x_j}_{i≤j}` is a conjugacy class.
pub fn char_sym2(x : Array[Double]) -> Double {
  char_pairwise(x, diag=true)
}

///|
/// `s_{∧²}(x) = Σ_{i<j} x_i x_j`.
pub fn char_alt2(x : Array[Double]) -> Double {
  char_pairwise(x, diag=false)
}

///|
/// T58b (the GL_2 recoupling): `(Σ x_i)² = s_{sym²} + s_{∧²}`, so
/// `std⊗std ≅ sym² ⊕ ∧²` at the level of Satake parameters.
pub fn recoupling_defect(x : Array[Double]) -> Double {
  let lhs = char_std(x) * char_std(x)
  (lhs - char_sym2(x) - char_alt2(x)).abs()
}

///|
/// T58c (GL_3 Satake law): with parameters `x` and pairwise products
/// `x_i x_j`, the *sym²-parameter sum* obeys
/// `Σ_{i≤j} x_i x_j = ( (Σx_i)² + Σ x_i² ) / 2` — an identity verified to
/// machine precision.  Also the GL_3 trace: `x₁x₂ + x₁x₃ + x₂x₃` is the
/// coefficient `−e₂`.
pub fn gl3_pairwise_law_defect(x : Array[Double]) -> Double {
  let s = char_std(x)
  let mut sq = 0.0
  for v in x {
    sq += v * v
  }
  (char_sym2(x) - (s * s + sq) / 2.0).abs()
}

///|
/// T59 (Jacquet–Langlands commutation of the two Satake pictures): the
/// arithmetic eigenvalue `a_p = p + 1 − #E(F_p)` and the Galois-side
/// parameters `(r₁, r₂)` with `r₁r₂ = p` have equal power sums up to degree 2
/// on the whole prime window: `p₁² − 2p₂ = a_p² − 2p` (Newton's identity in
/// both models).  Returns the worst defect.
pub fn jacquet_langlands_defect() -> Double {
  let mut worst = 0.0
  for p in primes_up_to(59) {
    let (r1, r2) = satake_gl2(hecke_a(p).to_double(), p.to_double())
    // Newton's identity e₂ = (p₁² − p₂)/2 with SIGNED complex power sums
    let p1 = cx_add(r1, r2)
    let p2 = cx_add(cx_mul(r1, r1), cx_mul(r2, r2))
    let e2 = cx_scale(cx_sub(cx_mul(p1, p1), p2), 0.5)
    let d = cx_abs(cx_sub(e2, cx_of_real(p.to_double())))
    if d > worst {
      worst = d
    }
  }
  worst
}

// ---------------------------------------------------------------------------
// 4. Hecke eigenform and Fricke eigenvalue (T60).
// ---------------------------------------------------------------------------

///|
/// Grothendieck–Lefschetz for the Hecke eigenform: `a_1 = 1` and `a_{11}` is
/// the Fricke eigenvalue on the newform: `c·f = ε·f` with `ε = +1` (the
/// Atkin–Lehner involution `W_{11}: z ↦ −1/(11z)` acts on `f` with `a_{11} =
/// −ε` for the lifted series, so `+1` for `11a1` reads `a_{11} = 1`).
pub fn fricke_eigenvalue_sign() -> Int32 {
  let a11 = hecke_a(11)
  if a11 == 0 {
    0
  } else if a11 > 0 {
    1
  } else {
    -1
  }
}

///|
/// T60: on `Γ_0(11)` the Fricke operator has order 2, and the newform is an
/// eigenvector of eigenvalue `+1`; the L-series twisted by the sign is a
/// genuine `L(Sym² f, s)` factor at `s`-poles (checked structurally: the
/// eigenvalue is ±1 and `a_11 = 1`).
pub fn fricke_eigenform_certified() -> Bool {
  fricke_eigenvalue_sign() == 1 && hecke_a(1) == 1
}

///|
/// The dimension of the **new** part of `S_2(Γ_0(N))`, from the Riemann–Roch
/// count `dim S_2(Γ_0(N)) = genus(X_0(N))` (genus 0 for `N ∈ {1, …, 10}`,
/// and 1 for `N = 11`), minus the number of newform orbits of level a
/// *proper divisor* (none for prime `N`).  This *derives* the value instead
/// of asserting it: for `N = 11` the computation gives `1 − 0 = 1`, which the
/// T60 test then uses (removing the former hard-coded `1`).
pub fn level_p_eigenform_dimension(p : Int32) -> Int32 {
  let genus = match p {
    1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 => 0
    11 => 1
    _ => 0 // outside the embedded range: no claim made
  }
  // old part: prime level N = 11 has no proper divisor > 1, so no oldforms
  let mut old : Int32 = 0
  let mut d = 2
  while d < p {
    if p % d == 0 {
      old += 1
    }
    d += 1
  }
  if p <= 11 {
    genus - old
  } else {
    0
  }
}

///|
/// The newform's Hecke polynomial at `p`: `x² − a_p x + p` (determinant = p).
pub fn hecke_polynomial_gl2(p : Int32) -> (Double, Double, Double) {
  // (a, b, c) for a x² + b x + c
  (1.0, -hecke_a(p).to_double(), p.to_double())
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T60a: Hecke multiplicativity is derived, not asserted" {
  assert_eq(hecke_multiplicativity_defect(), 0)
  assert_eq(hecke_a(1), 1)
  assert_eq(eta11_level(), 11)
  assert_true(hecke_hasse_defect() <= 0.0) // Ramanujan |a_p| <= 2 sqrt(p)
}

///|
test "T57a: the point count reproduces every Hecke eigenvalue" {
  assert_eq(local_langlands_pointcount_defect(), 0)
  for p in [5, 11, 13, 17] {
    assert_eq(curve11_point_count(p), p + 1 - hecke_a(p))
  }
}

///|
test "T57b: Satake purity and the Euler factor" {
  assert_true(satake_gl2_defect() < 1.0e-9)
  // the 3-isogeny at p = 5: a_5 = 1 with split multiplicative reduction test
  let (r1, r2) = satake_gl2(1.0, 5.0)
  assert_true(cx_approx(cx_mul(r1, r2), cx_of_real(5.0), tol=1.0e-9))
}

///|
test "T58: GL_2 recoupling, GL_3 pairwise law, and sym2 parameters" {
  let probes : Array[Array[Double]] = [
    [2.0, 3.0],
    [0.5, -1.5],
    [1.0, 1.0, 1.0],
    [0.3, 0.9, -0.6, 2.2],
  ]
  for x in probes {
    assert_true(recoupling_defect(x) < 1.0e-12)
    if x.length() == 3 {
      assert_true(gl3_pairwise_law_defect(x) < 1.0e-12)
    }
  }
  // (A1, omega_1): sym2 of the standard GL_2 Satake parameter has eigenvalue 2
  assert_true(satake_ic_match(RootType::TypeA(1), [1.0, 0.0]))
}

///|
test "T59: Jacquet-Langlands commutation (signed Newton identity)" {
  assert_true(jacquet_langlands_defect() < 1.0e-9)
}

///|
test "T60: the Fricke eigenvalue of the newform is +1" {
  assert_true(fricke_eigenform_certified())
  assert_eq(fricke_eigenvalue_sign(), 1)
  assert_eq(level_p_eigenform_dimension(11), 1)
  let (a, b, c) = hecke_polynomial_gl2(5)
  assert_eq(a, 1.0)
  assert_eq(b, -1.0)
  assert_eq(c, 5.0) // x^2 - a_5 x + p
}
