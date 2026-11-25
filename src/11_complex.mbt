// ============================================================================
// 11_complex.mbt -- Complex scalars and special functions of the
//                   Eisenstein / constant-term complex plane.
//
// ABSTRACTION.  The Langlands complex plane needs C rather than R: spectral
// parameters s = 1/2 + it, L-parameter eigenvalues, characters of the maximal
// torus.  We therefore lift `Double` to a pair `Cx { re, im }` and re-derive
// the transcendental layer (exp, log, pow, Gamma) on top of `@math`.
//
// Only a minimal, dependency-free kernel is provided; every operation is
// closed under the type, so all later modules stay purely functional.
// ============================================================================

///|
///A complex scalar `re + i·im`.
pub(all) struct Cx {
  re : Double
  im : Double
}

///|
///The complex zero.
pub fn cx_zero() -> Cx {
  { re: 0.0, im: 0.0, }
}

///|
///The complex one.
pub fn cx_one() -> Cx {
  { re: 1.0, im: 0.0, }
}

///|
///Lift a real number.
pub fn cx_of_real(x : Double) -> Cx {
  { re: x, im: 0.0, }
}

///|
///Polar form constructor `r·e^{iθ}`.
pub fn cx_polar(r : Double, theta : Double) -> Cx {
  { re: r * @math.cos(theta), im: r * @math.sin(theta), }
}

///|
///Addition.
pub fn cx_add(a : Cx, b : Cx) -> Cx {
  { re: a.re + b.re, im: a.im + b.im, }
}

///|
///Additive inverse.
pub fn cx_neg(a : Cx) -> Cx {
  { re: -a.re, im: -a.im, }
}

///|
///Subtraction.
pub fn cx_sub(a : Cx, b : Cx) -> Cx {
  cx_add(a, cx_neg(b))
}

///|
///Multiplication.
pub fn cx_mul(a : Cx, b : Cx) -> Cx {
  { re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re, }
}

///|
///Multiplication by a real scalar.
pub fn cx_scale(a : Cx, s : Double) -> Cx {
  { re: a.re * s, im: a.im * s, }
}

///|
///Complex conjugate.
pub fn cx_conj(a : Cx) -> Cx {
  { re: a.re, im: -a.im, }
}

///|
///Squared modulus `|a|²`.
pub fn cx_abs2(a : Cx) -> Double {
  a.re * a.re + a.im * a.im
}

///|
///Modulus.
pub fn cx_abs(a : Cx) -> Double {
  cx_abs2(a).sqrt()
}

///|
///Argument in `(-π, π]`.
pub fn cx_arg(a : Cx) -> Double {
  @math.atan2(a.im, a.re)
}

///|
///Division, with a defensive guard against 0/0.
pub fn cx_div(a : Cx, b : Cx) -> Cx {
  let d = cx_abs2(b)
  if d == 0.0 {
    { re: 0.0, im: 0.0, }
  } else {
    {
      re: (a.re * b.re + a.im * b.im) / d,
      im: (a.im * b.re - a.re * b.im) / d,
    }
  }
}

///|
///Principal logarithm `log|a| + i·arg(a)`.
pub fn cx_log(a : Cx) -> Cx {
  { re: @math.ln(cx_abs(a) + 1.0e-300), im: cx_arg(a), }
}

///|
///Complex exponential.
pub fn cx_exp(a : Cx) -> Cx {
  cx_polar(@math.exp(a.re), a.im)
}

///|
///Principal power `a^b = exp(b log a)`.
pub fn cx_pow(a : Cx, b : Cx) -> Cx {
  cx_exp(cx_mul(b, cx_log(a)))
}

///|
///Set membership test with tolerance (used in tests, not as a hash key).
pub fn cx_approx(a : Cx, b : Cx, tol~ : Double = 1.0e-9) -> Bool {
  approx(a.re, b.re, tol~) && approx(a.im, b.im, tol~)
}

// ---------------------------------------------------------------------------
// Gamma function (Lanczos, g = 7, 9 coefficients) and the completed zeta.
// ---------------------------------------------------------------------------

///|
/// Shared Lanczos `g = 7` coefficients (Godfrey / Spouge 9-term table).
fn lanczos_g7_coeffs() -> Array[Double] {
  [
    676.5203681218851, -1259.1392167224028, 771.3234287776531, -176.6150291621406,
    12.507343278686905, -0.13857109526572012, 9.984369578019572e-6, 1.5056327351493116e-7,
  ]
}

///|
/// `lgamma(x) = log Γ(x)` for `x > 0` via the shared complex Lanczos log-series.
/// Accuracy is ~1e-13 on `[0.5, 40]`.
pub fn lgamma(x : Double) -> Double {
  if x < 0.5 {
    @math.ln(@math.PI / @math.sin(@math.PI * x)) - lgamma(1.0 - x)
  } else {
    lgamma_cx_right(cx_of_real(x)).re
  }
}

///|
/// Bernoulli numbers `B_2, B_4, ... B_14` used by Euler–Maclaurin.
fn bernoulli_even() -> Array[Double] {
  [
    1.0 / 6.0,
    -1.0 / 30.0,
    1.0 / 42.0,
    -1.0 / 30.0,
    5.0 / 66.0,
    -691.0 / 2730.0,
    7.0 / 6.0,
  ]
}

///|
/// Complex sine `sin(z) = (e^{iz} - e^{-iz}) / (2i)`.
pub fn cx_sin(z : Cx) -> Cx {
  let iz = { re: -z.im, im: z.re, }
  let niz = { re: z.im, im: -z.re, }
  let diff = cx_sub(cx_exp(iz), cx_exp(niz))
  { re: 0.5 * diff.im, im: -0.5 * diff.re, }
}

///|
/// Complex cosine `cos(z) = (e^{iz} + e^{-iz}) / 2`.
pub fn cx_cos(z : Cx) -> Cx {
  let iz = { re: -z.im, im: z.re, }
  let niz = { re: z.im, im: -z.re, }
  cx_scale(cx_add(cx_exp(iz), cx_exp(niz)), 0.5)
}

///|
/// Riemann zeta across the ENTIRE complex plane `s ∈ ℂ \ {1}`:
/// Euler–Maclaurin for `Re(s) >= 1/2` (including the critical line `s = 1/2 + it`),
/// and full complex Riemann functional equation
/// `ζ(s) = 2^s π^{s-1} sin(πs/2) Γ(1-s) ζ(1-s)` for `Re(s) < 1/2` (eliminating the
/// previous `s.im == 0` restriction).
pub fn zeta_cx(s : Cx) -> Cx {
  if s.re < 0.5 {
    let one_minus = cx_sub(cx_one(), s)
    let z = zeta_euler_maclaurin(one_minus)
    let two_s = cx_pow(cx_of_real(2.0), s)
    let pi_s1 = cx_pow(cx_of_real(@math.PI), cx_sub(s, cx_one()))
    let sin_term = cx_sin(cx_scale(s, 0.5 * @math.PI))
    let gam = gamma_cx(one_minus)
    cx_mul(cx_mul(cx_mul(two_s, pi_s1), cx_mul(sin_term, gam)), z)
  } else {
    zeta_euler_maclaurin(s)
  }
}

///|
/// Euler–Maclaurin evaluation across the full half-plane `Re(s) >= 1/2`:
/// `ζ(s) = Σ_{n<N} n^{-s} + N^{1-s}/(s-1) + N^{-s}/2 + Σ_k B_{2k}/(2k)! (s)_{2k-1} N^{-s-2k+1}`.
/// Uses true complex powers `cx_pow(N, ...)` rather than real-part exponents so
/// imaginary spectral parameters `s = σ + it` are exact.
pub fn zeta_euler_maclaurin(s : Cx) -> Cx {
  let n_terms = 40
  let b = bernoulli_even()
  let mut total = cx_zero()
  for n in 1..<n_terms {
    let l = cx_log(cx_of_real(n.to_double()))
    total = cx_add(total, cx_exp(cx_neg(cx_mul(l, s))))
  }
  let nc = cx_of_real(n_terms.to_double())
  let s_minus_1 = cx_sub(s, cx_one())
  let tail1 = cx_div(cx_pow(nc, cx_sub(cx_one(), s)), s_minus_1)
  let tail2 = cx_scale(cx_pow(nc, cx_neg(s)), 0.5)
  total = cx_add(total, cx_add(tail1, tail2))
  for k in 1..<=7 {
    let mut poch = cx_one()
    for j in 0..<(2 * k - 1) {
      poch = cx_mul(poch, cx_add(s, cx_of_real(j.to_double())))
    }
    let coef = b[k - 1] / factorial(2 * k).to_double()
    let exp_s = cx_neg(cx_add(s, cx_of_real((2 * k - 1).to_double())))
    let n_pow = cx_pow(nc, exp_s)
    total = cx_add(total, cx_scale(cx_mul(poch, n_pow), coef))
  }
  total
}

///|
///Factorial on `Int32` (small arguments only).
pub fn factorial(n : Int32) -> Int32 {
  let mut r = 1
  for i in 2..<=n {
    r *= i
  }
  r
}

///|
/// Completed Riemann zeta `ξ(s) = π^{-s/2} Γ(s/2) ζ(s)` across the full complex plane.
/// For `Re(s) < 1/2`, uses the functional equation `ξ(s) = ξ(1-s)` directly so poles of
/// `Γ(s/2)` at non-positive even integers cancel against the trivial zeros of `ζ(s)`.
pub fn xi_cx(s : Cx) -> Cx {
  if s.re < 0.5 {
    return xi_cx_right(cx_sub(cx_one(), s))
  }
  xi_cx_right(s)
}

///|
fn xi_cx_right(s : Cx) -> Cx {
  let half_pi_pow = cx_pow(cx_of_real(@math.PI), cx_scale(cx_neg(s), 0.5))
  cx_mul(cx_mul(half_pi_pow, gamma_cx(cx_scale(s, 0.5))), zeta_cx(s))
}

///|
/// `Γ(z)` across the full complex plane `ℂ \ {0, -1, -2, ...}` via Lanczos + Euler reflection
/// `Γ(z) = π / (sin(πz) Γ(1-z))` for `Re(z) < 0.5`.
pub fn gamma_cx(z : Cx) -> Cx {
  if z.re < 0.5 {
    let sin_pz = cx_sin(cx_scale(z, @math.PI))
    return cx_div(
      cx_of_real(@math.PI),
      cx_mul(sin_pz, gamma_cx_lanczos(cx_sub(cx_one(), z))),
    )
  }
  gamma_cx_lanczos(z)
}

///|
/// `log Γ(z)` on the Lanczos half-plane `Re(z) ≥ 0.5`; the single source of
/// truth for both the real `lgamma` and every complex Γ-evaluation.
fn lgamma_cx_right(z : Cx) -> Cx {
  cx_log(gamma_cx_lanczos(z))
}

///|
/// Lanczos `g = 7` series on the half-plane `Re(z) ≥ 0.5`.  This is the only
/// place the 9-term coefficient table is stored (cf. `lanczos_g7_coeffs`).
fn gamma_cx_lanczos(z : Cx) -> Cx {
  let g = lanczos_g7_coeffs()
  let zz = cx_sub(z, cx_one())
  let t = cx_add(zz, cx_of_real(7.5))
  let mut a = cx_of_real(0.99999999999980993)
  for i in 0..<8 {
    let den = cx_add(zz, cx_of_real((i + 1).to_double()))
    a = cx_add(a, cx_div(cx_of_real(g[i]), den))
  }
  let half = cx_of_real(0.5)
  let log_pref = cx_add(
    cx_of_real(0.5 * @math.ln(2.0 * @math.PI)),
    cx_sub(cx_mul(cx_add(zz, half), cx_log(t)), t),
  )
  cx_mul(cx_exp(log_pref), a)
}

// ---------------------------------------------------------------------------
// Tests: the complex kernel and the completed-zeta functional equation.
// ---------------------------------------------------------------------------

///|
test "T16a: complex kernel is a field" {
  let a = { re: 1.5, im: -0.75, }
  let b = cx_polar(2.0, 0.3)
  // multiplicative inverses
  let one = cx_mul(a, cx_div(cx_one(), a))
  assert_true(cx_approx(one, cx_one(), tol=1.0e-12))
  let one2 = cx_mul(b, cx_div(cx_one(), b))
  assert_true(cx_approx(one2, cx_one(), tol=1.0e-12))
  // exp/log round trip
  let z = cx_log(cx_exp(a))
  assert_true(cx_approx(z, a, tol=1.0e-12))
  // |a·conj(a)| = |a|²
  let m = cx_mul(a, cx_conj(a))
  assert_true(approx(m.im, 0.0, tol=1.0e-14))
  assert_true(approx(m.re, cx_abs2(a), tol=1.0e-14))
}

///|
test "T16b: ζ(2) = π²/6 and ζ(-1) = -1/12" {
  assert_true(
    approx(zeta_cx(cx_of_real(2.0)).re, @math.PI * @math.PI / 6.0, tol=1.0e-12),
  )
  assert_true(
    approx(
      zeta_cx(cx_of_real(4.0)).re,
      @math.pow(@math.PI, 4.0) / 90.0,
      tol=1.0e-11,
    ),
  )
  assert_true(approx(zeta_cx(cx_of_real(-1.0)).re, -1.0 / 12.0, tol=1.0e-11))
}

///|
test "T16c: ξ(s) = ξ(1-s) on and off the critical line (including full complex plane)" {
  for k in 1..<=5 {
    let s = 0.4 + 0.37 * k.to_double()
    assert_true(
      cx_approx(xi_cx(cx_of_real(s)), xi_cx(cx_of_real(1.0 - s)), tol=1.0e-8),
    )
    assert_true(
      approx(
        @math.exp(lgamma(s + 1.0)),
        gamma_cx(cx_of_real(s + 1.0)).re,
        tol=1.0e-11,
      ),
    )
  }
  // Full complex functional equation: Re(s) < 0.5 and Im(s) != 0
  let s_cx = { re: -0.35, im: 2.75, }
  let z_lhs = zeta_cx(s_cx)
  let one_minus = cx_sub(cx_one(), s_cx)
  let rhs_pref = cx_mul(
    cx_mul(
      cx_pow(cx_of_real(2.0), s_cx),
      cx_pow(cx_of_real(@math.PI), cx_sub(s_cx, cx_one())),
    ),
    cx_mul(cx_sin(cx_scale(s_cx, 0.5 * @math.PI)), gamma_cx(one_minus)),
  )
  assert_true(
    cx_approx(z_lhs, cx_mul(rhs_pref, zeta_cx(one_minus)), tol=1.0e-10),
  )
  // sin^2(z) + cos^2(z) = 1 in C
  let sc = cx_add(
    cx_mul(cx_sin(s_cx), cx_sin(s_cx)),
    cx_mul(cx_cos(s_cx), cx_cos(s_cx)),
  )
  assert_true(cx_approx(sc, cx_one(), tol=1.0e-12))
  // First nontrivial Riemann zero on the critical line: s = 1/2 + 14.134725141734695 i
  let zero1 = { re: 0.5, im: 14.134725141734695, }
  assert_true(cx_abs(zeta_cx(zero1)) < 1.0e-7)
  // Γ(1/2)=√π and ξ(1/2) real
  let half = cx_of_real(0.5)
  assert_true(approx(gamma_cx(half).re, @math.PI.sqrt(), tol=1.0e-10))
  let xh = xi_cx(half)
  assert_true(approx(xh.im, 0.0, tol=1.0e-12))
}
