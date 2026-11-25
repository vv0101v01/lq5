// ============================================================================
// 12_eisenstein.mbt -- Eisenstein series and the constant-term functor.
//
// ABSTRACTION.
//   Mathematical object        Computer abstraction
//   -------------------------  ------------------------------------------------
//   weight λ ∈ a*_C           Array[Double] (rank r = dim of the torus)
//   Weyl group W (type A)     Array of permutations of 0..<r
//   intertwining operator     matrix on the W-indexed coefficient space
//   Eisenstein series Eis(f)  Σ_{w∈W} f_w y^{wλ}: the *section* functor
//   constant term CT           Σ_{w∈W} c(w,λ) P_w - Langlands' formula
//   functional equation        CT(w0)∘CT(w0) = Id  (⇔ c(w0,λ)c(w0,w0λ)=1)
//
// Everything is written for the minimal parabolic of GL_r / SL_r, i.e. Weyl
// group S_r, so the constant-term formula becomes a finite matrix identity
// that a test can *decide* rather than merely illustrate.
//
// Normalisation.  With ξ(s) = π^{-s/2}Γ(s/2)ζ(s) the SL_2 constant term of
//   E(z,s) = Σ_{γ∈Γ_∞\Γ} (Im γz)^s = (1/2) Σ_{(m,n)=1} y^s/|mz+n|^{2s}
// is  A_0(y,s) = y^s + φ(s) y^{1-s},  φ(s) = ξ(2s-1)/ξ(2s),
// and its residue at s = 1 is 1/ξ(2) = 3/π (Theorem T19).
// ============================================================================

///|
///Euler-φ like helper: greatest common divisor.
pub fn gcd_int(a : Int32, b : Int32) -> Int32 {
  let mut x = if a < 0 { -a } else { a }
  let mut y = if b < 0 { -b } else { b }
  while y != 0 {
    let t = x % y
    x = y
    y = t
  }
  x
}

// ---------------------------------------------------------------------------
// The local factor c(λ) = ξ(λ)/ξ(λ+1) of the c-function.
// ---------------------------------------------------------------------------

///|
/// The SL_2 local intertwining factor `c(t) = ξ(t)/ξ(t+1)`.
/// Laws (Theorem T17):
///   * `c(t)·c(-t) = 1`  (functional equation, local form);
///   * `c(0) = -1`       (both ξ(t) and ξ(t+1) have a pole at 0; the limit is -1);
///   * `c` has a simple pole at `t = 1`, with `(t-1)c(t) → 1/ξ(2) = 3/π`.
/// The pole at `t = 1` is the residual spectrum: a positive root with
/// `⟨λ,α^∨⟩ = 1` produces a non-zero residue of the Eisenstein series.
pub fn c_factor(t : Double) -> Double {
  let dist = (t - 1.0).abs()
  if dist < 1.0e-6 {
    // Sentinels for the pole: sign preserved, magnitude large but finite so
    // that downstream arithmetic stays in the field of Double.
    return if t > 1.0 { 1.0e6 } else { -1.0e6 }
  }
  let num = xi_cx(cx_of_real(t)).re
  let den = xi_cx(cx_of_real(t + 1.0)).re
  num / den
}

///|
/// Reflection defect `|c(t)c(-t) - 1|`: the numerical form of the
/// functional equation (Theorem T17).
pub fn c_reflection_defect(t : Double) -> Double {
  (c_factor(t) * c_factor(-t) - 1.0).abs()
}

///|
/// Pole witness `(t-1)·c(t)`, which converges to `1/ξ(2) = 3/π` as `t → 1`.
pub fn c_pole_witness(t : Double) -> Double {
  (t - 1.0) * c_factor(t)
}

// ---------------------------------------------------------------------------
// Weyl group of type A_{r-1} and the Gindikin-Karpelevich c-function.
// ---------------------------------------------------------------------------

///|
/// All permutations of `0..<r` (the Weyl group S_r), in lexicographic order.
pub fn weyl_group(r : Int32) -> Array[Array[Int32]] {
  let out : Array[Array[Int32]] = []
  let cur : Array[Int32] = []
  for i in 0..<r {
    cur.push(i)
  }
  permute_into(cur, 0, out)
  out
}

///|
///Recursive permutation generator (Heap-free, lexicographic).
fn permute_into(
  cur : Array[Int32],
  k : Int32,
  out : Array[Array[Int32]],
) -> Unit {
  let n = cur.length()
  if k == n {
    let copy : Array[Int32] = []
    for i in 0..<n {
      copy.push(cur[i])
    }
    out.push(copy)
    return
  }
  for i in k..<n {
    let tmp = cur[k]
    cur[k] = cur[i]
    cur[i] = tmp
    permute_into(cur, k + 1, out)
    let tmp2 = cur[k]
    cur[k] = cur[i]
    cur[i] = tmp2
  }
}

///|
/// Product of two permutations: `(w1·w2)(i) = w1(w2(i))`.
pub fn weyl_mul(w1 : Array[Int32], w2 : Array[Int32]) -> Array[Int32] {
  let out : Array[Int32] = []
  for i in 0..<w1.length() {
    out.push(w1[w2[i]])
  }
  out
}

///|
/// Inverse permutation.
pub fn weyl_inv(w : Array[Int32]) -> Array[Int32] {
  let n = w.length()
  let out : Array[Int32] = []
  for _ in 0..<n {
    out.push(0)
  }
  for i in 0..<n {
    out[w[i]] = i
  }
  out
}

///|
/// Action on weights: `(w·λ)_i = λ_{w^{-1}(i)}`.
pub fn weyl_act(w : Array[Int32], lambda : Array[Double]) -> Array[Double] {
  let n = w.length()
  let winv = weyl_inv(w)
  let out : Array[Double] = []
  for i in 0..<n {
    out.push(lambda[winv[i]])
  }
  out
}

///|
/// The longest element `w_0` (reversal): it sends every positive root to a
/// negative one, so `w_0·λ = -λ` up to a Weyl translate.
pub fn weyl_longest(r : Int32) -> Array[Int32] {
  let out : Array[Int32] = []
  for i in 0..<r {
    out.push(r - 1 - i)
  }
  out
}

///|
/// Gindikin-Karpelevich factor
/// `c(w,λ) = Π_{i<j, w(i)>w(j)} c(λ_i - λ_j)`.
/// The condition `w(i) > w(j)` is the computational form of `w α < 0`
/// for the positive root `α = e_i - e_j` (Theorem T18).
pub fn c_function(w : Array[Int32], lambda : Array[Double]) -> Double {
  let n = w.length()
  let mut t = 1.0
  for i in 0..<n {
    for j in (i + 1)..<n {
      if w[i] > w[j] {
        t *= c_factor(lambda[i] - lambda[j])
      }
    }
  }
  t
}

///|
/// Cocycle defect of the c-function:
/// `|c(w1·w2, λ) - c(w1, w2·λ)·c(w2, λ)|` (Theorem T18).  Vanishing for all
/// pairs means that `w ↦ c(w,·)` is a crossed homomorphism, which is what
/// lets `CT` be a *functor* (composition of parabolics) rather than a
/// function of the individual `w`.
pub fn c_cocycle_defect(
  w1 : Array[Int32],
  w2 : Array[Int32],
  lambda : Array[Double],
) -> Double {
  let lhs = c_function(weyl_mul(w1, w2), lambda)
  let rhs = c_function(w1, weyl_act(w2, lambda)) * c_function(w2, lambda)
  (lhs - rhs).abs()
}

///|
/// Functional-equation defect `|c(w0,λ)c(w0,w0λ) - 1|`.
pub fn c_functional_eq_defect(lambda : Array[Double]) -> Double {
  let r = lambda.length()
  let w0 = weyl_longest(r)
  (c_function(w0, lambda) * c_function(w0, weyl_act(w0, lambda)) - 1.0).abs()
}

// ---------------------------------------------------------------------------
// The constant-term functor on the W-indexed coefficient space.
// ---------------------------------------------------------------------------

///|
/// A minimal-parabolic Eisenstein datum: rank `r` and spectral parameter
/// `λ ∈ C^r`.  `sections` are coefficient vectors indexed by `W`: a formal
/// sum `Σ_w a_w y^{wλ}`.
pub(all) struct EisensteinDatum {
  rank : Int32
  lambda : Array[Double]
}

///|
/// Constructor.
pub fn eisenstein_datum(lambda : Array[Double]) -> EisensteinDatum {
  { rank: lambda.length(), lambda, }
}

///|
/// The constant-term matrix of the minimal parabolic:
///   `A(λ)[u][v] = c(u·v^{-1}, v·λ)`.
/// Interpretation: the entry says how much of the exponent `vλ` survives,
/// after averaging over the unipotent radical, into the exponent `uλ`.
/// Diagonal entries are `c(1, vλ) = 1` (the identity contributes nothing).
pub fn EisensteinDatum::constant_term_matrix(
  self : EisensteinDatum,
) -> Array[Double] {
  let ws = weyl_group(self.rank)
  let n = ws.length()
  let m = Array::make(n * n, 0.0)
  for a in 0..<n {
    let u = ws[a]
    for b in 0..<n {
      let v = ws[b]
      let w = weyl_mul(u, weyl_inv(v))
      let mu = weyl_act(v, self.lambda)
      m[a * n + b] = c_function(w, mu)
    }
  }
  m
}

///|
/// The Eisenstein spectral projection `Π(λ) = A(λ)/|W|`.
/// Theorem T20: `Π² = Π`, `tr Π = 1`: the constant-term functor carries a
/// canonical rank-one projection onto the Eisenstein line of the spectrum.
pub fn EisensteinDatum::spectral_projection(
  self : EisensteinDatum,
) -> Array[Double] {
  let n = weyl_group(self.rank).length()
  let a = self.constant_term_matrix()
  let out : Array[Double] = []
  for i in 0..<(n * n) {
    out.push(a[i] / n.to_double())
  }
  out
}

///|
/// `‖Π² - Π‖_∞`: the closure certificate of Theorem T20.
pub fn EisensteinDatum::projection_defect(self : EisensteinDatum) -> Double {
  let n = weyl_group(self.rank).length()
  let p = self.spectral_projection()
  let p2 = matrix_product(p, n, n, p, n)
  vec_inf_diff(p2, p)
}

///|
/// `tr Π`, which must equal `1` = the multiplicity of the Eisenstein
/// representation in the spectrum (Theorem T20b).
pub fn EisensteinDatum::projection_trace(self : EisensteinDatum) -> Double {
  let n = weyl_group(self.rank).length()
  let p = self.spectral_projection()
  let mut tr = 0.0
  for i in 0..<n {
    tr += p[i * n + i]
  }
  tr
}

///|
/// The scattering/constant-term determinant defect:
/// `|det A(λ)|` must vanish, since `A = |W|·Π` with `Π` of rank one.
pub fn EisensteinDatum::scattering_det(self : EisensteinDatum) -> Double {
  let n = weyl_group(self.rank).length()
  determinant(self.constant_term_matrix(), n)
}

// ---------------------------------------------------------------------------
// The archimedean SL_2 Eisenstein series: numerics.
// ---------------------------------------------------------------------------

///|
/// Truncated Eisenstein series
/// `E_x,y(z,s) = (1/2) Σ_{gcd(m,n)=1, |m|,|n|≤N} y^s / |mz+n|^{2s}`,
/// `z = x + iy`.  Converges for `Re s > 1`; the truncation error is `O(N^{1-2s})`.
pub fn eisenstein_sl2(
  x : Double,
  y : Double,
  s : Double,
  n_terms : Int32,
) -> Double {
  let mut total = 0.0
  for m in -n_terms..<=n_terms {
    for n in -n_terms..<=n_terms {
      if m == 0 && n == 0 {
        continue
      }
      if gcd_int(m, n) != 1 {
        continue
      }
      let re = m.to_double() * x + n.to_double()
      let im = m.to_double() * y
      let abs2 = re * re + im * im
      total += @math.pow(y, s) / @math.pow(abs2, s)
    }
  }
  0.5 * total
}

///|
/// Langlands' constant term for `SL_2`:
/// `A_0(y,s) = y^s + φ(s) y^{1-s}`, `φ(s) = ξ(2s-1)/ξ(2s) = c(2s-1)`.
pub fn constant_term_sl2(y : Double, s : Double) -> Double {
  @math.pow(y, s) + c_factor(2.0 * s - 1.0) * @math.pow(y, 1.0 - s)
}

///|
/// The same quantity obtained by *integrating* the series over `x ∈ [0,1]`,
/// i.e. by applying the constant-term functor numerically.  This is the
/// computation that Langlands performed by hand in 1966.
pub fn constant_term_numeric(
  y : Double,
  s : Double,
  n_terms : Int32,
  n_quad : Int32,
) -> Double {
  let q = gauss_legendre(n_quad)
  let mut acc = 0.0
  for i in 0..<q.nodes.length() {
    // map [-1,1] -> [0,1]
    let x = 0.5 * (q.nodes[i] + 1.0)
    acc += q.weights[i] * 0.5 * eisenstein_sl2(x, y, s, n_terms)
  }
  acc
}

///|
/// `|φ(s)φ(1-s) - 1|`: the functional equation of the intertwining operator,
/// which is nothing but the local reflection law `c(t)c(-t) = 1` at `t = 2s-1`.
pub fn phi_functional_defect(s : Double) -> Double {
  c_reflection_defect(2.0 * s - 1.0)
}

///|
/// Residue of `E(z,s)` at `s = 1`, predicted by the pole of `c` at `t = 1`:
/// `lim_{s→1}(s-1)φ(s) = ½·lim_{t→1}(t-1)c(t) = 1/(2ξ(2)) = 3/π`
/// (Theorem T19b).  Expressed through `c_pole_witness` so that the residue
/// exists in exactly one place in the code base.
pub fn residue_at_one(eps~ : Double = 1.0e-6) -> Double {
  0.5 * c_pole_witness(1.0 + 2.0 * eps)
}

///|
/// Laplace eigenvalue check on the constant term: for `Δ = y²(∂_x²+∂_y²)`
/// we must have `Δ A_0 = s(s-1) A_0` (Theorem T19b).  Computed by central
/// differences with the exact closed-form `A_0`.
pub fn constant_term_laplacian(
  y : Double,
  s : Double,
  h~ : Double = 1.0e-4,
) -> Double {
  let f = fn(_x, yy) { constant_term_sl2(yy, s) }
  let dxx = (f(0.0 + h, y) - 2.0 * f(0.0, y) + f(0.0 - h, y)) / (h * h)
  let dyy = (f(0.0, y + h) - 2.0 * f(0.0, y) + f(0.0, y - h)) / (h * h)
  y * y * (dxx + dyy)
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T17: c(t)·c(-t) = 1 (Eisenstein functional equation, local form)" {
  for k in 1..<=6 {
    let t = 0.31 * k.to_double()
    assert_true(c_reflection_defect(t) < 1.0e-8)
    assert_true(c_reflection_defect(-t) < 1.0e-8)
  }
  // c(0) = -1 (removable crossing of the two poles) and the *simple* pole of
  // c sits at t = 1, with residue 3/π.
  assert_true((c_factor(1.0e-3) + 1.0).abs() < 5.0e-3)
  // (t-1)c(t) -> 1/ξ(2) = 6/π; in the s-variable (t = 2s-1) this halves to 3/π.
  assert_true((c_pole_witness(1.0 + 1.0e-4) - 6.0 / @math.PI).abs() < 5.0e-3)
}

///|
test "T18: Gindikin-Karpelevich c-function is a crossed homomorphism" {
  let lambda : Array[Double] = [1.7, -0.4, 2.3]
  let ws = weyl_group(3)
  for a in 0..<ws.length() {
    for b in 0..<ws.length() {
      assert_true(c_cocycle_defect(ws[a], ws[b], lambda) < 1.0e-8)
    }
  }
  // w0 sends every positive root to a negative one: full product.
  let full = c_factor(1.7 + 0.4) * c_factor(1.7 - 2.3) * c_factor(-0.4 - 2.3)
  assert_true(approx(c_function(weyl_longest(3), lambda), full, tol=1.0e-9))
  // rank 4 functional equation
  let lam4 : Array[Double] = [0.9, -1.3, 2.1, 0.4]
  assert_true(c_functional_eq_defect(lam4) < 1.0e-8)
}

///|
test "T19: Langlands constant term computed by integration" {
  let y = 1.37
  let s = 2.0
  let numeric = constant_term_numeric(y, s, 40, 40)
  let formula = constant_term_sl2(y, s)
  assert_true((numeric - formula).abs() < 2.0e-3)
  // residue at s = 1 equals 3/π = 1/ξ(2)
  assert_true((residue_at_one() - 3.0 / @math.PI).abs() < 1.0e-3)
  // φ(s)φ(1-s) = 1
  assert_true(phi_functional_defect(1.7) < 1.0e-8)
  // Δ A_0 = s(1-s) A_0
  for k in 1..<=3 {
    let ss = 1.3 + 0.4 * k.to_double()
    let lhs = constant_term_laplacian(1.21, ss)
    let rhs = ss * (ss - 1.0) * constant_term_sl2(1.21, ss)
    assert_true((lhs - rhs).abs() < 1.0e-4 * (1.0 + rhs.abs()))
  }
}

///|
test "T20: constant-term matrix is |W| times a rank-one projection" {
  let lam2 : Array[Double] = [1.4, -0.9]
  let lam3 : Array[Double] = [1.7, -0.4, 2.3]
  let lam4 : Array[Double] = [0.9, -1.3, 2.1, 0.4]
  for lam in [lam2, lam3, lam4] {
    let d = eisenstein_datum(lam)
    let n = weyl_group(d.rank).length()
    assert_true(d.projection_defect() < 1.0e-8)
    assert_true((d.projection_trace() - 1.0).abs() < 1.0e-9)
    assert_true(d.scattering_det().abs() < 1.0e-8)
    // A itself is not a projection: A² = |W|·A ≠ A.
    let a = d.constant_term_matrix()
    let a2 = matrix_product(a, n, n, a, n)
    assert_true(vec_inf_diff(a2, vec_scale(a, n.to_double())) < 1.0e-8)
    // diagonal amplitudes are 1
    for i in 0..<n {
      assert_true((a[i * n + i] - 1.0).abs() < 1.0e-12)
    }
  }
}
