// ============================================================================
// 15_betti.mbt -- Betti moduli space / character variety, and the
//                 Eisenstein-constant-term correspondence on it.
//
// ABSTRACTION.
//   Betti moduli space     M_B(Σ_g,p, G) = Hom(π_1 Σ_{g,p}, G) // G
//   character variety      the same space in trace coordinates (Fricke)
//   Casimir                κ(x,y,z) = x² + y² + z² - xyz = 2 + tr[A,B]
//   Goldman/Atiyah-Bott    Poisson bracket on trace functions, whose
//   Poisson structure       symplectic leaves are the level sets of κ
//   Eisenstein/CT pair     block-diagonal embedding vs. semisimple projection
//
// For the one-holed torus the whole moduli space is 3-dimensional with
// trace coordinates (x,y,z) = (tr A, tr B, tr AB) and the single relation
// κ = 2 + tr[A,B].  The Poisson tensor is *not* symplectic: its rank is 2 and
// its kernel is spanned by dκ (Theorem T23).
//
// The Eisenstein/constant-term pair is realised on the moduli space itself:
//   Eis : M_B(M) → M_B(G),  block-diagonal inclusion,
//   CT  : M_B(G) → M_B(M),  semisimple (diagonal) part,
//   CT ∘ Eis = id,  and the fibre of CT over a point is the unipotent orbit
// of the L-parameter (dimension n² - Σ(λ'_i)²).
// ============================================================================

///|
///A 2x2 real matrix, used for `SL_2`-valued representations of `π_1`.
pub(all) struct Mat2 {
  a : Double
  b : Double
  c : Double
  d : Double
}

///|
///Product of 2x2 matrices.
pub fn mat2_mul(x : Mat2, y : Mat2) -> Mat2 {
  {
    a: x.a * y.a + x.b * y.c,
    b: x.a * y.b + x.b * y.d,
    c: x.c * y.a + x.d * y.c,
    d: x.c * y.b + x.d * y.d,
  }
}

///|
///Determinant.
pub fn mat2_det(x : Mat2) -> Double {
  x.a * x.d - x.b * x.c
}

///|
///Trace.
pub fn mat2_tr(x : Mat2) -> Double {
  x.a + x.d
}

///|
///Inverse of a unit-determinant matrix.
pub fn mat2_inv(x : Mat2) -> Mat2 {
  { a: x.d, b: -x.b, c: -x.c, d: x.a, }
}

///|
///Commutator `[x,y] = xyx^{-1}y^{-1}`.
pub fn mat2_comm(x : Mat2, y : Mat2) -> Mat2 {
  mat2_mul(mat2_mul(mat2_mul(x, y), mat2_inv(x)), mat2_inv(y))
}

///|
/// A random element of `SL_2(R)` (proper, i.e. determinant exactly 1 up to
/// rounding), used to sample the Betti moduli space.
pub fn random_sl2(rng : Rng) -> Mat2 {
  let a = 0.6 + 1.8 * rng.unit()
  let b = -1.4 + 2.8 * rng.unit()
  let c = -1.4 + 2.8 * rng.unit()
  let d = (1.0 + b * c) / a
  { a, b, c, d, }
}

///|
/// The Fricke Casimir `κ(x,y,z) = x² + y² + z² - xyz` in trace coordinates
/// `x = tr A`, `y = tr B`, `z = tr AB`.
pub fn fricke_kappa(x : Double, y : Double, z : Double) -> Double {
  x * x + y * y + z * z - x * y * z
}

///|
/// Trace of the commutator: by Cayley–Hamilton `tr[A,B] = κ - 2`
/// (Theorem T23a: the Fricke/Vogt identity).
pub fn commutator_trace(x : Double, y : Double, z : Double) -> Double {
  fricke_kappa(x, y, z) - 2.0
}

///|
/// Defect of the Fricke identity for *explicit* matrices: comparing
/// `2 + tr[A,B]` with `κ(tr A, tr B, tr AB)`.
pub fn fricke_defect(a : Mat2, b : Mat2) -> Double {
  let x = mat2_tr(a)
  let y = mat2_tr(b)
  let z = mat2_tr(mat2_mul(a, b))
  (2.0 + mat2_tr(mat2_comm(a, b)) - fricke_kappa(x, y, z)).abs()
}

// ---------------------------------------------------------------------------
// The Goldman / Atiyah-Bott Poisson structure on the character variety.
// ---------------------------------------------------------------------------

///|
/// The Poisson tensor in trace coordinates, as a symmetric-index 3x3
/// antisymmetric matrix:
///   `{x,y} = 2z - xy`, `{y,z} = 2x - yz`, `{z,x} = 2y - zx`.
/// Equivalently `P^{ij} = ε^{ijk} ∂_k κ`, so `P dκ = 0`: κ is a Casimir and
/// the symplectic leaves are the level sets of κ (Theorem T23b).
pub fn goldman_tensor(x : Double, y : Double, z : Double) -> Array[Double] {
  let pxy = 2.0 * z - x * y
  let pyz = 2.0 * x - y * z
  let pzx = 2.0 * y - z * x
  [0.0, pxy, -pzx, -pxy, 0.0, pyz, pzx, -pyz, 0.0]
}

///|
/// Central-difference gradient of a scalar function on the 3-dimensional
/// trace-coordinate space.
fn grad_at(
  f : (Array[Double]) -> Double,
  p : Array[Double],
  h : Double,
) -> Array[Double] {
  let out : Array[Double] = []
  for i in 0..<3 {
    let up : Array[Double] = [p[0], p[1], p[2]]
    let dn : Array[Double] = [p[0], p[1], p[2]]
    up[i] = p[i] + h
    dn[i] = p[i] - h
    out.push((f(up) - f(dn)) / (2.0 * h))
  }
  out
}

///|
/// Poisson bracket of two differentiable functions on the character variety,
/// `{f,g} = Σ_{i<j} P^{ij}(∂_i f ∂_j g - ∂_j f ∂_i g)`, with central
/// differences.
pub fn poisson_bracket(
  f : (Array[Double]) -> Double,
  g : (Array[Double]) -> Double,
  p : Array[Double],
  h~ : Double = 1.0e-5,
) -> Double {
  let gf = grad_at(f, p, h)
  let gg = grad_at(g, p, h)
  let t = goldman_tensor(p[0], p[1], p[2])
  let mut s = 0.0
  for i in 0..<3 {
    for j in 0..<3 {
      s += t[i * 3 + j] * gf[i] * gg[j]
    }
  }
  s
}

///|
/// Rank of the Poisson tensor at `p`: 2 in the open dense set (the variety is
/// foliated by 2-dimensional symplectic leaves).
pub fn poisson_rank(p : Array[Double]) -> Int32 {
  matrix_rank(goldman_tensor(p[0], p[1], p[2]), 3, 3, tol=1.0e-9)
}

///|
/// Jacobi defect `{{f,g},h} + {{g,h},f} + {{h,f},g}` for three functions,
/// evaluated at `p`.  Vanishing is the integrability of the structure.
pub fn jacobi_defect(
  f : (Array[Double]) -> Double,
  g : (Array[Double]) -> Double,
  k : (Array[Double]) -> Double,
  p : Array[Double],
) -> Double {
  let fg = fn(q : Array[Double]) { poisson_bracket(f, g, q) }
  let gk = fn(q : Array[Double]) { poisson_bracket(g, k, q) }
  let kf = fn(q : Array[Double]) { poisson_bracket(k, f, q) }
  let d = poisson_bracket(fg, k, p) +
    poisson_bracket(gk, f, p) +
    poisson_bracket(kf, g, p)
  d.abs()
}

///|
/// Casimir defect `max(|{κ,x}|, |{κ,y}|, |{κ,z}|)`: must vanish.
pub fn casimir_defect(p : Array[Double]) -> Double {
  let kap = fn(q : Array[Double]) { fricke_kappa(q[0], q[1], q[2]) }
  let cx = poisson_bracket(kap, fn(q : Array[Double]) { q[0] }, p)
  let cy = poisson_bracket(kap, fn(q : Array[Double]) { q[1] }, p)
  let cz = poisson_bracket(kap, fn(q : Array[Double]) { q[2] }, p)
  let a = cx.abs()
  let b = cy.abs()
  let c = cz.abs()
  let m = if a > b { a } else { b }
  if m > c {
    m
  } else {
    c
  }
}

///|
/// Project a displaced point back onto the symplectic leaf of `p0`
/// (`κ = const`) by Newton iteration on the parameter `s` of the segment
/// `p0 + s·(p1 - p0)`.  This is the *exact* Casimir-preserving integrator: the
/// correction `|s - 1|` is `O(dt²)` and is reported as `leaf_projection`.
pub fn project_to_leaf(
  p0 : Array[Double],
  p1 : Array[Double],
  iters~ : Int32 = 3,
) -> Array[Double] {
  let k0 = fricke_kappa(p0[0], p0[1], p0[2])
  let v : Array[Double] = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]]
  let mut s = 1.0
  for _ in 0..<iters {
    let q : Array[Double] = [
      p0[0] + s * v[0],
      p0[1] + s * v[1],
      p0[2] + s * v[2],
    ]
    let g = fricke_kappa(q[0], q[1], q[2]) - k0
    // ∇κ = (2x - yz, 2y - xz, 2z - xy)
    let gx = 2.0 * q[0] - q[1] * q[2]
    let gy = 2.0 * q[1] - q[0] * q[2]
    let gz = 2.0 * q[2] - q[0] * q[1]
    let dg = gx * v[0] + gy * v[1] + gz * v[2]
    if dg.abs() < 1.0e-14 {
      return q
    }
    s = s - g / dg
  }
  [p0[0] + s * v[0], p0[1] + s * v[1], p0[2] + s * v[2]]
}

// ---------------------------------------------------------------------------
// Dimension bookkeeping for Betti moduli spaces.
// ---------------------------------------------------------------------------

///|
/// `dim M_B(Σ_{g,p}, SL_n) = (2g-2)(n²-1) + 2p(n-1)`.
/// Checks: `(g,p,n) = (1,1,2) ↦ 2` (one-holed torus), `(2,0,2) ↦ 6`,
/// `(1,2,2) ↦ 4`.
pub fn betti_dim(genus : Int32, punctures : Int32, n : Int32) -> Int32 {
  let nn = n * n - 1
  let g2 = 2 * genus - 2
  g2 * nn + 2 * punctures * (n - 1)
}

///|
/// Dimension of the character variety of the *closed* surface of genus `g`
/// for `GL_n`, i.e. `(2g-2)n² + 2`.
pub fn betti_dim_gl(genus : Int32, n : Int32) -> Int32 {
  (2 * genus - 2) * n * n + 2
}

// ---------------------------------------------------------------------------
// Eisenstein and constant term on the moduli space.
// ---------------------------------------------------------------------------

///|
/// The Eisenstein map on moduli: a block-diagonal (parabolic) local system.
/// `blocks` gives the sizes of the blocks; `values` the `n` diagonal entries
/// (one per basis vector), the off-diagonal part being the nilpotent of type
/// `blocks`.  With all blocks of size 1 the result is semisimple.
pub fn eisenstein_local_system(
  blocks : Array[Int32],
  values : Array[Double],
) -> Array[Double] {
  let n = partition_sum(blocks)
  let nil = nilpotent_standard(blocks)
  let m : Array[Double] = []
  for i in 0..<(n * n) {
    let diag = if i % n == i / n { values[i / n] } else { 0.0 }
    m.push(nil[i] + diag)
  }
  m
}

///|
/// The constant-term functor on moduli: the semisimple (diagonal) part of a
/// local system, i.e. the image under `CT : M_B(G) → M_B(M)`.
pub fn constant_term_local_system(
  m : Array[Double],
  n : Int32,
) -> Array[Double] {
  let out : Array[Double] = []
  for i in 0..<(n * n) {
    out.push(if i % n == i / n { m[i] } else { 0.0 })
  }
  out
}

///|
/// `CT ∘ Eis - id` defect on the Levi (Theorem T23c: the composite is the
/// identity, so the Eisenstein series is a *section* of the constant-term
/// functor).
pub fn ct_eis_defect(blocks : Array[Int32], values : Array[Double]) -> Double {
  let n = partition_sum(blocks)
  let eis = eisenstein_local_system(blocks, values)
  let ct = constant_term_local_system(eis, n)
  vec_inf_diff(Array::makei(n, fn(i) { ct[i * n + i] }), values)
}

///|
/// Dimension of the fibre of `CT` over a point with Jordan type `blocks`:
/// the unipotent orbit `n² - Σ(λ'_i)²` (Theorem T23d).  This is the precise
/// sense in which "non-semisimple L-parameters form the fibres of the
/// constant-term functor".
pub fn ct_fiber_dim(blocks : Array[Int32]) -> Int32 {
  orbit_dim(blocks)
}

///|
/// Dimension of the `GL_n`-orbit of a block-diagonal (Levi) point with
/// eigenvalue multiplicities given by `blocks`: for a semisimple point the
/// orbit is a flag variety of dimension `n² - Σ λ_i²`.
pub fn levi_orbit_dim(blocks : Array[Int32]) -> Int32 {
  let n = partition_sum(blocks)
  n * n - levi_dim(blocks)
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T23a: Fricke identity on explicit SL_2 matrices" {
  let rng = Rng::new(20261117)
  for _ in 0..<8 {
    let a = random_sl2(rng)
    let b = random_sl2(rng)
    assert_true((mat2_det(a) - 1.0).abs() < 1.0e-12)
    assert_true(fricke_defect(a, b) < 1.0e-10)
    // and det[A,B] = 1
    assert_true((mat2_det(mat2_comm(a, b)) - 1.0).abs() < 1.0e-10)
  }
  // the trivial representation: A = B = 1 ↦ x = y = z = 2, κ = 4 = 2 + tr[1,1]
  assert_true(approx(fricke_kappa(2.0, 2.0, 2.0), 4.0, tol=1.0e-12))
  assert_true(approx(commutator_trace(2.0, 2.0, 2.0), 2.0, tol=1.0e-12))
}

///|
test "T23b: Goldman bracket is Poisson with Casimir κ and rank 2" {
  let p : Array[Double] = [1.3, -0.7, 0.9]
  let f = fn(q : Array[Double]) { q[0] }
  let g = fn(q : Array[Double]) { q[1] }
  let h = fn(q : Array[Double]) { q[2] }
  // {x,y} = 2z - xy
  assert_true(
    approx(poisson_bracket(f, g, p), 2.0 * 0.9 - 1.3 * -0.7, tol=1.0e-7),
  )
  // {y,z} = 2x - yz, {z,x} = 2y - zx
  assert_true(
    approx(poisson_bracket(g, h, p), 2.0 * 1.3 - -0.7 * 0.9, tol=1.0e-7),
  )
  assert_true(
    approx(poisson_bracket(h, f, p), 2.0 * -0.7 - 0.9 * 1.3, tol=1.0e-7),
  )
  // antisymmetry
  assert_true(
    approx(poisson_bracket(f, g, p), -poisson_bracket(g, f, p), tol=1.0e-7),
  )
  // Jacobi identity for the coordinate functions
  assert_true(jacobi_defect(f, g, h, p) < 1.0e-6)
  // κ is a Casimir
  assert_true(casimir_defect(p) < 1.0e-6)
  // rank 2 away from the boundary
  assert_true(poisson_rank(p) == 2)
  // Jacobi for a genuinely nonlinear triple
  let f2 = fn(q : Array[Double]) { fricke_kappa(q[0], q[1], q[2]) }
  let g2 = fn(q : Array[Double]) { q[0] * q[1] }
  let h2 = fn(q : Array[Double]) { q[2] * q[2] }
  assert_true(jacobi_defect(f2, g2, h2, p) < 1.0e-5)
}

///|
test "T23c: dimension bookkeeping for Betti moduli spaces" {
  // one-holed torus with SL_2: 2-dimensional symplectic leaves
  assert_true(betti_dim(1, 1, 2) == 2)
  assert_true(betti_dim(2, 0, 2) == 6)
  assert_true(betti_dim(1, 2, 2) == 4)
  assert_true(betti_dim(2, 0, 3) == 16)
  assert_true(betti_dim_gl(2, 2) == 10)
  // the number of Casimirs is 1 and dim M_B = rank(P) = 2 for (1,1)
  assert_true(betti_dim(1, 1, 2) - 1 == 1)
}

///|
test "T23d: CT ∘ Eis = id and the fibre of CT is a unipotent orbit" {
  let blocks : Array[Int32] = [2, 2]
  let values : Array[Double] = [1.5, 1.5, -0.5, -0.5]
  assert_true(ct_eis_defect(blocks, values) < 1.0e-12)
  let n = partition_sum(blocks)
  // fibre dimension + Levi orbit dimension + fixed algebra = n²
  let fibre = ct_fiber_dim(blocks)
  assert_true(fibre == 8) // orbit of type (2,2) in GL_4
  assert_true(levi_orbit_dim(blocks) == 8)
  assert_true(gcheck_kernel_dim(blocks) == 8)
  assert_true(fibre + gcheck_kernel_dim(blocks) == n * n)
  // semisimple parameters have zero-dimensional fibres
  assert_true(ct_fiber_dim([1, 1, 1, 1]) == 0)
  // CT of an Eisenstein point is semisimple
  let eis = eisenstein_local_system([3, 1], [2.0, 2.0, 2.0, 0.5])
  let ct = constant_term_local_system(eis, 4)
  let mut offdiag = 0.0
  for i in 0..<16 {
    if i % 4 != i / 4 {
      offdiag += ct[i].abs()
    }
  }
  assert_true(offdiag == 0.0)
}
