// ============================================================================
// 02_supermanifold.mbt -- M^{m|n} and its super-volume.
//
// ABSTRACTION.  A superchart (x^1..x^m | θ^1..θ^n) is a pair
//   (even body points in R^m, odd fibre Λ(n)).
// We take the even body to be the Poincare upper half-plane H^2 (constant
// curvature -1), so "Poincare-Riemann supermanifold" is a concrete instance
// rather than a slogan.  The super-metric is block diagonal
//   G = diag(g_even, J_odd),  J_odd symplectic,
// and its Berezinian Ber(G) = det(g_even)/det(J_odd) gives the super volume.
// ============================================================================

///|
pub(all) struct SuperManifold {
  dim : GradedDim
  points : Array[Array[Double]] // body coordinates [x, y], y > 0
  odd_form : Array[Double] // n x n antisymmetric matrix
}

///|
/// Standard symplectic form on an even-dimensional odd fibre.
pub fn symplectic_form(n : Int32) -> Array[Double] {
  Array::makei(n * n, fn(k) {
    let i = k / n
    let j = k % n
    if i % 2 == 0 && j == i + 1 {
      1.0
    } else if i % 2 == 1 && j == i - 1 {
      -1.0
    } else {
      0.0
    }
  })
}

///|
pub fn SuperManifold::poincare_halfplane(
  points : Array[Array[Double]],
  n : Int32,
) -> SuperManifold {
  SuperManifold::{
    dim: GradedDim::{ m: 2, n, },
    points,
    odd_form: symplectic_form(n),
  }
}

///|
pub fn SuperManifold::sdim(self : SuperManifold) -> Int32 {
  self.dim.super_dimension()
}

///|
pub fn SuperManifold::fibre_dim(self : SuperManifold) -> Int32 {
  pow2(self.dim.n)
}

///|
/// Riemannian metric of H^2 at p = (x,y): g = y^{-2} I_2.
pub fn poincare_metric(p : Array[Double]) -> Array[Double] {
  let s = 1.0 / (p[1] * p[1])
  [s, 0.0, 0.0, s]
}

///|
/// Gaussian curvature of H^2 is the constant -1.
pub fn poincare_curvature() -> Double {
  -1.0
}

///|
/// Hyperbolic distance: d(p,q) = arcosh(1 + |p-q|^2 / (2 y_p y_q)).
pub fn poincare_distance(p : Array[Double], q : Array[Double]) -> Double {
  let dx = p[0] - q[0]
  let dy = p[1] - q[1]
  let t = 1.0 + (dx * dx + dy * dy) / (2.0 * p[1] * q[1])
  @math.acosh(t)
}

///|
/// sqrt(det g) = 1/y^2: the Riemannian volume density.
pub fn poincare_volume_density(p : Array[Double]) -> Double {
  1.0 / (p[1] * p[1])
}

///|
/// sqrt(Ber G) with G = diag(g_even, J_odd).
/// Because det(J_odd) = 1 for the standard symplectic form, this reduces to
/// the ordinary Riemannian volume density (Theorem T5).
pub fn SuperManifold::super_volume_density(
  self : SuperManifold,
  i : Int32,
) -> Double {
  let g = poincare_metric(self.points[i])
  let det_g = determinant(g, 2)
  let det_j = determinant(self.odd_form, self.dim.n)
  (det_g / det_j).sqrt()
}

///|
pub fn SuperManifold::distance_matrix(self : SuperManifold) -> Array[Double] {
  let n = self.points.length()
  Array::makei(n * n, fn(k) {
    poincare_distance(self.points[k / n], self.points[k % n])
  })
}

///|
/// Berezinian of a pure block-diagonal super-Jacobian (A | B):
/// Ber = det(A_even) / det(B_odd).
pub fn super_jacobian_ber(
  a : Array[Double],
  m : Int32,
  b : Array[Double],
  n : Int32,
) -> Double {
  determinant(a, m) / determinant(b, n)
}

///|
/// Scaling law: under x -> λx, θ -> μθ, Ber(G) scales as λ^{2m}/μ^{2n}.
/// Invariance of the super volume forces μ = λ^{m/n} (Theorem T5b).
pub fn ber_scaling_law(
  m : Int32,
  n : Int32,
  lambda : Double,
  mu : Double,
) -> Double {
  @math.pow(lambda, (2 * m).to_double()) / @math.pow(mu, (2 * n).to_double())
}

///|
/// Computational proof of the odd change-of-variables rule:
/// ∫ dθ_n...dθ_1 f(Bθ) = det(B) ∫ dθ_n...dθ_1 f(θ).
/// Take f = θ_1...θ_n; substitute θ_i -> Σ_j B_ij θ_j in Λ(n) and read the
/// top coefficient by Berezin integration.  It equals det(B).
pub fn berezin_linear_substitution(b : Array[Double], n : Int32) -> Double {
  let mut acc = SuperNumber::one(n)
  for i in 0..<n {
    let mut image = SuperNumber::zero(n)
    for j in 0..<n {
      image = image.add(SuperNumber::theta(n, j).scale(b[i * n + j]))
    }
    acc = acc.mul(image)
  }
  acc.berezin()
}

///|
test "T5: Poincare metric, distance symmetries and curvature" {
  let p = [0.0, 1.0]
  let q = [1.0, 2.0]
  let r = [-0.5, 3.0]
  assert_true(approx(poincare_distance(p, p), 0.0, tol=1.0e-12))
  assert_true(approx(poincare_distance(p, q), poincare_distance(q, p)))
  // triangle inequality
  assert_true(
    poincare_distance(p, r) <=
    poincare_distance(p, q) + poincare_distance(q, r) + 1.0e-9,
  )
  // scaling (x,y) -> (λx, λy) is an isometry of H^2
  let lam = 3.0
  let pp = [lam * p[0], lam * p[1]]
  let qq = [lam * q[0], lam * q[1]]
  assert_true(approx(poincare_distance(pp, qq), poincare_distance(p, q)))
  assert_true(approx(poincare_curvature(), -1.0))
}

///|
test "T5: super volume reduces to Riemannian volume; scaling law" {
  let sm = SuperManifold::poincare_halfplane([[0.0, 2.0], [1.0, 4.0]], 2)
  assert_eq(sm.sdim(), 0)
  assert_eq(sm.fibre_dim(), 4)
  assert_true(
    approx(
      sm.super_volume_density(0),
      poincare_volume_density(sm.points[0]),
      tol=1.0e-12,
    ),
  )
  // μ = λ^{m/n} keeps Ber(G) invariant
  assert_true(approx(ber_scaling_law(2, 2, 2.0, 2.0), 1.0))
  assert_true(approx(ber_scaling_law(2, 4, 3.0, 3.0.sqrt()), 1.0, tol=1.0e-9))
  let dm = sm.distance_matrix()
  assert_true(approx(dm[1], dm[2]))
  assert_true(approx(dm[0], 0.0))
}

///|
test "T5: Berezinian = det ratio, and odd substitution gives det(B)" {
  let a = [2.0, 0.0, 0.0, 3.0]
  let b = [2.0, 1.0, 1.0, 3.0]
  assert_true(approx(super_jacobian_ber(a, 2, b, 2), 6.0 / 5.0))
  assert_true(
    approx(berezin_linear_substitution(b, 2), determinant(b, 2), tol=1.0e-12),
  )
  let b3 = [1.0, 2.0, 0.0, 3.0, 1.0, 1.0, 0.0, 2.0, 1.0]
  assert_true(
    approx(berezin_linear_substitution(b3, 3), determinant(b3, 3), tol=1.0e-12),
  )
}
