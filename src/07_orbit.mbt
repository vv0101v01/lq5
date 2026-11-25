// ============================================================================
// 07_orbit.mbt -- Second-order super-mass tensor and its orbit space.
//
// ABSTRACTION.  The "mass tensor" of an agent field is the second moment
//   G_even = ⟨x xᵀ⟩,  x = (cos θ, sin θ) ∈ S^1 ⊂ R²,
// together with the odd (symplectic) block
//   G_odd  = ⟨ x_i ∧ x_{i+1} ⟩ J  (lattice chiral current).
// G = diag(G_even, G_odd) is a graded symmetric/symplectic supermetric.
//
// ORBIT SPACE.  Under SO(2) conjugation (relabelling of the even frame),
// G_even ↦ Rᵀ G_even R.  The quotient is R²/S₂ (unordered eigenvalue pairs),
// a Weyl chamber.  Theorem T12: invariants are (tr, det, Ber); the quotient
// metric is the optimal-matching distance between SORTED eigenvalues; and
// the anisotropy (λ₁-λ₂)/(λ₁+λ₂) equals |⟨exp(2iθ)⟩|, i.e. it detects
// NEMATIC (anti-phase / two-cluster) order that the polar order parameter
// R = |⟨exp(iθ)⟩| cannot see.  That is an emergent observable.
// ============================================================================

///|
pub(all) struct SuperMassTensor {
  even : Array[Double] // 2x2 symmetric positive semidefinite
  odd : Array[Double] // 2x2 antisymmetric
}

///|
/// Build G from the instantaneous agent phases.
pub fn SuperMassTensor::from_phases(theta : Array[Double]) -> SuperMassTensor {
  let n = theta.length()
  let mut xx = 0.0
  let mut xy = 0.0
  let mut yy = 0.0
  for x in theta {
    let c = @math.cos(x)
    let s = @math.sin(x)
    xx += c * c
    xy += c * s
    yy += s * s
  }
  let nf = n.to_double()
  xx /= nf
  xy /= nf
  yy /= nf
  // chiral current: ⟨sin(θ_{i+1}-θ_i)⟩ on the ring
  let mut w = 0.0
  for i in 0..<n {
    w += @math.sin(theta[(i + 1) % n] - theta[i])
  }
  w /= nf
  SuperMassTensor::{ even: [xx, xy, xy, yy], odd: [0.0, w, -w, 0.0], }
}

///|
/// Eigenvalues of a symmetric 2x2 matrix, sorted descending.
pub fn eigen2(g : Array[Double]) -> Array[Double] {
  let tr = g[0] + g[3]
  let det = g[0] * g[3] - g[1] * g[2]
  let disc = (tr * tr / 4.0 - det).sqrt()
  let l1 = tr / 2.0 + disc
  let l2 = tr / 2.0 - disc
  if l1 >= l2 {
    [l1, l2]
  } else {
    [l2, l1]
  }
}

///|
pub fn SuperMassTensor::trace(self : SuperMassTensor) -> Double {
  self.even[0] + self.even[3]
}

///|
pub fn SuperMassTensor::det_even(self : SuperMassTensor) -> Double {
  determinant(self.even, 2)
}

///|
pub fn SuperMassTensor::det_odd(self : SuperMassTensor) -> Double {
  determinant(self.odd, 2)
}

///|
/// Ber(G) = det(G_even)/det(G_odd) (regularised when the chiral current is 0).
pub fn SuperMassTensor::berezinian(
  self : SuperMassTensor,
  eps~ : Double = 1.0e-8,
) -> Double {
  let d = self.det_odd()
  self.det_even() / (if d.abs() < eps { eps } else { d })
}

///|
pub fn SuperMassTensor::canonical_form(self : SuperMassTensor) -> Array[Double] {
  eigen2(self.even)
}

///|
pub fn SuperMassTensor::anisotropy(self : SuperMassTensor) -> Double {
  let l = eigen2(self.even)
  (l[0] - l[1]) / (l[0] + l[1] + 1.0e-12)
}

///|
/// Full invariant vector of the orbit: [tr, det_even, Ber, λ₁, λ₂, anisotropy, w].
pub fn SuperMassTensor::invariants(self : SuperMassTensor) -> Array[Double] {
  let l = eigen2(self.even)
  [
    self.trace(),
    self.det_even(),
    self.berezinian(),
    l[0],
    l[1],
    self.anisotropy(),
    self.odd[1],
  ]
}

///|
/// Conjugation by a rotation: G ↦ Rᵀ G R.
pub fn SuperMassTensor::conjugate(
  self : SuperMassTensor,
  r : Array[Double],
) -> SuperMassTensor {
  SuperMassTensor::{ even: pullback_metric(self.even, r, 2), odd: self.odd, }
}

///|
/// Quotient metric on R²/S₂: optimal matching of sorted eigenvalues.
pub fn orbit_distance(a : SuperMassTensor, b : SuperMassTensor) -> Double {
  let la = eigen2(a.even)
  let lb = eigen2(b.even)
  ((la[0] - lb[0]) * (la[0] - lb[0]) + (la[1] - lb[1]) * (la[1] - lb[1])).sqrt()
}

///|
/// Nematic (second harmonic) order parameter `r_2`: the head-tail symmetric
/// coherence that a rank-2 tensor — not a vector — can detect.
pub fn nematic_order(theta : Array[Double]) -> Double {
  order_harmonic(theta, 2)
}

///|
test "T12: orbit invariants are SO(2)-invariant; quotient metric is optimal matching" {
  let theta = [0.2, 0.9, 1.7, 2.4, 3.9, 5.5]
  let g = SuperMassTensor::from_phases(theta)
  let r = rotation(0.83)
  let gr = g.conjugate(r)
  assert_true(approx(orbit_distance(g, gr), 0.0, tol=1.0e-10))
  assert_true(approx(g.trace(), gr.trace(), tol=1.0e-10))
  assert_true(approx(g.det_even(), gr.det_even(), tol=1.0e-10))
  let other = SuperMassTensor::from_phases([0.0, 0.0, 0.0, 0.0])
  assert_true(orbit_distance(g, other) > 0.0)
  // symmetry and positivity of the distance
  assert_true(approx(orbit_distance(g, other), orbit_distance(other, g)))
  assert_true(approx(orbit_distance(g, g), 0.0))
}

///|
test "T12: anisotropy = |<exp(2iθ)>| detects nematic order invisible to R" {
  // uniform on the circle: no polar, no nematic order
  let uniform = [0.0, 1.5707963267948966, 3.141592653589793, 4.71238898038469]
  let gu = SuperMassTensor::from_phases(uniform)
  assert_true(approx(order_parameter(uniform), 0.0, tol=1.0e-12))
  assert_true(approx(gu.anisotropy(), 0.0, tol=1.0e-12))
  // two anti-phase clusters: R = 0 but nematic order = 1
  let clusters = [0.0, 0.0, 3.141592653589793, 3.141592653589793]
  let gc = SuperMassTensor::from_phases(clusters)
  assert_true(approx(order_parameter(clusters), 0.0, tol=1.0e-12))
  assert_true(approx(gc.anisotropy(), 1.0, tol=1.0e-10))
  assert_true(approx(gc.anisotropy(), nematic_order(clusters), tol=1.0e-10))
  // complete synchrony: R = 1 and nematic = 1
  let sync = Array::make(6, 0.7)
  let gs = SuperMassTensor::from_phases(sync)
  assert_true(approx(order_parameter(sync), 1.0, tol=1.0e-12))
  assert_true(approx(gs.anisotropy(), 1.0, tol=1.0e-10))
  assert_true(approx(gs.trace(), 1.0, tol=1.0e-10))
}

///|
test "T12: Berezinian is a graded invariant of the orbit" {
  let g = SuperMassTensor::from_phases([0.1, 0.4, 1.2, 2.0, 3.4, 5.9])
  let inv = g.invariants()
  assert_eq(inv.length(), 7)
  assert_true(inv[3] >= inv[4])
  assert_true(inv[0] > 0.0)
  // det_odd = w²
  assert_true(approx(g.det_odd(), g.odd[1] * g.odd[1], tol=1.0e-12))
  let d = g.det_odd()
  let expected = g.det_even() / (if d.abs() < 1.0e-8 { 1.0e-8 } else { d })
  assert_true(approx(g.berezinian(), expected, tol=1.0e-12))
}
