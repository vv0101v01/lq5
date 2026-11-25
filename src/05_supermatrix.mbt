// ============================================================================
// 05_supermatrix.mbt -- Noncommutative geometry and super-matrix equations.
//
// ABSTRACTION.  Mat(m|n) is the Z_2-graded associative algebra of endomorphisms
// of the super vector space R^{m|n}.  A matrix is stored in row-major order
// with even rows/cols first; its parity is a declared grading label and the
// product sets |AB| = |A| + |B| (mod 2).  The graded commutator
//   [A,B] = AB - (-1)^{|A||B|} BA
// turns Mat(m|n) into a Lie superalgebra (Theorem T8).
//
// Noncommutative geometry enters through a finite SPECTRAL TRIPLE (A,H,D):
// A = C^2 (functions on two "agents"), D = [[0,t],[t,0]].  Connes' distance
// formula d(0,1) = sup{|f0-f1| : ||[D,f]|| <= 1} = 1/|t| (Theorem T9):
// geometry is the inverse of communication strength.
// ============================================================================

///|
pub(all) struct SuperMatrix {
  re : Int32 // even rows
  ro : Int32 // odd rows
  ce : Int32 // even cols
  co : Int32 // odd cols
  parity : Int32 // 0 even, 1 odd (grading label)
  data : Array[Double] // row-major, (re+ro)*(ce+co)
}

///|
pub fn SuperMatrix::new(
  re : Int32,
  ro : Int32,
  ce : Int32,
  co : Int32,
  parity : Int32,
  data : Array[Double],
) -> SuperMatrix {
  SuperMatrix::{ re, ro, ce, co, parity, data, }
}

///|
pub fn SuperMatrix::rows(self : SuperMatrix) -> Int32 {
  self.re + self.ro
}

///|
pub fn SuperMatrix::cols(self : SuperMatrix) -> Int32 {
  self.ce + self.co
}

///|
pub fn SuperMatrix::get(self : SuperMatrix, i : Int32, j : Int32) -> Double {
  self.data[i * self.cols() + j]
}

///|
/// Extract one of the four blocks: rb/cb = 0 (even) or 1 (odd).
pub fn SuperMatrix::block(
  self : SuperMatrix,
  rb : Int32,
  cb : Int32,
) -> Array[Double] {
  let ro0 = if rb == 0 { 0 } else { self.re }
  let rn = if rb == 0 { self.re } else { self.ro }
  let co0 = if cb == 0 { 0 } else { self.ce }
  let cn = if cb == 0 { self.ce } else { self.co }
  Array::makei(rn * cn, fn(t) {
    self.data[(ro0 + t / cn) * self.cols() + co0 + t % cn]
  })
}

///|
pub fn SuperMatrix::add(self : SuperMatrix, other : SuperMatrix) -> SuperMatrix {
  SuperMatrix::{
    re: self.re,
    ro: self.ro,
    ce: self.ce,
    co: self.co,
    parity: self.parity,
    data: Array::makei(self.data.length(), fn(k) {
      self.data[k] + other.data[k]
    }),
  }
}

///|
pub fn SuperMatrix::sub(self : SuperMatrix, other : SuperMatrix) -> SuperMatrix {
  self.add(other.scale(-1.0))
}

///|
pub fn SuperMatrix::scale(self : SuperMatrix, s : Double) -> SuperMatrix {
  SuperMatrix::{
    re: self.re,
    ro: self.ro,
    ce: self.ce,
    co: self.co,
    parity: self.parity,
    data: Array::makei(self.data.length(), fn(k) { s * self.data[k] }),
  }
}

///|
pub fn SuperMatrix::mul(self : SuperMatrix, other : SuperMatrix) -> SuperMatrix {
  let r = self.rows()
  let k = self.cols()
  let c = other.cols()
  let data = Array::makei(r * c, fn(t) {
    let i = t / c
    let j = t % c
    let mut s = 0.0
    for h in 0..<k {
      s += self.data[i * k + h] * other.data[h * c + j]
    }
    s
  })
  SuperMatrix::{
    re: self.re,
    ro: self.ro,
    ce: other.ce,
    co: other.co,
    parity: (self.parity + other.parity) % 2,
    data,
  }
}

///|
pub fn SuperMatrix::norm(self : SuperMatrix) -> Double {
  vec_inf_norm(self.data)
}

///|
/// Supertrace str(M) = tr(M_EE) - tr(M_OO) for a square m|n matrix.
pub fn SuperMatrix::supertrace(self : SuperMatrix) -> Double {
  let n = self.rows()
  let mut s = 0.0
  for i in 0..<n {
    let sign = if i < self.re { 1.0 } else { -1.0 }
    s += sign * self.data[i * n + i]
  }
  s
}

///|
/// Graded transpose with sign flip on the odd-odd block:
/// (M^st)_{ij} = (-1)^{p(i)p(j)} M_{ji}.
pub fn SuperMatrix::supertranspose(self : SuperMatrix) -> SuperMatrix {
  let r = self.rows()
  let c = self.cols()
  let data = Array::makei(r * c, fn(t) {
    let i = t / c
    let j = t % c
    let pi = if i < self.re { 0 } else { 1 }
    let pj = if j < self.ce { 0 } else { 1 }
    let s = if pi * pj % 2 == 0 { 1.0 } else { -1.0 }
    s * self.data[j * c + i]
  })
  SuperMatrix::{
    re: self.ce,
    ro: self.co,
    ce: self.re,
    co: self.ro,
    parity: self.parity,
    data,
  }
}

///|
/// Supercommutator [A,B] = AB - (-1)^{|A||B|} BA.
pub fn SuperMatrix::bracket(
  self : SuperMatrix,
  other : SuperMatrix,
) -> SuperMatrix {
  let s = if self.parity * other.parity % 2 == 0 { 1.0 } else { -1.0 }
  self.mul(other).sub(other.mul(self).scale(s))
}

///|
/// Super-Jacobiator (Theorem T8).
pub fn super_jacobiator_m(
  a : SuperMatrix,
  b : SuperMatrix,
  c : SuperMatrix,
) -> SuperMatrix {
  let t1 = a
    .bracket(b.bracket(c))
    .scale(if a.parity * c.parity % 2 == 0 { 1.0 } else { -1.0 })
  let t2 = b
    .bracket(c.bracket(a))
    .scale(if b.parity * a.parity % 2 == 0 { 1.0 } else { -1.0 })
  let t3 = c
    .bracket(a.bracket(b))
    .scale(if c.parity * b.parity % 2 == 0 { 1.0 } else { -1.0 })
  t1.add(t2).add(t3)
}

///|
/// [a,[b,c]] - [[a,b],c] - (-1)^{|a||b|}[b,[a,c]] = 0 (graded derivation).
pub fn derivation_defect_m(
  a : SuperMatrix,
  b : SuperMatrix,
  c : SuperMatrix,
) -> SuperMatrix {
  let s = if a.parity * b.parity % 2 == 0 { 1.0 } else { -1.0 }
  a
  .bracket(b.bracket(c))
  .sub(a.bracket(b).bracket(c))
  .sub(b.bracket(a.bracket(c)).scale(s))
}

///|
/// Berezinian of an even square supermatrix M = [[A,B],[C,D]]:
/// Ber(M) = det(A - B D^{-1} C) / det(D).
pub fn SuperMatrix::berezinian(self : SuperMatrix) -> Double {
  let m = self.re
  let n = self.ro
  let a = self.block(0, 0)
  if n == 0 {
    return determinant(a, m)
  }
  let b = self.block(0, 1)
  let c = self.block(1, 0)
  let d = self.block(1, 1)
  let d_inv = inverse(d, n)
  let bd = matrix_product(b, m, n, d_inv, n)
  let bdc = matrix_product(bd, m, n, c, m)
  let schur = Array::makei(m * m, fn(k) { a[k] - bdc[k] })
  determinant(schur, m) / determinant(d, n)
}

// ---------------------------------------------------------------------------
// Graded Sylvester equation:  A X + σ X B = C,  σ = (-1)^{|A||X|}.
// Vectorised: (I ⊗ A + σ Bᵀ ⊗ I) vec(X) = vec(C).  Solvable iff
// λ_i(A) + σ λ_j(B) ≠ 0 for all i,j (Theorem T10).
// ---------------------------------------------------------------------------

///|
pub fn sylvester_sigma(a : SuperMatrix, x_parity : Int32) -> Double {
  if a.parity * x_parity % 2 == 0 {
    1.0
  } else {
    -1.0
  }
}

///|
pub fn graded_sylvester_solve(
  a : SuperMatrix,
  b : SuperMatrix,
  c : SuperMatrix,
  x_parity : Int32,
) -> SuperMatrix {
  let r = a.rows()
  let cc = b.cols()
  let n = r * cc
  let sigma = sylvester_sigma(a, x_parity)
  let k = Array::make(n * n, 0.0)
  let rhs = Array::make(n, 0.0)
  for i in 0..<r {
    for j in 0..<cc {
      let row = i * cc + j
      for h in 0..<r {
        k[row * n + h * cc + j] += a.get(i, h)
      }
      for l in 0..<cc {
        k[row * n + i * cc + l] += sigma * b.get(l, j)
      }
      rhs[row] = c.get(i, j)
    }
  }
  let x = gauss_solve(k, rhs, n)
  SuperMatrix::{
    re: a.re,
    ro: a.ro,
    ce: b.ce,
    co: b.co,
    parity: x_parity,
    data: x,
  }
}

///|
pub fn graded_sylvester_residual(
  a : SuperMatrix,
  b : SuperMatrix,
  c : SuperMatrix,
  x : SuperMatrix,
) -> Double {
  let sigma = sylvester_sigma(a, x.parity)
  a.mul(x).add(x.mul(b).scale(sigma)).sub(c).norm()
}

// ---------------------------------------------------------------------------
// Finite spectral triple (Connes)： two-point noncommutative space.
// ---------------------------------------------------------------------------

///|
pub(all) struct SpectralTriple {
  coupling : Double // t = |D_01|
}

///|
pub fn SpectralTriple::new(coupling : Double) -> SpectralTriple {
  SpectralTriple::{ coupling, }
}

///|
/// ||[D, diag(f0,f1)]|| = |f0 - f1| |t| : the Lipschitz seminorm.
pub fn SpectralTriple::commutator_norm(
  self : SpectralTriple,
  f : Array[Double],
) -> Double {
  (f[0] - f[1]).abs() * self.coupling.abs()
}

///|
/// d_Connes(0,1) = 1/|t| (Theorem T9).
pub fn SpectralTriple::connes_distance(self : SpectralTriple) -> Double {
  1.0 / self.coupling.abs()
}

///|
/// Numerically maximise |f0-f1| under the Lipschitz constraint (Algorithm A6).
pub fn SpectralTriple::connes_distance_grid(
  self : SpectralTriple,
  steps : Int32,
) -> Double {
  let bound = 2.0 * self.connes_distance()
  let mut best = 0.0
  for j in 0..<=steps {
    let d = bound * j.to_double() / steps.to_double()
    let f = [d, 0.0]
    if self.commutator_norm(f) <= 1.0 + 1.0e-12 && d > best {
      best = d
    }
  }
  best
}

///|
/// Emergent metric radius of an agent graph with coupling K.
pub fn emergent_metric_radius(coupling : Double) -> Double {
  1.0 / coupling.abs()
}

///|
test "T8: supertrace cyclicity, supertranspose involution, Berezinian" {
  // even matrix: EE + OO blocks;  odd matrix: EO + OE blocks (2|1)
  let e = SuperMatrix::new(2, 1, 2, 1, 0, [
    2.0, 0.5, 0.0, 0.25, 3.0, 0.0, 0.0, 0.0, 5.0,
  ])
  let o = SuperMatrix::new(2, 1, 2, 1, 1, [
    0.0, 0.0, 1.5, 0.0, 0.0, -0.75, 0.4, 1.2, 0.0,
  ])
  let eo = e.mul(o)
  let oe = o.mul(e)
  // str(EO) = (-1)^{|E||O|} str(OE);  both sides vanish here by block parity
  assert_true(eo.supertrace().abs() < 1.0e-12)
  assert_true(oe.supertrace().abs() < 1.0e-12)
  assert_true(
    e.supertranspose().supertranspose().norm() < 1.0e-12 ||
    vec_inf_diff(e.supertranspose().supertranspose().data, e.data) < 1.0e-12,
  )
  // Ber(diag(A,D)) = det(A)/det(D)
  assert_true(approx(e.berezinian(), 1.175, tol=1.0e-12))
  let with_odd = SuperMatrix::new(2, 1, 2, 1, 0, [
    2.0, 0.0, 1.0, 0.0, 3.0, 0.5, 0.25, 0.75, 4.0,
  ])
  // independent closed form: Schur complement with D = [4]
  let a = with_odd.block(0, 0)
  let bb = with_odd.block(0, 1)
  let c = with_odd.block(1, 0)
  let schur = Array::makei(4, fn(k) { a[k] - bb[k / 2] * c[k % 2] / 4.0 })
  assert_true(
    approx(with_odd.berezinian(), determinant(schur, 2) / 4.0, tol=1.0e-12),
  )
}

///|
test "T8: Super-Jacobi Identity for Mat(2|1)" {
  let x = SuperMatrix::new(2, 1, 2, 1, 0, [
    0.7, 0.2, 0.1, -0.3, 1.1, 0.4, 0.2, -0.5, 0.9,
  ])
  let y = SuperMatrix::new(2, 1, 2, 1, 1, [
    0.3, -0.8, 0.6, 1.2, 0.1, -0.4, 0.7, 0.3, 0.5,
  ])
  let z = SuperMatrix::new(2, 1, 2, 1, 0, [
    1.4, 0.5, -0.2, 0.1, -0.9, 0.3, -0.6, 0.8, 1.3,
  ])
  let z_odd = SuperMatrix::new(2, 1, 2, 1, 1, [
    0.2, 0.9, 1.1, -0.4, 0.3, 0.8, 0.6, -1.0, 0.2,
  ])
  assert_true(super_jacobiator_m(x, y, z).norm() < 1.0e-11)
  assert_true(super_jacobiator_m(x, y, z_odd).norm() < 1.0e-11)
  assert_true(derivation_defect_m(x, y, z).norm() < 1.0e-11)
  assert_true(derivation_defect_m(x, y, z_odd).norm() < 1.0e-11)
  // odd * odd = even, and the bracket of two odds is the ANTICOMMUTATOR
  let b = y.bracket(z_odd)
  assert_eq(b.parity, 0)
  let anti = y.mul(z_odd).add(z_odd.mul(y))
  assert_true(vec_inf_diff(b.data, anti.data) < 1.0e-12)
}

///|
test "T10: graded Sylvester equation solves A X + σ X B = C" {
  let a = SuperMatrix::new(3, 0, 3, 0, 0, [
    5.0, 0.2, 0.0, 0.1, 7.0, 0.3, 0.0, 0.4, 9.0,
  ])
  let b = SuperMatrix::new(3, 0, 3, 0, 0, [
    1.0, 0.1, 0.0, 0.0, 2.0, 0.2, 0.1, 0.0, 3.0,
  ])
  let c = SuperMatrix::new(3, 0, 3, 0, 0, [
    1.0, 2.0, 0.5, -0.5, 0.3, 1.2, 0.7, 0.2, -1.0,
  ])
  let x0 = graded_sylvester_solve(a, b, c, 0)
  assert_true(graded_sylvester_residual(a, b, c, x0) < 1.0e-10)
  let x1 = graded_sylvester_solve(a, b, c, 1)
  assert_true(sylvester_sigma(a, 1) == 1.0) // |a| = 0 so σ = 1 always here
  assert_true(graded_sylvester_residual(a, b, c, x1) < 1.0e-10)
  // genuinely odd A flips σ to -1
  let a_odd = SuperMatrix::new(3, 0, 3, 0, 1, a.data)
  let x_odd = graded_sylvester_solve(a_odd, b, c, 1)
  assert_true(sylvester_sigma(a_odd, 1) == -1.0)
  assert_true(graded_sylvester_residual(a_odd, b, c, x_odd) < 1.0e-10)
}

///|
test "T9: Connes distance = inverse coupling (emergent geometry)" {
  let st = SpectralTriple::new(2.5)
  assert_true(approx(st.connes_distance(), 0.4))
  let g = st.connes_distance_grid(4000)
  assert_true(approx(g, 0.4, tol=1.0e-3))
  assert_true(approx(emergent_metric_radius(2.5), st.connes_distance()))
  let f = [0.3, 0.0]
  assert_true(approx(st.commutator_norm(f), 0.75))
}
