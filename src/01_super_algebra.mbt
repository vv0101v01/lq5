// ============================================================================
// 01_super_algebra.mbt -- Z_2-graded commutative algebra Λ(n).
//
// ABSTRACTION.  A "superpoint" M^{m|n} is a locally ringed space whose sheaf
// of functions is C^∞(R^m) ⊗ Λ(θ_1,...,θ_n).  We discretise the odd fibre:
// a Grassmann monomial θ_{i_1}...θ_{i_k} (indices increasing) is encoded by
// the Int32 bitmask (1<<i_1)|...|(1<<i_k).  A supernumber is then just the
// coefficient vector of length 2^n.  Nilpotency θ_i^2 = 0 becomes the single
// bit test (a & b) == 0, and the Koszul sign rule becomes popcount arithmetic.
//
// This is the "computer abstraction" from which Theorem T1..T4 are derived.
// ============================================================================

///|
pub(all) struct SuperNumber {
  n : Int32 // number of odd generators (n <= 16 in practice)
  coeff : Array[Double] // index = monomial bitmask, length 2^n
}

///|
pub fn SuperNumber::zero(n : Int32) -> SuperNumber {
  SuperNumber::{ n, coeff: Array::make(pow2(n), 0.0), }
}

///|
pub fn SuperNumber::one(n : Int32) -> SuperNumber {
  let x = SuperNumber::zero(n)
  x.coeff[0] = 1.0
  x
}

///|
pub fn SuperNumber::body(n : Int32, x : Double) -> SuperNumber {
  let s = SuperNumber::zero(n)
  s.coeff[0] = x
  s
}

///|
/// The odd generator θ_i (i = 0 .. n-1).
pub fn SuperNumber::theta(n : Int32, i : Int32) -> SuperNumber {
  let s = SuperNumber::zero(n)
  s.coeff[1 << i] = 1.0
  s
}

///|
pub fn SuperNumber::add(self : SuperNumber, other : SuperNumber) -> SuperNumber {
  SuperNumber::{
    n: self.n,
    coeff: Array::makei(pow2(self.n), fn(k) { self.coeff[k] + other.coeff[k] }),
  }
}

///|
pub fn SuperNumber::sub(self : SuperNumber, other : SuperNumber) -> SuperNumber {
  self.add(other.scale(-1.0))
}

///|
pub fn SuperNumber::scale(self : SuperNumber, s : Double) -> SuperNumber {
  SuperNumber::{
    n: self.n,
    coeff: Array::makei(pow2(self.n), fn(k) { s * self.coeff[k] }),
  }
}

///|
/// Koszul sign: θ_A θ_B = (-1)^{#{(i in A, j in B): j < i}} θ_{A∪B},
/// and 0 if A ∩ B is nonempty.  Returns 0.0 for a vanishing product.
pub fn koszul_sign(a : Int32, b : Int32) -> Double {
  if (a & b) != 0 {
    return 0.0
  }
  let mut inversions = 0
  for i in 0..<31 {
    if (a & (1 << i)) != 0 {
      inversions += popcount(b & ((1 << i) - 1))
    }
  }
  if inversions % 2 == 0 {
    1.0
  } else {
    -1.0
  }
}

///|
/// Graded-commutative product in Λ(n).
pub fn SuperNumber::mul(self : SuperNumber, other : SuperNumber) -> SuperNumber {
  let size = pow2(self.n)
  let out = SuperNumber::zero(self.n)
  for a in 0..<size {
    if self.coeff[a] != 0.0 {
      for b in 0..<size {
        if other.coeff[b] != 0.0 {
          let s = koszul_sign(a, b)
          if s != 0.0 {
            out.coeff[a | b] += s * self.coeff[a] * other.coeff[b]
          }
        }
      }
    }
  }
  out
}

///|
/// Left odd derivation ∂/∂θ_i of parity 1 (graded Leibniz rule).
pub fn SuperNumber::odd_derivative(
  self : SuperNumber,
  i : Int32,
) -> SuperNumber {
  let size = pow2(self.n)
  let bit = 1 << i
  let out = SuperNumber::zero(self.n)
  for mask in 0..<size {
    if (mask & bit) != 0 && self.coeff[mask] != 0.0 {
      let before = popcount(mask & (bit - 1))
      let s = if before % 2 == 0 { 1.0 } else { -1.0 }
      out.coeff[mask ^ bit] += s * self.coeff[mask]
    }
  }
  out
}

///|
/// Berezin integral: coefficient of the top monomial θ_1...θ_n.
/// Convention ∫ dθ_n ... dθ_1 θ_1...θ_n = 1.
pub fn SuperNumber::berezin(self : SuperNumber) -> Double {
  self.coeff[pow2(self.n) - 1]
}

///|
pub fn SuperNumber::body_part(self : SuperNumber) -> Double {
  self.coeff[0]
}

///|
pub fn SuperNumber::soul_norm(self : SuperNumber) -> Double {
  let mut s = 0.0
  for k in 1..<pow2(self.n) {
    s += self.coeff[k].abs()
  }
  s
}

///|
/// exp(body + soul) = exp(body) * Σ_{k=0}^{n} soul^k/k!  (nilpotent soul).
pub fn SuperNumber::exp(self : SuperNumber) -> SuperNumber {
  let n = self.n
  let body = self.body_part()
  let soul = SuperNumber::{
    n,
    coeff: Array::makei(pow2(n), fn(k) {
      if k == 0 {
        0.0
      } else {
        self.coeff[k]
      }
    }),
  }
  let mut term = SuperNumber::one(n)
  let mut acc = SuperNumber::one(n)
  for k in 1..<(n + 2) {
    term = term.mul(soul).scale(1.0 / k.to_double())
    acc = acc.add(term)
  }
  acc.scale(@math.exp(body))
}

// ---------------------------------------------------------------------------
// Homogeneous graded elements: the carrier of the graded Lie algebra.
// ---------------------------------------------------------------------------

///|
pub(all) struct GradedElement {
  deg : Int32 // 0 or 1
  num : SuperNumber
}

///|
pub fn GradedElement::theta(n : Int32, i : Int32) -> GradedElement {
  GradedElement::{ deg: 1, num: SuperNumber::theta(n, i), }
}

///|
pub fn GradedElement::even(n : Int32, x : Double) -> GradedElement {
  GradedElement::{ deg: 0, num: SuperNumber::body(n, x), }
}

///|
pub fn GradedElement::add(
  self : GradedElement,
  other : GradedElement,
) -> GradedElement {
  GradedElement::{ deg: self.deg, num: self.num.add(other.num), }
}

///|
pub fn GradedElement::scale(self : GradedElement, s : Double) -> GradedElement {
  GradedElement::{ deg: self.deg, num: self.num.scale(s), }
}

///|
pub fn GradedElement::sub(
  self : GradedElement,
  other : GradedElement,
) -> GradedElement {
  self.add(other.scale(-1.0))
}

///|
pub fn GradedElement::mul(
  self : GradedElement,
  other : GradedElement,
) -> GradedElement {
  GradedElement::{
    deg: (self.deg + other.deg) % 2,
    num: self.num.mul(other.num),
  }
}

///|
/// Supercommutator [a,b] = ab - (-1)^{|a||b|} ba.
pub fn GradedElement::bracket(
  self : GradedElement,
  other : GradedElement,
) -> GradedElement {
  let s = if self.deg * other.deg % 2 == 0 { 1.0 } else { -1.0 }
  GradedElement::{
    deg: (self.deg + other.deg) % 2,
    num: self.num.mul(other.num).sub(other.num.mul(self.num).scale(s)),
  }
}

///|
pub fn GradedElement::norm(self : GradedElement) -> Double {
  vec_inf_norm(self.num.coeff)
}

///|
/// Super-Jacobiator
///   (-1)^{|a||c|}[a,[b,c]] + (-1)^{|b||a|}[b,[c,a]] + (-1)^{|c||b|}[c,[a,b]].
/// Vanishing is Theorem T3.
pub fn super_jacobiator(
  a : GradedElement,
  b : GradedElement,
  c : GradedElement,
) -> GradedElement {
  let t1 = a
    .bracket(b.bracket(c))
    .scale(if a.deg * c.deg % 2 == 0 { 1.0 } else { -1.0 })
  let t2 = b
    .bracket(c.bracket(a))
    .scale(if b.deg * a.deg % 2 == 0 { 1.0 } else { -1.0 })
  let t3 = c
    .bracket(a.bracket(b))
    .scale(if c.deg * b.deg % 2 == 0 { 1.0 } else { -1.0 })
  t1.add(t2).add(t3)
}

///|
/// Defect of ad_a being a graded derivation:
///   [a,[b,c]] - [[a,b],c] - (-1)^{|a||b|}[b,[a,c]].
pub fn derivation_defect(
  a : GradedElement,
  b : GradedElement,
  c : GradedElement,
) -> GradedElement {
  let s = if a.deg * b.deg % 2 == 0 { 1.0 } else { -1.0 }
  a
  .bracket(b.bracket(c))
  .sub(a.bracket(b).bracket(c))
  .sub(b.bracket(a.bracket(c)).scale(s))
}

///|
test "T1: Λ(n) is 2^n-dimensional and θ_i^2 = 0" {
  for n in 1..<5 {
    let x = SuperNumber::theta(n, 0)
    assert_eq(x.coeff.length(), pow2(n))
    assert_eq(x.mul(x).soul_norm(), 0.0)
    assert_eq(x.mul(x).body_part(), 0.0)
  }
}

///|
test "T2: graded commutativity and the Koszul sign" {
  let n = 4
  let a = SuperNumber::theta(n, 0)
  let b = SuperNumber::theta(n, 2)
  let ab = a.mul(b)
  let ba = b.mul(a)
  assert_true(approx(ab.coeff[(1 << 0) | (1 << 2)], 1.0))
  assert_true(approx(ba.coeff[(1 << 0) | (1 << 2)], -1.0))
  assert_true(approx(ab.add(ba).soul_norm(), 0.0))
  assert_eq(koszul_sign(1, 1), 0.0)
}

///|
test "T3: Berezin integral is a total derivative (super Stokes theorem)" {
  let n = 3
  let f = SuperNumber::{
    n,
    coeff: Array::makei(pow2(n), fn(k) { (k + 1).to_double() * 0.125 }),
  }
  let prod = SuperNumber::theta(n, 0)
    .mul(SuperNumber::theta(n, 1))
    .mul(SuperNumber::theta(n, 2))
  assert_true(approx(prod.berezin(), 1.0))
  for i in 0..<n {
    let d = f.odd_derivative(i)
    assert_true(approx(d.berezin(), 0.0))
  }
}

///|
test "T3: super-Jacobi identity and ad as a graded derivation" {
  let n = 3
  let a = GradedElement::theta(n, 0).add(GradedElement::even(n, 0.5))
  let b = GradedElement::theta(n, 1).add(GradedElement::even(n, -0.25))
  let c = GradedElement::theta(n, 2).add(GradedElement::even(n, 1.5))
  // homogeneous pieces: verify all 8 parity combinations
  let odds = [
    GradedElement::theta(n, 0),
    GradedElement::theta(n, 1),
    GradedElement::theta(n, 2),
  ]
  let evens = [
    GradedElement::even(n, 0.5),
    GradedElement::even(n, -0.25),
    GradedElement::even(n, 1.5),
  ]
  for pa in 0..<2 {
    for pb in 0..<2 {
      for pc in 0..<2 {
        let x = if pa == 0 { evens[0] } else { odds[0] }
        let y = if pb == 0 { evens[1] } else { odds[1] }
        let z = if pc == 0 { evens[2] } else { odds[2] }
        assert_true(super_jacobiator(x, y, z).norm() < 1.0e-12)
        assert_true(derivation_defect(x, y, z).norm() < 1.0e-12)
      }
    }
  }
  // inhomogeneous (general) elements still satisfy linearity of the identity
  assert_true(super_jacobiator(a, b, c).norm() < 1.0e-12)
}

///|
test "T4: exp terminates on nilpotent soul" {
  let n = 3
  let x = SuperNumber::theta(n, 0).mul(SuperNumber::theta(n, 1))
  let e = x.exp()
  assert_true(approx(e.body_part(), 1.0))
  assert_true(approx(e.coeff[(1 << 0) | (1 << 1)], 1.0))
  // (θ1θ2)^2 = 0, so no higher terms appear
  assert_true(approx(e.soul_norm(), 1.0))
  let c = SuperNumber::body(n, 2.0)
  assert_true(approx(c.exp().body_part(), @math.exp(2.0), tol=1.0e-12))
}
