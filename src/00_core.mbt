// ============================================================================
// 00_core.mbt -- Shared algebraic / numeric primitives.
//
// Abstraction level 0: a mathematical structure is represented by
//   (carrier, operations, laws).  Here the carrier is fixed-length arrays
//   of Double; the operations are matrix/vector algebra; the laws are
//   checked in the inline tests of the later modules.
//
// MoonBit 0.10 syntax used here (verified against moonc v0.10.20, 2026-11-05):
//   * `#derive(...)` as an ATTRIBUTE above a type (new in v0.10.19)
//   * `pub extend T with Trait::{...}` (suppresses implicit-impl warnings)
//   * `Int32` / `Double` explicit numeric widths
//   * `pub(all) struct` for fully public constructors and fields
//   * `~` labelled arguments with defaults
// ============================================================================

///|
#derive(Debug, Eq)
pub(all) struct GradedDim {
  m : Int32 // even (bosonic) dimension
  n : Int32 // odd (fermionic) dimension
}

///|
pub extend GradedDim with Eq::{equal, not_equal}

///|
pub extend GradedDim with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// The superdimension (Euler characteristic of the structure sheaf).
pub fn GradedDim::super_dimension(self : GradedDim) -> Int32 {
  self.m - self.n
}

///|
/// dim of the Grassmann fibre Λ(n) = 2^n (Theorem T1).
pub fn GradedDim::grassmann_dim(self : GradedDim) -> Int32 {
  pow2(self.n)
}

///|
pub fn pow2(k : Int32) -> Int32 {
  let mut r = 1
  for _ in 0..<k {
    r *= 2
  }
  r
}

///|
pub fn popcount(mask : Int32) -> Int32 {
  mask.popcnt()
}

///|
/// Hardware-intrinsic 64-bit population count (`%u64.popcnt` in MoonBit 0.10.22).
pub fn popcount64(mask : UInt64) -> Int32 {
  mask.popcnt()
}

///|
/// Z_2 parity of a Grassmann monomial.
pub fn parity(mask : Int32) -> Int32 {
  popcount(mask) % 2
}

///|
/// Z_2 parity of a 64-bit Grassmann monomial (`n <= 64`).
pub fn parity64(mask : UInt64) -> Int32 {
  popcount64(mask) % 2
}

///|
/// IEEE-754 total-order sorting (`Double::total_compare`, new in MoonBit 0.10.19+),
/// deterministic even in the presence of `-0.0` or `NaN`.
pub fn sort_total(xs : Array[Double]) -> Array[Double] {
  let out = xs.copy()
  out.sort_by(Double::total_compare)
  out
}

///|
/// Relative numerical comparison (avoids IEEE NaN-as-key issues).
pub fn approx(a : Double, b : Double, tol~ : Double = 1.0e-9) -> Bool {
  (a - b).abs() <= tol * (1.0 + a.abs() + b.abs())
}

///|
pub fn clamp(x : Double, lo : Double, hi : Double) -> Double {
  if x < lo {
    lo
  } else if x > hi {
    hi
  } else {
    x
  }
}

///|
pub fn dot(x : Array[Double], y : Array[Double]) -> Double {
  let mut s = 0.0
  for i in 0..<x.length() {
    s += x[i] * y[i]
  }
  s
}

///|
pub fn vec_inf_norm(x : Array[Double]) -> Double {
  let mut s = 0.0
  for v in x {
    if v.abs() > s {
      s = v.abs()
    }
  }
  s
}

///|
pub fn vec_inf_diff(x : Array[Double], y : Array[Double]) -> Double {
  let mut s = 0.0
  for i in 0..<x.length() {
    let d = (x[i] - y[i]).abs()
    if d > s {
      s = d
    }
  }
  s
}

///|
pub fn vec_add(x : Array[Double], y : Array[Double]) -> Array[Double] {
  Array::makei(x.length(), fn(i) { x[i] + y[i] })
}

///|
pub fn vec_hadamard(x : Array[Double], y : Array[Double]) -> Array[Double] {
  Array::makei(x.length(), fn(i) { x[i] * y[i] })
}

///|
pub fn vec_pow_i(x : Array[Double], n : Int32) -> Array[Double] {
  Array::makei(x.length(), fn(i) {
    let mut r = 1.0
    for _ in 0..<n {
      r *= x[i]
    }
    r
  })
}

///|
pub fn vec_scale(x : Array[Double], s : Double) -> Array[Double] {
  Array::makei(x.length(), fn(i) { s * x[i] })
}

///|
pub fn normalize(x : Array[Double]) -> Array[Double] {
  let mut s = 0.0
  for v in x {
    s += v
  }
  if s == 0.0 {
    Array::make(x.length(), 1.0 / x.length().to_double())
  } else {
    Array::makei(x.length(), fn(i) { x[i] / s })
  }
}

///|
/// Deterministic LCG: reproducible agent initialisation without external deps.
pub(all) struct Rng {
  mut seed : Int32
}

///|
pub fn Rng::new(seed : Int32) -> Rng {
  Rng::{ seed, }
}

///|
pub fn Rng::unit(self : Rng) -> Double {
  self.seed = 1103515245 * self.seed + 12345
  let x = (self.seed >> 8) & 0x00FFFFFF
  x.to_double() / 16777216.0
}

///|
pub fn Rng::range(self : Rng, lo : Double, hi : Double) -> Double {
  lo + (hi - lo) * self.unit()
}

// ---------------------------------------------------------------------------
// Dense matrix algebra (row-major, length n*n for square matrices).
// ---------------------------------------------------------------------------

///|
pub fn mat_identity(n : Int32) -> Array[Double] {
  Array::makei(n * n, fn(k) { if k / n == k % n { 1.0 } else { 0.0 } })
}

///|
/// Matrix-vector product `a·v` for an `n × n` row-major `a`.
/// (Canonical location: every other module delegates here.)
pub fn matvec(a : Array[Double], v : Array[Double], n : Int32) -> Array[Double] {
  Array::makei(n, fn(i) {
    let mut s = 0.0
    for j in 0..<n {
      s += a[i * n + j] * v[j]
    }
    s
  })
}

///|
/// Gaussian-elimination rank of an `r × c` matrix, with a relative pivot
/// threshold.  Integer-valued matrices (the only ones used here) give exact
/// ranks for any `tol < 1`.
pub fn matrix_rank(
  a_in : Array[Double],
  r : Int32,
  c : Int32,
  tol~ : Double = 1.0e-9,
) -> Int32 {
  let a : Array[Double] = []
  for i in 0..<(r * c) {
    a.push(a_in[i])
  }
  let mut rank = 0
  let mut row = 0
  for col in 0..<c {
    if row >= r {
      break
    }
    // partial pivot
    let mut best = row
    let mut bestv = a[row * c + col].abs()
    for i in (row + 1)..<r {
      if a[i * c + col].abs() > bestv {
        best = i
        bestv = a[i * c + col].abs()
      }
    }
    if bestv <= tol {
      continue
    }
    if best != row {
      for j in 0..<c {
        let t = a[row * c + j]
        a[row * c + j] = a[best * c + j]
        a[best * c + j] = t
      }
    }
    let piv = a[row * c + col]
    for i in (row + 1)..<r {
      let f = a[i * c + col] / piv
      if f != 0.0 {
        for j in col..<c {
          a[i * c + j] -= f * a[row * c + j]
        }
      }
    }
    row += 1
    rank += 1
  }
  rank
}

///|
/// Nullity of an `r × c` matrix: `c - rank`.
pub fn matrix_nullity(a : Array[Double], r : Int32, c : Int32) -> Int32 {
  c - matrix_rank(a, r, c)
}

///|
/// Products of matrices in the same `n × n` shape, `k` times.
pub fn mat_power(a : Array[Double], n : Int32, k : Int32) -> Array[Double] {
  let mut out = mat_identity(n)
  for _ in 0..<k {
    out = matrix_product(out, n, n, a, n)
  }
  out
}

///|
/// Product of a (ra x ca) matrix and a (ca x cb) matrix.
pub fn matrix_product(
  a : Array[Double],
  ra : Int32,
  ca : Int32,
  b : Array[Double],
  cb : Int32,
) -> Array[Double] {
  Array::makei(ra * cb, fn(k) {
    let i = k / cb
    let j = k % cb
    let mut s = 0.0
    for t in 0..<ca {
      s += a[i * ca + t] * b[t * cb + j]
    }
    s
  })
}

///|
pub fn transpose(a : Array[Double], r : Int32, c : Int32) -> Array[Double] {
  Array::makei(c * r, fn(k) {
    let i = k / r
    let j = k % r
    a[j * c + i]
  })
}

///|
/// Determinant via Gaussian elimination with partial pivoting.
pub fn determinant(a : Array[Double], n : Int32) -> Double {
  let m = a.copy()
  let mut det = 1.0
  for k in 0..<n {
    let mut piv = k
    for i in (k + 1)..<n {
      if m[i * n + k].abs() > m[piv * n + k].abs() {
        piv = i
      }
    }
    if m[piv * n + k].abs() < 1.0e-14 {
      return 0.0
    }
    if piv != k {
      for j in 0..<n {
        let t = m[k * n + j]
        m[k * n + j] = m[piv * n + j]
        m[piv * n + j] = t
      }
      det = -det
    }
    let p = m[k * n + k]
    det *= p
    for i in (k + 1)..<n {
      let f = m[i * n + k] / p
      for j in k..<n {
        m[i * n + j] -= f * m[k * n + j]
      }
    }
  }
  det
}

///|
/// Solve A x = b for a nonsingular A (Gaussian elimination + back substitution).
pub fn gauss_solve(
  a : Array[Double],
  b : Array[Double],
  n : Int32,
) -> Array[Double] {
  let m = a.copy()
  let rhs = b.copy()
  for k in 0..<n {
    let mut piv = k
    for i in (k + 1)..<n {
      if m[i * n + k].abs() > m[piv * n + k].abs() {
        piv = i
      }
    }
    if piv != k {
      for j in 0..<n {
        let t = m[k * n + j]
        m[k * n + j] = m[piv * n + j]
        m[piv * n + j] = t
      }
      let t = rhs[k]
      rhs[k] = rhs[piv]
      rhs[piv] = t
    }
    for i in (k + 1)..<n {
      let f = m[i * n + k] / m[k * n + k]
      for j in k..<n {
        m[i * n + j] -= f * m[k * n + j]
      }
      rhs[i] -= f * rhs[k]
    }
  }
  let x = Array::make(n, 0.0)
  let mut i = n - 1
  while i >= 0 {
    let mut s = rhs[i]
    for j in (i + 1)..<n {
      s -= m[i * n + j] * x[j]
    }
    x[i] = s / m[i * n + i]
    i -= 1
  }
  x
}

///|
pub fn inverse(a : Array[Double], n : Int32) -> Array[Double] {
  let cols : Array[Array[Double]] = []
  for j in 0..<n {
    let e = Array::makei(n, fn(i) { if i == j { 1.0 } else { 0.0 } })
    cols.push(gauss_solve(a, e, n))
  }
  // columns of the inverse are the solutions above
  Array::makei(n * n, fn(k) {
    let i = k / n
    let j = k % n
    cols[j][i]
  })
}

///|
pub fn mat_inf_norm(a : Array[Double], n : Int32) -> Double {
  let mut s = 0.0
  for i in 0..<n {
    let mut r = 0.0
    for j in 0..<n {
      r += a[i * n + j].abs()
    }
    if r > s {
      s = r
    }
  }
  s
}

///|
test "core: graded dimension, popcount, determinant, inverse" {
  let d = GradedDim::{ m: 2, n: 3, }
  assert_eq(d.super_dimension(), -1)
  assert_eq(d.grassmann_dim(), 8)
  assert_eq(popcount(0b101101), 4)
  assert_eq(parity(0b101101), 0)
  assert_eq(popcount64(0x8000_0000_0000_000FUL), 5)
  assert_eq(parity64(0x8000_0000_0000_000FUL), 1)
  let ts = sort_total([3.0, -1.0, 0.0 / 0.0, -0.0, 2.0])
  assert_true(ts[0] != ts[0]) // NaN sorts first under Double::total_compare
  assert_true(approx(ts[1], -1.0))
  assert_true(approx(ts[4], 3.0))
  let a = [1.0, 2.0, 3.0, 4.0]
  assert_true(approx(determinant(a, 2), -2.0))
  let inv = inverse([4.0, 0.0, 0.0, 2.0], 2)
  assert_true(approx(inv[0], 0.25))
  assert_true(approx(inv[3], 0.5))
  assert_true(approx(gauss_solve([2.0, 0.0, 0.0, 4.0], [2.0, 8.0], 2)[1], 2.0))
}

///|
test "core: Rng is deterministic and in [0,1)" {
  let r1 = Rng::new(42)
  let r2 = Rng::new(42)
  let mut ok = true
  for _ in 0..<100 {
    let x = r1.unit()
    ok = ok && x >= 0.0 && x < 1.0 && x == r2.unit()
  }
  assert_true(ok)
}
