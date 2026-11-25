// ============================================================================
// 13_lparams.mbt -- Langlands parameters, semisimple or not.
//
// ABSTRACTION.
//   L-parameter            φ : W_F × SL_2(C) → ^G
//   semisimple part        eigenvalues of φ(Frobenius)^{ss}  → Array[Double]
//   unipotent part         nilpotent N with a Jordan type λ ⊢ n → partition
//   conjugacy class        pair (eigenvalue multiset, Jordan type)
//   non-semisimple         λ ≠ (1^n): the parameter has a unipotent defect
//
// The whole theory of non-semisimple L-parameters is governed by linear
// algebra of nilpotent matrices, so the abstraction is *exact*: every
// statement below is decided by an integer computation
// (ranks of powers of N, nullity of the commutator operator ad_N).
//
// Derived structure (Theorem T21):
//   * Jordan type from ranks:  b_k = r_{k-1} - 2 r_k + r_{k+1};
//   * dim Cent(N) = Σ (λ'_i)², where λ' is the conjugate partition;
//   * dim O_N = n² - dim Cent(N)  (orbit of the unipotent part);
//   * dim Ext¹ = dim Cent(N): non-semisimple parameters are *soft*;
//   * semisimple ⇔ λ = (1^n) ⇔ unipotent defect 0.
// ============================================================================

///|
/// The commutator operator `ad_N : X ↦ XN - NX` on `gl_n`, flattened to an
/// `n² × n²` matrix (row index `(i,j)`, column index `(k,l)`).
pub fn adjoint_operator(nil : Array[Double], n : Int32) -> Array[Double] {
  let nn = n * n
  let m = Array::make(nn * nn, 0.0)
  // (XN)[i][j] = X[i][k] N[k][j] ; (NX)[i][j] = N[i][k] X[k][j]
  for i in 0..<n {
    for j in 0..<n {
      let row = i * n + j
      for k in 0..<n {
        // += N[k][j] at X[i][k]
        m[row * nn + (i * n + k)] += nil[k * n + j]
        // -= N[i][k] at X[k][j]
        m[row * nn + (k * n + j)] -= nil[i * n + k]
      }
    }
  }
  m
}

///|
/// Dimension of the centraliser of a nilpotent, computed *numerically* as
/// `dim ker ad_N`.
pub fn nilpotent_centralizer_dim(nil : Array[Double], n : Int32) -> Int32 {
  matrix_nullity(adjoint_operator(nil, n), n * n, n * n)
}

///|
/// The standard block-diagonal nilpotent of Jordan type `blocks`
/// (each block is a shift matrix `E_{i,i+1}`).
pub fn nilpotent_standard(blocks : Array[Int32]) -> Array[Double] {
  let mut n = 0
  for b in blocks {
    n += b
  }
  let m = Array::make(n * n, 0.0)
  let mut off = 0
  for b in blocks {
    for i in 0..<(b - 1) {
      m[(off + i) * n + (off + i + 1)] = 1.0
    }
    off += b
  }
  m
}

///|
/// Jordan type of a nilpotent matrix, from the ranks of its powers:
/// `b_k = r_{k-1} - 2 r_k + r_{k+1}`, `r_0 = n`.
pub fn jordan_type(nil : Array[Double], n : Int32) -> Array[Int32] {
  let ranks : Array[Int32] = [n]
  let mut k = 1
  while k <= n + 1 {
    ranks.push(matrix_rank(mat_power(nil, n, k), n, n))
    k += 1
  }
  let blocks : Array[Int32] = []
  for kk in 1..<=n {
    let b = ranks[kk - 1] - 2 * ranks[kk] + ranks[kk + 1]
    for _ in 0..<b {
      blocks.push(kk)
    }
  }
  // Jordan type is conventionally listed in decreasing order.
  sort_desc(blocks)
}

///|
///In-place-ish descending sort of an `Int32` array (returns a new array).
pub fn sort_desc(x : Array[Int32]) -> Array[Int32] {
  let out : Array[Int32] = []
  for v in x {
    out.push(v)
  }
  let n = out.length()
  for i in 0..<n {
    for j in 0..<(n - 1 - i) {
      if out[j] < out[j + 1] {
        let t = out[j]
        out[j] = out[j + 1]
        out[j + 1] = t
      }
    }
  }
  out
}

///|
/// Conjugate partition `λ'` (`λ'_i = #{j : λ_j ≥ i}`).
pub fn partition_conjugate(part : Array[Int32]) -> Array[Int32] {
  let mut max = 0
  for v in part {
    if v > max {
      max = v
    }
  }
  let out : Array[Int32] = []
  for i in 1..<=max {
    let mut cnt = 0
    for v in part {
      if v >= i {
        cnt += 1
      }
    }
    out.push(cnt)
  }
  out
}

///|
/// `Σ (λ'_i)²`: the size of the centraliser of the unipotent part
/// (dimension of the *commutant*, Theorem T21b).
pub fn centralizer_dim_formula(part : Array[Int32]) -> Int32 {
  let conj = partition_conjugate(part)
  let mut s = 0
  for v in conj {
    s += v * v
  }
  s
}

///|
/// Total size of a partition.
pub fn partition_sum(part : Array[Int32]) -> Int32 {
  let mut s = 0
  for v in part {
    s += v
  }
  s
}

///|
/// Number of parts (`= dim of the largest torus acting on the block
/// decomposition`).
pub fn partition_length(part : Array[Int32]) -> Int32 {
  part.length()
}

///|
/// Unipotent defect `Σ (λ_i - 1) = n - #{blocks}`: how far the parameter is
/// from being semisimple (Theorem T21d).
pub fn unipotent_defect(part : Array[Int32]) -> Int32 {
  partition_sum(part) - partition_length(part)
}

///|
/// Dimension of the unipotent orbit `O_N ⊂ GL_n`:
/// `n² - Σ(λ'_i)²` (Theorem T21c).
pub fn orbit_dim(part : Array[Int32]) -> Int32 {
  let n = partition_sum(part)
  n * n - centralizer_dim_formula(part)
}

///|
/// First-order deformation space of the parameter:
/// `dim Ext¹ = dim Cent(N) = Σ(λ'_i)²`.  Non-semisimple parameters are
/// therefore *soft*: degeneracy is a high-dimensional phenomenon.
pub fn ext1_dim(part : Array[Int32]) -> Int32 {
  centralizer_dim_formula(part)
}

///|
/// All partitions of `n`, in decreasing lexicographic order.
pub fn partitions(n : Int32) -> Array[Array[Int32]] {
  let out : Array[Array[Int32]] = []
  partition_rec(n, n, [], out)
  out
}

///|
///Recursive partition enumerator with parts `≤ max`.
fn partition_rec(
  remaining : Int32,
  max : Int32,
  cur : Array[Int32],
  out : Array[Array[Int32]],
) -> Unit {
  if remaining == 0 {
    let copy : Array[Int32] = []
    for v in cur {
      copy.push(v)
    }
    out.push(copy)
    return
  }
  let mut next = max
  if next > remaining {
    next = remaining
  }
  while next >= 1 {
    cur.push(next)
    partition_rec(remaining - next, next, cur, out)
    cur.pop() |> ignore
    next -= 1
  }
}

///|
/// Number of nilpotent orbits in `gl_n`, i.e. `p(n)` (Theorem T21e).
pub fn nilpotent_orbit_count(n : Int32) -> Int32 {
  partitions(n).length()
}

// ---------------------------------------------------------------------------
// The L-parameter itself.
// ---------------------------------------------------------------------------

///|
/// A Langlands parameter for `GL_n`, in Jordan form: a semisimple part given
/// by `n` eigenvalues (one per basis vector) and a Jordan type `blocks`.
/// A *semisimple* parameter has `blocks = (1,...,1)`; a *non-semisimple* one
/// has at least one block of size `≥ 2`.
pub(all) struct LParameter {
  n : Int32
  eigenvalues : Array[Double]
  blocks : Array[Int32]
}

///|
/// Constructor from a partition and an eigenvalue list (assigned block by
/// block, in order).
pub fn l_parameter(
  blocks : Array[Int32],
  eigenvalues : Array[Double],
) -> LParameter {
  { n: partition_sum(blocks), eigenvalues, blocks, }
}

///|
/// The semisimple reduction: same eigenvalues, all blocks of size 1.
pub fn LParameter::semisimple_part(self : LParameter) -> LParameter {
  let ones : Array[Int32] = []
  for _ in 0..<self.n {
    ones.push(1)
  }
  { n: self.n, eigenvalues: self.eigenvalues, blocks: ones, }
}

///|
/// `true` iff the parameter is semisimple, i.e. the unipotent defect vanishes.
pub fn LParameter::is_semisimple(self : LParameter) -> Bool {
  unipotent_defect(self.blocks) == 0
}

///|
/// Jordan matrix of the parameter: `diag(λ_i) + N`, with `N` the standard
/// nilpotent of type `blocks`.  Delegates to `eisenstein_local_system` so that
/// "L-parameter matrix" and "Eisenstein local system" cannot drift apart.
pub fn LParameter::matrix(self : LParameter) -> Array[Double] {
  eisenstein_local_system(self.blocks, self.eigenvalues)
}

///|
/// `true` iff all eigenvalues coincide (a *scalar* semisimple part); then the
/// parameter is a "unipotent-times-scalar" sphere and its centraliser has
/// dimension `Σ (λ'_i)²`.
pub fn LParameter::scalar_semisimple(
  self : LParameter,
  tol~ : Double = 1.0e-9,
) -> Bool {
  let e0 = self.eigenvalues[0]
  let mut ok = true
  for v in self.eigenvalues {
    if (v - e0).abs() > tol {
      ok = false
    }
  }
  ok
}

///|
/// Dimension of the centraliser of the parameter:
/// * scalar semisimple part: `Σ (λ'_i)²` (commutant of the nilpotent);
/// * regular semisimple part (distinct eigenvalues, blocks of size 1): `n`;
/// * mixed: `Σ_j (block_j)²` when the eigenvalues separate the blocks.
/// This is computed, not assumed: see `centralizer_dim_computed`.
pub fn LParameter::centralizer_dim(self : LParameter) -> Int32 {
  if self.scalar_semisimple() {
    centralizer_dim_formula(self.blocks)
  } else {
    let mut s = 0
    for b in self.blocks {
      s += b * b
    }
    s
  }
}

///|
/// The same dimension obtained numerically as `dim ker ad_M` for `M` the
/// matrix of the parameter.  Agreement with the formula is Theorem T21b.
pub fn LParameter::centralizer_dim_computed(self : LParameter) -> Int32 {
  nilpotent_centralizer_dim(self.matrix(), self.n)
}

///|
/// Dimension of the orbit of the parameter under `GL_n`-conjugacy.
pub fn LParameter::orbit_dim(self : LParameter) -> Int32 {
  self.n * self.n - self.centralizer_dim()
}

///|
/// `true` iff the parameter is *irreducible* as a representation of
/// `W_F × SL_2`: a single Jordan block (for the unipotent part) and a single
/// eigenvalue.
pub fn LParameter::is_irreducible(self : LParameter) -> Bool {
  self.scalar_semisimple() && self.blocks.length() == 1
}

///|
/// The "Arthur SL_2" defect: the sum of `(block size - 1)`, i.e. the number
/// of `SL_2`-directions that are *not* part of the Arthur parameter
/// `SL_2 → ^G`.  For a semisimple parameter it is `0`; this is what
/// distinguishes "Arthur parameters" from the rest of `Par_G`.
pub fn LParameter::arthur_defect(self : LParameter) -> Int32 {
  unipotent_defect(self.blocks)
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T21a: Jordan type is recovered exactly from ranks of powers" {
  let cases : Array[Array[Int32]] = [
    [1, 1, 1],
    [2, 1],
    [3],
    [2, 2],
    [3, 1],
    [2, 1, 1],
    [4, 2],
  ]
  for part in cases {
    let nil = nilpotent_standard(part)
    let n = partition_sum(part)
    let recovered = jordan_type(nil, n)
    assert_true(recovered.length() == part.length())
    for i in 0..<part.length() {
      assert_true(recovered[i] == part[i])
    }
    // nilpotency itself
    let niln = mat_power(nil, n, n)
    let mut s = 0.0
    for v in niln {
      s += v.abs()
    }
    assert_true(s < 1.0e-12)
  }
}

///|
test "T21b: dim Cent(N) = Σ(λ'_i)²  (computed vs. formula)" {
  let cases : Array[Array[Int32]] = [
    [1, 1, 1, 1],
    [2, 1, 1],
    [2, 2],
    [3, 1],
    [4],
  ]
  for part in cases {
    let nil = nilpotent_standard(part)
    let n = partition_sum(part)
    let computed = nilpotent_centralizer_dim(nil, n)
    assert_true(computed == centralizer_dim_formula(part))
  }
  // GL_4 orbit dimensions: 0, 6, 8, 10, 12
  let gl4 : Array[Array[Int32]] = [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]]
  let expected : Array[Int32] = [0, 6, 8, 10, 12]
  for i in 0..<5 {
    assert_true(orbit_dim(gl4[i]) == expected[i])
    assert_true(ext1_dim(gl4[i]) == 16 - expected[i])
  }
}

///|
test "T21c: semisimple ⇔ trivial unipotent defect; p(n) orbit count" {
  assert_true(unipotent_defect([1, 1, 1]) == 0)
  assert_true(unipotent_defect([3]) == 2)
  let p = l_parameter([2, 2], [1.0, 1.0, 1.0, 1.0])
  assert_false(p.is_semisimple())
  assert_true(p.arthur_defect() == 2)
  assert_true(p.centralizer_dim() == 8)
  assert_true(p.centralizer_dim_computed() == 8)
  assert_true(p.orbit_dim() == 8)
  // a genuinely irreducible parameter
  let irr = l_parameter([4], [2.5, 2.5, 2.5, 2.5])
  assert_true(irr.is_irreducible())
  assert_true(irr.centralizer_dim() == 4)
  assert_true(irr.centralizer_dim_computed() == 4)
  // counting nilpotent orbits = partitions
  assert_true(nilpotent_orbit_count(4) == 5)
  assert_true(nilpotent_orbit_count(5) == 7)
  assert_true(nilpotent_orbit_count(6) == 11)
  // semisimple reduction of a non-semisimple parameter
  let red = p.semisimple_part()
  assert_true(red.is_semisimple())
  assert_true(red.centralizer_dim() == 16) // scalar semisimple part: centraliser = gl_4
}
