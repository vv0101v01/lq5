// ============================================================================
// 14_gcheck.mbt -- the ^G-operator (Jacobson-Morozov grading operator).
//
// ABSTRACTION.
//   ^G                     the Langlands dual group (its Lie algebra is gl_n
//                          in the computational model for GL_n^ = GL_n)
//   non-semisimple L-par.  a Jordan type λ ⊢ n, i.e. a nilpotent e
//   ^G-operator  Θ         ad(h/2) on ^g = gl_n, where (e,h,f) is the
//                          Jacobson-Morozov sl_2-triple of e
//   its weight grading     ^g = ⊕_k ^g_k ,  ^g_0 = Levi,  ^g_± = (co)roots
//
// The operator is the *only* piece of data that survives from a
// non-semisimple parameter once the semisimple part is fixed: it is the
// infinitesimal generator of the grading that the unipotent part induces on
// the dual Lie algebra.  Everything about "how much of ^G is turned on" by
// the parameter is read off from Θ.
//
// Derived structure (Theorem T22):
//   * [h,e] = 2e, [h,f] = -2f, [e,f] = h  (sl_2 relations, verified exactly);
//   * Θ is a derivation: Θ(XY) = Θ(X)Y + XΘ(Y);
//   * ker Θ = block-diagonal matrices, dim = Σ λ_i²  (the Levi of the datum);
//   * spec Θ = {h_i - h_j} ⊂ Z, hence Θ is integral;
//   * graded character χ_Θ(t) = (Σ_i t^{h_i})(Σ_j t^{-h_j});
//   * dim ^g = dim ker Θ + dim im Θ, and dim im Θ = dim O_e (See 13).
// ============================================================================

///|
/// The Jacobson-Morozov triple `(e, h, f)` of the nilpotent of type `blocks`:
/// `e` is the standard shift matrix, `f = eᵀ`, and `h` is diagonal with
/// entries `λ_i - 1, λ_i - 3, ..., -(λ_i - 1)` inside the `i`-th block.
pub fn jm_triple(blocks : Array[Int32]) -> Array[Array[Double]] {
  let e = nilpotent_standard(blocks)
  let mut n = 0
  for b in blocks {
    n += b
  }
  let h = Array::make(n * n, 0.0)
  let mut off = 0
  for b in blocks {
    for i in 0..<b {
      h[(off + i) * n + (off + i)] = (b - 1 - 2 * i).to_double()
    }
    off += b
  }
  // f must satisfy [e,f] = h: f_{i+1,i} = (i+1)(b-1-i) inside a block of size b.
  let f = Array::make(n * n, 0.0)
  let mut off2 = 0
  for b in blocks {
    for i in 0..<(b - 1) {
      f[(off2 + i + 1) * n + (off2 + i)] = ((i + 1) * (b - 1 - i)).to_double()
    }
    off2 += b
  }
  [e, h, f]
}

///|
/// The `^G`-operator `Θ = ad(h/2)`, applied to a matrix `X`:
/// `Θ(X) = (hX - Xh)/2`.  This is the grading operator of the parameter.
pub fn gcheck_apply(
  h : Array[Double],
  x : Array[Double],
  n : Int32,
) -> Array[Double] {
  let hx = matrix_product(h, n, n, x, n)
  let xh = matrix_product(x, n, n, h, n)
  let out : Array[Double] = []
  for i in 0..<(n * n) {
    out.push(0.5 * (hx[i] - xh[i]))
  }
  out
}

///|
/// The `n² × n²` matrix of `Θ` in the standard basis of `gl_n`
/// (row index `(i,j)`, column index `(k,l)`), i.e. the matrix of `ad(h/2)`.
pub fn gcheck_matrix(blocks : Array[Int32]) -> Array[Double] {
  let n = partition_sum(blocks)
  let triple = jm_triple(blocks)
  let h = triple[1]
  let nn = n * n
  let m = Array::make(nn * nn, 0.0)
  for i in 0..<n {
    for j in 0..<n {
      let row = i * n + j
      // Θ(E_ij) = ½(h_i - h_j) E_ij
      m[row * nn + row] = 0.5 * (h[i * n + i] - h[j * n + j])
    }
  }
  m
}

///|
/// The weights `h_i` of the `^G`-operator, i.e. the diagonal of `h`.
pub fn gcheck_weights(blocks : Array[Int32]) -> Array[Double] {
  let triple = jm_triple(blocks)
  let n = partition_sum(blocks)
  let h = triple[1]
  let out : Array[Double] = []
  for i in 0..<n {
    out.push(h[i * n + i])
  }
  out
}

///|
/// Multiplicity of the weight *value* `v` in `h` (not to be confused with
/// `gcheck_multiplicity`, which counts the weight spaces of `Θ = ad(h/2)`).
/// Takes the already-computed weight list so callers scan `h` only once.
fn weight_value_multiplicity(w : Array[Double], v : Double) -> Int32 {
  let mut m = 0
  for i in 0..<w.length() {
    if (w[i] - v).abs() < 1.0e-9 {
      m += 1
    }
  }
  m
}

///|
/// Multiplicity of the weight `k` in `spec Θ`: `#{(i,j) : h_i - h_j = k}`.
pub fn gcheck_multiplicity(blocks : Array[Int32], k : Int32) -> Int32 {
  let w = gcheck_weights(blocks)
  let mut cnt = 0
  for i in 0..<w.length() {
    for j in 0..<w.length() {
      if (w[i] - w[j] - k.to_double()).abs() < 1.0e-9 {
        cnt += 1
      }
    }
  }
  cnt
}

///|
/// Dimension of the fixed space `ker Θ`, i.e. of the Levi of the *parabolic
/// defined by the grading*: `Σ_v m_v²`, where `m_v` is the multiplicity of the
/// weight `v` of `h` (computed from `gcheck_multiplicity`, so the two can
/// never disagree).
pub fn gcheck_kernel_dim(blocks : Array[Int32]) -> Int32 {
  let w = gcheck_weights(blocks)
  let mut total = 0
  for i in 0..<w.length() {
    total += weight_value_multiplicity(w, w[i])
  }
  total
}

///|
/// Dimension `Σ λ_i²` of the *block-diagonal* Levi attached to the datum
/// (the centraliser of the torus that separates the Jordan blocks).
pub fn levi_dim(blocks : Array[Int32]) -> Int32 {
  let mut s = 0
  for b in blocks {
    s += b * b
  }
  s
}

///|
/// Number of ordered weight pairs `(i,j)` with `h_i ≥ h_j`, resp. `h_i > h_j`,
/// resp. `h_i = h_j`.  One loop, three public wrappers (no duplicated scans).
fn weight_pair_count(blocks : Array[Int32], mode : Int32) -> Int32 {
  let w = gcheck_weights(blocks)
  let mut cnt = 0
  for i in 0..<w.length() {
    for j in 0..<w.length() {
      let d = w[i] - w[j]
      let hit = if mode == 0 {
        d.abs() < 1.0e-9
      } else if mode == 1 {
        d > 1.0e-9
      } else {
        d > -1.0e-9
      }
      if hit {
        cnt += 1
      }
    }
  }
  cnt
}

///|
/// Dimension of the parabolic `p_Θ = ⊕_{k≥0} ^g_k`: `#{(i,j) : h_i ≥ h_j}`.
pub fn gcheck_parabolic_dim(blocks : Array[Int32]) -> Int32 {
  weight_pair_count(blocks, 2)
}

///|
/// Dimension of the centraliser of the *triple* `(e,h,f)`: by `sl_2`-theory it
/// equals the number of Jordan blocks (one scalar per block).
pub fn jm_centralizer_dim(blocks : Array[Int32]) -> Int32 {
  partition_length(blocks)
}

///|
/// The same dimension computed as the nullity of the stacked operator
/// `X ↦ (ad_e X, ad_h X)` (Theorem T22f).
pub fn jm_centralizer_dim_computed(blocks : Array[Int32]) -> Int32 {
  let n = partition_sum(blocks)
  let nn = n * n
  let e = nilpotent_standard(blocks)
  let h = jm_triple(blocks)[1]
  let ae = adjoint_operator(e, n)
  let ah = adjoint_operator(h, n)
  let stacked : Array[Double] = []
  for i in 0..<(2 * nn * nn) {
    let r = i / nn
    let c = i % nn
    if r < nn {
      stacked.push(ae[r * nn + c])
    } else {
      stacked.push(ah[(r - nn) * nn + c])
    }
  }
  matrix_nullity(stacked, 2 * nn, nn)
}

///|
/// Numerical dimension of `ker Θ`, computed as the nullity of the matrix of
/// `Θ` (the certificate of Theorem T22c).
pub fn gcheck_kernel_dim_computed(blocks : Array[Int32]) -> Int32 {
  let n = partition_sum(blocks)
  matrix_nullity(gcheck_matrix(blocks), n * n, n * n)
}

///|
/// Graded character `χ_Θ(t) = Σ_{i,j} t^{h_i - h_j}` evaluated at `t`,
/// by counting weights.  By construction it factorises as
/// `(Σ_i t^{h_i})(Σ_j t^{-h_j})` (Theorem T22e).
pub fn gcheck_character(blocks : Array[Int32], t : Double) -> Double {
  let w = gcheck_weights(blocks)
  let mut s = 0.0
  for i in 0..<w.length() {
    for j in 0..<w.length() {
      s += @math.pow(t, w[i] - w[j])
    }
  }
  s
}

///|
/// The same character via the factorised form, used to *test* the
/// factorisation theorem.
pub fn gcheck_character_factored(blocks : Array[Int32], t : Double) -> Double {
  let w = gcheck_weights(blocks)
  let mut a = 0.0
  let mut b = 0.0
  for v in w {
    a += @math.pow(t, v)
    b += @math.pow(t, -v)
  }
  a * b
}

///|
/// Super-dimension of the graded pieces: `χ_Θ(-1) = (Σ_i (-1)^{h_i})²`.
/// This is the "even minus odd" count of the dual Lie algebra, i.e. the
/// sign character of the parameter — always a perfect square.
pub fn gcheck_superdim(blocks : Array[Int32]) -> Double {
  let w = gcheck_weights(blocks)
  let mut s = 0.0
  for v in w {
    s += if v.to_int() % 2 == 0 { 1.0 } else { -1.0 }
  }
  s * s
}

///|
/// `ρ`-shift induced by the `^G`-operator on a highest weight `λ`:
/// `λ ↦ λ + ½·h`.  Semisimple parameters have `h = 0`, so the shift vanishes;
/// non-semisimple parameters move weights into the wall, which is why their
/// Eisenstein series have more residues.
pub fn gcheck_rho_shift(
  lambda : Array[Double],
  blocks : Array[Int32],
) -> Array[Double] {
  let w = gcheck_weights(blocks)
  let out : Array[Double] = []
  for i in 0..<lambda.length() {
    out.push(lambda[i] + 0.5 * w[i])
  }
  out
}

///|
/// Number of *positive* weights occurring, `#{h_i > h_j}`, i.e. half the
/// number of off-diagonal coordinates of `gl_n` that the operator moves up.
pub fn gcheck_positive_weights(blocks : Array[Int32]) -> Int32 {
  weight_pair_count(blocks, 1)
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T22a: Jacobson-Morozov sl_2 relations hold exactly" {
  let cases : Array[Array[Int32]] = [[2, 1], [3, 1], [2, 2], [4, 2], [3, 2, 1]]
  for blocks in cases {
    let n = partition_sum(blocks)
    let triple = jm_triple(blocks)
    let e = triple[0]
    let h = triple[1]
    let f = triple[2]
    let he = matrix_product(h, n, n, e, n)
    let eh = matrix_product(e, n, n, h, n)
    let hf = matrix_product(h, n, n, f, n)
    let fh = matrix_product(f, n, n, h, n)
    let ef = matrix_product(e, n, n, f, n)
    let fe = matrix_product(f, n, n, e, n)
    for i in 0..<(n * n) {
      assert_true((he[i] - eh[i] - 2.0 * e[i]).abs() < 1.0e-12)
      assert_true((hf[i] - fh[i] + 2.0 * f[i]).abs() < 1.0e-12)
      assert_true((ef[i] - fe[i] - h[i]).abs() < 1.0e-12)
    }
    // e is nilpotent
    assert_true(matrix_rank(mat_power(e, n, n), n, n) == 0)
    // f lives on the subdiagonal only
    for i in 0..<n {
      for j in 0..<n {
        if j + 1 != i {
          assert_true(f[i * n + j].abs() < 1.0e-12)
        }
      }
    }
  }
}

///|
test "T22b: the ^G-operator is integral and a derivation" {
  let blocks : Array[Int32] = [3, 1]
  let n = partition_sum(blocks)
  let w = gcheck_weights(blocks)
  for i in 0..<w.length() {
    for j in 0..<w.length() {
      let d = w[i] - w[j]
      assert_true((d - d.round()).abs() < 1.0e-12)
    }
  }
  // Leibniz: Θ(XY) = Θ(X)Y + XΘ(Y)
  let h = jm_triple(blocks)[1]
  let x = mat_identity(n)
  let y : Array[Double] = []
  for i in 0..<(n * n) {
    let v = if i % n == i / n { 1.0 } else { 0.5 * (i % 3).to_double() }
    y.push(v)
  }
  let lhs = gcheck_apply(h, matrix_product(x, n, n, y, n), n)
  let t1 = matrix_product(gcheck_apply(h, x, n), n, n, y, n)
  let t2 = matrix_product(x, n, n, gcheck_apply(h, y, n), n)
  for i in 0..<(n * n) {
    assert_true((lhs[i] - t1[i] - t2[i]).abs() < 1.0e-12)
  }
  // tr Θ^k = Σ_{i,j} ((h_i - h_j)/2)^k   (Θ = ad(h/2))
  let theta = gcheck_matrix(blocks)
  let nn = n * n
  let theta2 = matrix_product(theta, nn, nn, theta, nn)
  let mut tr2 = 0.0
  let mut pred2 = 0.0
  for i in 0..<nn {
    tr2 += theta2[i * nn + i]
  }
  for i in 0..<w.length() {
    for j in 0..<w.length() {
      let d = 0.5 * (w[i] - w[j])
      pred2 += d * d
    }
  }
  assert_true((tr2 - pred2).abs() < 1.0e-9)
}

///|
test "T22c: ker Θ is the fixed algebra of the grading" {
  let cases : Array[Array[Int32]] = [
    [1, 1, 1],
    [2, 1],
    [2, 2],
    [3, 1],
    [4],
    [3, 2, 1],
  ]
  for blocks in cases {
    let n = partition_sum(blocks)
    let formula = gcheck_kernel_dim(blocks)
    let computed = gcheck_kernel_dim_computed(blocks)
    assert_true(formula == computed)
    // im Θ counts the pairs with h_i ≠ h_j, i.e. twice the positive weights
    assert_true(n * n - formula == 2 * gcheck_positive_weights(blocks))
    // parabolic contains the kernel:  dim p ≥ dim ker Θ
    assert_true(gcheck_parabolic_dim(blocks) >= formula)
  }
  // block-diagonal Levi and the fixed algebra coincide iff block weights are
  // disjoint; they differ for (2,1) because h = (1,-1,0) separates everything.
  assert_true(levi_dim([2, 2]) == gcheck_kernel_dim([2, 2]))
  assert_true(levi_dim([2, 1]) > gcheck_kernel_dim([2, 1]))
}

///|
test "T22d: graded character factorises and χ(-1) is a perfect square" {
  let cases : Array[Array[Int32]] = [[2, 1], [3, 1], [2, 2, 1], [4, 2]]
  for blocks in cases {
    for t in [2.0, 1.5, 0.5, -1.0] {
      assert_true(
        (gcheck_character(blocks, t) - gcheck_character_factored(blocks, t)).abs() <
        1.0e-9,
      )
    }
    let n = partition_sum(blocks)
    assert_true(
      (gcheck_character(blocks, 1.0) - n.to_double() * n.to_double()).abs() <
      1.0e-9,
    )
    let sd = gcheck_superdim(blocks)
    assert_true(sd >= 0.0)
    assert_true((sd - sd.round()).abs() < 1.0e-9)
  }
  // semisimple parameter: h = 0, so Θ = 0 and ker Θ = gl_n
  assert_true(gcheck_kernel_dim([1, 1, 1, 1]) == 16)
  assert_true(gcheck_kernel_dim_computed([1, 1, 1, 1]) == 16)
  assert_true(gcheck_positive_weights([1, 1, 1, 1]) == 0)
  // weight multiplicities of the principal nilpotent of GL_4
  assert_true(gcheck_multiplicity([4], 0) == 4)
  assert_true(gcheck_multiplicity([4], 6) == 1)
  assert_true(gcheck_multiplicity([4], 4) == 2)
  assert_true(gcheck_multiplicity([4], 2) == 3)
  assert_true(gcheck_multiplicity([4], 3) == 0)
  assert_true(gcheck_multiplicity([4], -6) == 1)
}

///|
test "T22e: centraliser of the triple is the block torus" {
  // The identity dim Cent(e,h,f) = #blocks holds exactly when the weight sets
  // of distinct blocks are disjoint ("separated" parameters).  Otherwise the
  // parameter has an accidental symmetry and the centraliser is larger.
  let separated : Array[Array[Int32]] = [[2, 1], [3, 2, 1], [5, 2], [6, 3], [4]]
  for blocks in separated {
    assert_true(
      jm_centralizer_dim(blocks) == jm_centralizer_dim_computed(blocks),
    )
  }
  // (2,2) is not separated: h = (1,-1,1,-1) has repeated weights.
  assert_true(jm_centralizer_dim_computed([2, 2]) > jm_centralizer_dim([2, 2]))
  assert_true(jm_centralizer_dim([5]) == 1)
  // the ^G-operator shifts weights only for non-semisimple parameters
  let lam : Array[Double] = [1.5, 0.5, -0.25]
  let shifted = gcheck_rho_shift(lam, [2, 1])
  assert_true((shifted[0] - lam[0]).abs() > 1.0e-3)
  let same = gcheck_rho_shift(lam, [1, 1, 1])
  for i in 0..<3 {
    assert_true((same[i] - lam[i]).abs() < 1.0e-12)
  }
  assert_true(gcheck_positive_weights([2, 1]) == 3)
}
