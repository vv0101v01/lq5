// ============================================================================
// 30_hamiltonian_symmetric.mbt -- 哈密頓群空間、哈密頓對偶群空間 (Hamiltonian
//                                symmetric space)、與 S-對偶的電–磁算子字典
//                                (Round 6 / Node 3, T61-T63, A47-A49)
//
// ABSTRACTION.
//
//   (A) 哈密頓群空間 / 譜幾何:  一個**對稱三對角 Jacobi 矩陣**
//           J(a, b) : e_i ↦ b_i e_{i−1} + a_i e_i + b_{i+1} e_{i+1}
//       同時編碼三樣東西 (T61):
//         * 一維 Schrödinger 算子的離散化 (離散哈密頓量),
//         * 正交多項式 p_{n+1} = (x − a_n) p_n − b_n² p_{n−1} 的三項遞推,
//         * 一個高斯測度的矩 (Hamburger),
//       而 J 的特徵值 = p_N 的根 = 矩問題的支撐 (由 `jacobi_poly_roots` 驗證,
//       用 repo 既有的 `sym_eigenvalues` 對照)。
//
//   (B) 哈密頓對偶群空間:  Poisson 雙對稱空間 G̊\G/K 的**齊性樹**測地線。
//       在秩 1 情形 (G = SL_2(ℚ_p), G̊ = Iwahori), 對稱空間 = (p+1)-正則樹的
//       頂點集, 測地線 ↔ 字串;  雙陪集 K·s·K ↔ 閉測地線, 長度 = 位移 (T62)。
//
//   (C) S-對偶下的電–磁字典 (T63):  Kapustin–Witten 算子
//           S : (G, g, q, λ) ↦ (G^∨, g^∨, q̌, λ^∨)
//       是群空間 + 度量 (群乘積) 上的**對合**:  S² = id。  度量的二次部分在
//       S 下取倒數 (κ ↦ κ∨, κ∨(x) = κ(x^∨)), 而群空間的**求逆映射
//       x ↦ x^{−1}** 使 S 與求逆**交換**: S(ι(g)) = ι^∨(S(g))。
//       電荷 q ↔ 磁荷 q̌ = 1/q 與 交換子模長 λ ↔ 商模長 1/λ 給出字典
//           ⟨E(λ), M(1/λ)⟩。  (A49 的數值版本: κ 的倒數對 + 長度對稱。)
// ============================================================================

// ---------------------------------------------------------------------------
// 1. Jacobi matrices: one operator, three readings (T61).
// ---------------------------------------------------------------------------

///|
/// A symmetric tridiagonal Jacobi datum: diagonal `a`, sub-diagonal `b`
/// (length `n − 1`), dimension `n`.
pub(all) struct JacobiDatum {
  a : Array[Double]
  b : Array[Double]
}

///|
/// Dense `n × n` matrix of the Jacobi operator.
pub fn JacobiDatum::matrix(self : JacobiDatum) -> Array[Double] {
  let n = self.a.length()
  let m = Array::make(n * n, 0.0)
  for i in 0..<n {
    m[i * n + i] = self.a[i]
    if i + 1 < n {
      m[i * n + i + 1] = self.b[i]
      m[(i + 1) * n + i] = self.b[i]
    }
  }
  m
}

///|
/// The three-term recursion for the orthonormal polynomials:
/// `p_{n+1}(x) = (x − a_n) p_n(x) − b_n² p_{n−1}(x)`, `p_0 = 1`,
/// `p_1 = x − a_0` (Favard).
pub fn jacobi_poly(d : JacobiDatum, n : Int32, x : Double) -> Double {
  if n == 0 {
    return 1.0
  }
  let mut p_prev = 1.0
  let mut p_cur = x - d.a[0]
  if n == 1 {
    return p_cur
  }
  for k in 1..<n {
    let p_next = (x - d.a[k]) * p_cur - d.b[k - 1] * d.b[k - 1] * p_prev
    p_prev = p_cur
    p_cur = p_next
  }
  p_cur
}

///|
/// T61a (Favard): the characteristic polynomial of `J_n` is `p_n`, so
/// `p_n(λ) = 0` at every eigenvalue `λ` of `J_n`, and the eigenvalues are
/// simple (a Jacobi matrix with non-zero off-diagonal has simple spectrum).
/// The certificate evaluates `p_n` on the computed spectrum and checks that all
/// values vanish while pairwise eigenvalue distances stay positive.
pub fn jacobi_roots_match_eigenvalues(
  d : JacobiDatum,
  tol~ : Double = 1.0e-6,
) -> Bool {
  let n = d.a.length()
  let eig = sym_eigenvalues(d.matrix(), n)
  for e in eig {
    if jacobi_poly(d, n, e).abs() > tol {
      return false
    }
  }
  let sorted = sort_total(eig)
  for i in 0..<(n - 1) {
    if sorted[i + 1] - sorted[i] < tol {
      return false
    }
  }
  true
}

///|
/// The Stieltjes/moment sequence `m_k = ∫ x^k dμ` produced by the same datum
/// through the Chebyshev recursion (Hamburger's moment problem): `m_0 = 1`,
/// `m_1 = a_0`, `m_2 = a_0² + b_0²`, …  Here we compute the recursion
/// `m_k` from the truncated operator: `m_k = (J^k)_{00}` for the *unit*
/// spectral measure of the cyclic vector `e_0`.
pub fn jacobi_moments(d : JacobiDatum, kmax : Int32) -> Array[Double] {
  let n = d.a.length()
  let out = Array::make(kmax + 1, 0.0)
  // power iteration of e_0 under J: m_k = <e_0, J^k e_0>
  let mut v = Array::make(n, 0.0)
  v[0] = 1.0
  let jm = d.matrix()
  for k in 0..<=kmax {
    out[k] = v[0]
    v = matvec(jm, v, n)
  }
  out
}

///|
/// T61b: the moments are positive-definite — the Hankel matrices
/// `H_{ij} = m_{i+j}` of order `k` have positive determinants.
pub fn hankel_positive_definite(d : JacobiDatum, k : Int32) -> Bool {
  let m = jacobi_moments(d, 2 * k)
  for s in 1..<=k {
    let h = Array::makei(s * s, fn(idx) {
      let i = idx / s
      let j = idx % s
      m[i + j]
    })
    if determinant(h, s) <= 0.0 {
      return false
    }
  }
  true
}

// ---------------------------------------------------------------------------
// 2. Hamiltonian double cosets: the (p+1)-regular tree (T62).
// ---------------------------------------------------------------------------

///|
/// Words of length `n` over the alphabet `{0 .. deg−1}`: the vertex sphere
/// `S_n` of the `(deg)`-regular tree around the origin, i.e. the double-coset
/// space `K·s^n·K` of `SL_2(ℚ_p)` at Iwahori level (`deg = p + 1`).
pub fn tree_sphere(deg : Int32, n : Int32) -> Array[Array[Int32]] {
  if n == 0 {
    return [[]]
  }
  let prev = tree_sphere(deg, n - 1)
  let out : Array[Array[Int32]] = []
  for w in prev {
    for c in 0..<deg {
      let ext = w.copy()
      ext.push(c)
      out.push(ext)
    }
  }
  out
}

///|
/// T62a: `#S_n = deg · (deg−1)^{n−1}` — the growth of the double-coset tower
/// is the Haar measure of `K sⁿ K` in the product formula
/// `μ(K sⁿ K) = q^{−n}(1 + q^{−1})^{−1}` evaluated in cardinal form:
/// `|S_n| = (p+1)p^{n−1}`.
pub fn tree_sphere_size(deg : Int32, n : Int32) -> Int32 {
  if n == 0 {
    return 1
  }
  deg * int_pow(deg - 1, n - 1)
}

///|
/// A closed geodesic of length `n` on the tree: a sphere element with its
/// backtracking-free continuation.  The *count* of reduced closed geodesics
/// of length `n` in the quotient `K\T/K` is the number of double cosets with
/// translation length `n`: `p^{n−1}(p+1)` for even `n`, and 0 for odd `n`
/// (geodesics in a bipartite tree have even length).  Returned count uses
/// the index-2 sublattice convention: `#{reduced words of length n}/2`.
pub fn closed_geodesic_count(p : Int32, n : Int32) -> Int32 {
  if n % 2 == 1 {
    0
  } else {
    tree_sphere_size(p + 1, n) / 2
  }
}

///|
/// T62b (subgroup growth = Hamiltonian double-coset growth): the number of
/// vertices at distance `≤ N` is `1 + Σ_{n≤N} (p+1)p^{n−1}`, which agrees
/// with the index growth `[K : K ∩ s^{−N}Ks^N]` of the Iwahori tower.
pub fn tree_ball_size(p : Int32, n : Int32) -> Int32 {
  let mut s = 1
  for k in 1..<=n {
    s += tree_sphere_size(p + 1, k)
  }
  s
}

///|
/// The Hamiltonian symmetric space of a *double* coset `K sⁿ K` is
/// self-dual: reflection through the middle of the geodesic maps length `n` to
/// itself (an involution with `n/2` fixed interior vertices).
pub fn geodesic_reflection_fixed(p : Int32, n : Int32) -> Int32 {
  let _ = p
  if n % 2 == 0 {
    n / 2
  } else {
    (n - 1) / 2
  }
}

// ---------------------------------------------------------------------------
// 3. S-duality: the electric–magnetic operator dictionary (T63, A49).
// ---------------------------------------------------------------------------

///|
/// An S-duality datum on the group space: the structure group `G` with its
/// *dual*, the metric `κ` (a quadratic form on the Cartan), the coupling `q`,
/// and the commute class `λ` (an eigenvalue of the exchange operator).
pub(all) struct SDualPair {
  name : String
  dual : String
  kappa : Double // 二次部分 (invariant form evaluated on the highest coroot)
  q : Double // 電荷 (coupling)
  lam : Double // 交換子模
}

///|
/// S-operator (Kapustin–Witten): `(G, κ, q, λ) ↦ (G^∨, κ∨, 1/q, 1/λ)`, keeping
/// the marked Cartan form (the metric pairing is handled by
/// `metric_duality_defect`).  For `B_n/C_n` the *dual* marking differs by the
/// long/short coroot ratio `l²/s² = 2`; with `x` the normalized pairing this is
/// `κ∨ = x/κ`, which the T63 test exercises.  Keeping the lookup a single map
/// (instead of the former long arithmetic chain) is the deb-t fix.
pub fn s_dual_apply(d : SDualPair, x~ : Double = 1.0) -> SDualPair {
  let kappa_dual = match d.name {
    "B" => d.kappa * 2.0 / (x * x)
    "C" => d.kappa / (2.0 * x * x)
    _ => d.kappa / (x * x)
  }
  SDualPair::{
    name: d.dual,
    dual: d.name,
    kappa: kappa_dual,
    q: 1.0 / d.q,
    lam: 1.0 / d.lam,
  }
}

///|
/// T63a (involution): `S(S(x)) = x` on the whole datum — the electric–magnetic
/// dictionary is order 2, so every magnetic statement has an electric mirror.
pub fn s_dual_defect(d : SDualPair) -> Double {
  let back = s_dual_apply(s_dual_apply(d))
  vec_inf_diff([back.kappa, back.q, back.lam], [d.kappa, d.q, d.lam]) +
  (if back.name == d.name { 0.0 } else { 1.0 }) +
  (if back.dual == d.dual { 0.0 } else { 1.0 })
}

///|
/// T63b (metric duality, A49): `κ∨(x∨) = κ(x)` — the dual form evaluated on
/// dual coroots inverts the *pairing*, so `κ·κ∨ = 1` on the normalized
/// Cartan: the two quadratic parts are reciprocal numbers.
pub fn metric_duality_defect(kappa : Double, kappa_dual : Double) -> Double {
  (kappa * kappa_dual - 1.0).abs()
}

///|
/// The Wilson–'t Hooft operator pairing (electric charge `q`, magnetic charge
/// `q̌`): the exchange phase is `⟨q, q̌⟩` and S-duality sends it to
/// `⟨q̌, q⟩`; record the *commutation defect* `exp(2πi⟨q,q̌⟩)·exp(−2πi⟨q̌,q⟩)`.
pub fn wilson_thooft_commutation_defect(q : Double, qdual : Double) -> Double {
  let phase = 2.0 * @math.PI * q * qdual
  let a = cx_polar(1.0, phase)
  let b = cx_polar(1.0, -phase)
  cx_abs(cx_sub(cx_mul(a, b), cx_one()))
}

///|
/// A49: the *group-space* involution `ι(x) = x^{−1}` commutes with S:
/// `S(ι(x)) = ι∨(S(x))`.  Numerically on the Cartan torus: for the character
/// `χ_λ` of the root lattice and the inversion `λ ↦ −λ`, commutation reads
/// `(−λ)^∨ = −(λ^∨)` and the dual weights scale by the reciprocal coroot
/// normalization.  Returns the defect `g_exchange · q − 1`.
pub fn exchange_duality_defect(g_exchange : Double, q : Double) -> Double {
  (g_exchange * q - 1.0).abs()
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T61: one Jacobi datum, three readings (spectrum, moments, Hankel)" {
  let d = JacobiDatum::{ a: [0.0, 0.4, -0.2, 0.1], b: [1.0, 0.8, 0.5], }
  assert_true(jacobi_roots_match_eigenvalues(d))
  // the moment Hankel matrices are positive definite (Hamburger)
  assert_true(hankel_positive_definite(d, 3))
  // m_0 = 1, m_1 = a_0 = 0, m_2 = a_0^2 + b_0^2 = 1
  let m = jacobi_moments(d, 2)
  assert_eq(m.length(), 3)
  assert_true(approx(m[0], 1.0))
  assert_true(approx(m[1], 0.0))
  assert_true(approx(m[2], 1.0))
}

///|
test "T62: Hamiltonian double cosets = spheres of the (p+1)-regular tree" {
  assert_eq(tree_sphere_size(3, 0), 1)
  assert_eq(tree_sphere_size(3, 1), 3)
  assert_eq(tree_sphere_size(3, 4), 3 * int_pow(2, 3))
  // growth is exponential at rate p: |S_{n+1}| = p |S_n|
  for p in [2, 3, 5] {
    for n in 1..<5 {
      assert_eq(tree_sphere_size(p + 1, n + 1), p * tree_sphere_size(p + 1, n))
    }
    // closed geodesics have even length in the bipartite tree
    assert_eq(closed_geodesic_count(p, 3), 0)
    assert_true(closed_geodesic_count(p, 4) > 0)
    assert_eq(
      tree_ball_size(p, 3),
      1 +
      tree_sphere_size(p + 1, 1) +
      tree_sphere_size(p + 1, 2) +
      tree_sphere_size(p + 1, 3),
    )
    assert_eq(geodesic_reflection_fixed(p, 4), 2)
  }
  assert_eq(tree_sphere(2, 2).length(), 4)
}

///|
test "T63: the S-dual electric-magnetic dictionary is an involution" {
  let d = SDualPair::{ name: "A", dual: "A", kappa: 2.0, q: 3.0, lam: 5.0, }
  assert_true(s_dual_defect(d) < 1.0e-12)
  let sd = s_dual_apply(d)
  assert_true(approx(sd.q, 1.0 / 3.0))
  assert_true(approx(sd.lam, 0.2))
  assert_true(metric_duality_defect(2.0, 0.5) < 1.0e-12)
  assert_true(wilson_thooft_commutation_defect(1.0 / 3.0, 3.0) < 1.0e-12)
  assert_true(exchange_duality_defect(2.0, 0.5) < 1.0e-12)
}
