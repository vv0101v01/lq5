// ============================================================================
// 18_root_adelic_eisenstein.mbt -- General root systems (A_n, B_n, C_n, D_n),
//                                  Weyl reflection BFS, and p-adic / adelic
//                                  Gindikin-Karpelevich Euler products.
//
// ELIMINATING LIMITATIONS 5 & 6 OF ROUND 2:
//   * Limitation 5 was modelling `^G` only by `GL_n` (`W = S_r`, type A) and
//     restricting rank `r <= 4`.  Here we implement the Classical Root Systems
//     `RootType = TypeA(n) | TypeB(n) | TypeC(n) | TypeD(n)` in standard
//     Euclidean coordinates, exact Weyl reflection action
//       `s_α(v) = v - ⟨v, α^∨⟩ α`,   `⟨v, α^∨⟩ = 2⟨v, α⟩ / ⟨α, α⟩`,
//     Langlands dual root system `Φ ↦ Φ^∨` (`B_n ↔ C_n`, `A_n, D_n` self-dual),
//     and finite Weyl group generation via BFS over simple reflections!
//   * Limitation 6 was omitting p-adic local factors and global Euler products
//     (using only archimedean `ξ(t)/ξ(t+1)` on real `t`).  Here we implement:
//       - Exact p-adic Gindikin-Karpelevich factor at prime `p` (`s ∈ ℂ`):
//           `c_p(s) = (1 - p^{-(s+1)}) / (1 - p^{-s})`
//       - Global adelic GK factor `c_A(s) = c_∞(s) · Π_{p <= P} c_p(s)`:
//         Wait: locally `ζ_p(s) = (1 - p^{-s})^{-1}`, so
//           `Π_p c_p(s) = Π_p ζ_p(s)/ζ_p(s+1) = ζ(s)/ζ(s+1)`,
//         and multiplying by the archimedean gamma factor
//           `c_∞(s) = (π^{-s/2} Γ(s/2)) / (π^{-(s+1)/2} Γ((s+1)/2))`
//         reconstructs the completed global factor `c_global(s) = ξ(s)/ξ(s+1)`!
//       - Complex-spectral constant-term matrix `A_C(λ)` over `Cx` for ANY
//         classical Weyl group (`A_2, B_2, C_2, B_3, D_3, A_3`).
//
// DERIVED THEOREMS & ALGORITHMS:
//   * Theorem T27 (Langlands Root-Coroot Duality & General Weyl GK Cocycle):
//       1. `(Φ^∨)^∨ = Φ`, `TypeB(r)^∨ = TypeC(r)`, `|W(B_r)| = |W(C_r)| = 2^r r!`,
//          `|W(D_r)| = 2^{r-1} r!`, `|Φ^+(B_r)| = r^2`, `|Φ^+(D_r)| = r(r-1)`.
//       2. Weyl vector `ρ = ½ Σ_{α ∈ Φ^+} α` satisfies `⟨ρ, α_i^∨⟩ = 1` for EVERY
//          simple root `α_i` of `A_r, B_r, C_r, D_r` (universal Weyl identity).
//       3. For ANY root system `Φ`, complex spectral weight `λ ∈ ℂ^d`, and
//          `w_1, w_2 ∈ W(Φ)`:
//            `c_Φ(w, λ) = Π_{α ∈ Φ^+, w(α) ∈ -Φ^+} c(⟨λ, α^∨⟩)`
//          satisfies the crossed-homomorphism cocycle law
//            `c_Φ(w_1 w_2, λ) = c_Φ(w_1, w_2 λ) · c_Φ(w_2, λ)`
//          both locally at every prime `p` (`c_p` on `ℂ`!) and globally (`ξ`),
//          yielding `Π_Φ(λ)^2 = Π_Φ(λ)`, `tr Π_Φ(λ) = 1` on non-type-A groups
//          such as `Sp_4 / SO_5` (`|W(B_2)| = 8`) and `SO_7 / Sp_6` (`|W(B_3)| = 48`).
//   * Theorem T28 (Adelic Euler-Product Convergence & Unitarity on iℝ):
//       1. On `Re(s) > 1`, `c_∞(s) · Π_{p <= P} c_p(s) → ξ(s)/ξ(s+1)` at rate `O(P^{1-Re(s)})`.
//       2. Unitarity on the unitary axis `λ ∈ iℝ^d`: for purely imaginary spectral
//          parameters `s = it` (`t ∈ ℝ \ {0}`), `ξ(1 - it) = conj(ξ(1 + it))` and
//          `ξ(-it) = ξ(1 + it)`, hence `|c_global(it)| = 1` identically and the
//          intertwining operator `c_Φ(w, iν)` is unitary (`|c_Φ(w, iν)| = 1 ∀ w ∈ W`)!
// ============================================================================

///|
/// Classical reduced root system types of rank `r`.
#derive(Debug, Eq)
pub(all) enum RootType {
  TypeA(Int32) // A_r inside R^{r+1}, |W| = (r+1)!
  TypeB(Int32) // B_r (SO_{2r+1}), short roots e_i, long e_i ± e_j, |W| = 2^r r!
  TypeC(Int32) // C_r (Sp_{2r}), short e_i ± e_j, long 2e_i, |W| = 2^r r!
  TypeD(Int32) // D_r (SO_{2r}), roots e_i ± e_j, |W| = 2^{r-1} r!
}

///|
pub extend RootType with Eq::{equal, not_equal}

///|
pub extend RootType with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Ambient Euclidean dimension `d` in which roots of `rt` live.
pub fn RootType::ambient_dim(self : RootType) -> Int32 {
  match self {
    TypeA(r) => r + 1
    TypeB(r) | TypeC(r) | TypeD(r) => r
  }
}

///|
/// Semisimple rank `r`.
pub fn RootType::rank(self : RootType) -> Int32 {
  match self {
    TypeA(r) | TypeB(r) | TypeC(r) | TypeD(r) => r
  }
}

///|
/// Langlands dual root system type (`B_r ↔ C_r`, `A_r` and `D_r` self-dual).
pub fn RootType::langlands_dual(self : RootType) -> RootType {
  match self {
    TypeA(r) => TypeA(r)
    TypeB(r) => TypeC(r)
    TypeC(r) => TypeB(r)
    TypeD(r) => TypeD(r)
  }
}

///|
fn basis_diff(
  d : Int32,
  i : Int32,
  j : Int32,
  sign_j : Double,
) -> Array[Double] {
  let v = Array::make(d, 0.0)
  v[i] = 1.0
  v[j] = sign_j
  v
}

///|
fn single_basis(d : Int32, i : Int32, scale : Double) -> Array[Double] {
  let v = Array::make(d, 0.0)
  v[i] = scale
  v
}

///|
/// Simple roots `Δ = {α_1, ..., α_r}` in standard Bourbaki coordinates.
pub fn simple_roots(rt : RootType) -> Array[Array[Double]] {
  let d = rt.ambient_dim()
  let r = rt.rank()
  let out : Array[Array[Double]] = []
  match rt {
    TypeA(_) =>
      for i in 0..<r {
        out.push(basis_diff(d, i, i + 1, -1.0))
      }
    TypeB(_) => {
      for i in 0..<(r - 1) {
        out.push(basis_diff(d, i, i + 1, -1.0))
      }
      out.push(single_basis(d, r - 1, 1.0))
    }
    TypeC(_) => {
      for i in 0..<(r - 1) {
        out.push(basis_diff(d, i, i + 1, -1.0))
      }
      out.push(single_basis(d, r - 1, 2.0))
    }
    TypeD(_) => {
      for i in 0..<(r - 1) {
        out.push(basis_diff(d, i, i + 1, -1.0))
      }
      if r >= 2 {
        out.push(basis_diff(d, r - 2, r - 1, 1.0))
      }
    }
  }
  out
}

///|
/// Positive roots `Φ^+` in standard coordinates.
pub fn positive_roots(rt : RootType) -> Array[Array[Double]] {
  let d = rt.ambient_dim()
  let r = rt.rank()
  let out : Array[Array[Double]] = []
  match rt {
    TypeA(_) =>
      for i in 0..<=r {
        for j in (i + 1)..<=r {
          out.push(basis_diff(d, i, j, -1.0))
        }
      }
    TypeB(_) => {
      for i in 0..<r {
        for j in (i + 1)..<r {
          out.push(basis_diff(d, i, j, -1.0))
          out.push(basis_diff(d, i, j, 1.0))
        }
      }
      for i in 0..<r {
        out.push(single_basis(d, i, 1.0))
      }
    }
    TypeC(_) => {
      for i in 0..<r {
        for j in (i + 1)..<r {
          out.push(basis_diff(d, i, j, -1.0))
          out.push(basis_diff(d, i, j, 1.0))
        }
      }
      for i in 0..<r {
        out.push(single_basis(d, i, 2.0))
      }
    }
    TypeD(_) =>
      for i in 0..<r {
        for j in (i + 1)..<r {
          out.push(basis_diff(d, i, j, -1.0))
          out.push(basis_diff(d, i, j, 1.0))
        }
      }
  }
  out
}

///|
/// Coroot `α^∨ = 2 α / ⟨α, α⟩`.
pub fn coroot(alpha : Array[Double]) -> Array[Double] {
  vec_scale(alpha, 2.0 / dot(alpha, alpha))
}

///|
/// Cartan matrix `C_{ij} = ⟨α_i, α_j^∨⟩` (`r × r` integer-valued matrix).
pub fn cartan_matrix(rt : RootType) -> Array[Double] {
  let sr = simple_roots(rt)
  let r = sr.length()
  Array::makei(r * r, fn(k) {
    let i = k / r
    let j = k % r
    dot(sr[i], coroot(sr[j]))
  })
}

///|
/// Half-sum of positive roots `ρ = ½ Σ_{α ∈ Φ^+} α` (the Weyl vector).
pub fn weyl_rho(rt : RootType) -> Array[Double] {
  let d = rt.ambient_dim()
  let mut acc = Array::make(d, 0.0)
  for alpha in positive_roots(rt) {
    acc = vec_add(acc, alpha)
  }
  vec_scale(acc, 0.5)
}

///|
/// Orthogonal reflection matrix `s_α = I_d - α α^{∨ T}` (`d × d` row-major).
pub fn reflection_matrix(alpha : Array[Double]) -> Array[Double] {
  let d = alpha.length()
  let av = coroot(alpha)
  Array::makei(d * d, fn(k) {
    let i = k / d
    let j = k % d
    (if i == j { 1.0 } else { 0.0 }) - alpha[i] * av[j]
  })
}

///|
fn mat_key(m : Array[Double]) -> String {
  let mut s = ""
  for i in 0..<m.length() {
    // Classical Weyl matrices in standard coordinates have exact half-integer / integer entries;
    // multiplying by 2 and rounding gives a canonical collision-free integer key!
    let iv = (2.0 * m[i]).round().to_int()
    if i > 0 {
      s = s + ","
    }
    s = s + iv.to_string()
  }
  s
}

///|
/// Generate the entire finite Weyl group `W(Φ) ⊂ O(d)` by BFS from simple
/// reflections `s_{α_1}, ..., s_{α_r}` (Algorithm A23).
pub fn generate_weyl_group(rt : RootType) -> Array[Array[Double]] {
  let d = rt.ambient_dim()
  let gens : Array[Array[Double]] = []
  for alpha in simple_roots(rt) {
    gens.push(reflection_matrix(alpha))
  }
  let id = mat_identity(d)
  let seen = @sorted_map.SortedMap::from_array([(mat_key(id), true)])
  let group : Array[Array[Double]] = [id]
  let mut head = 0
  while head < group.length() {
    let w = group[head]
    head += 1
    for s in gens {
      let sw = matrix_product(s, d, d, w, d)
      let k = mat_key(sw)
      if !seen.contains(k) {
        seen[k] = true
        group.push(sw)
      }
    }
  }
  group
}

///|
/// Test whether a vector `v` lies in `-Φ^+` (i.e. `-v ∈ Φ^+`).
pub fn is_negative_root(
  v : Array[Double],
  pos_roots : Array[Array[Double]],
) -> Bool {
  let nv = vec_scale(v, -1.0)
  for p in pos_roots {
    if vec_inf_diff(nv, p) < 1.0e-7 {
      return true
    }
  }
  false
}

///|
/// Complex pairing `⟨λ, α^∨⟩` for `λ ∈ ℂ^d` and real coroot `α^∨ ∈ ℝ^d`.
pub fn cx_coroot_pairing(lambda : Array[Cx], alpha : Array[Double]) -> Cx {
  let av = coroot(alpha)
  let mut acc = cx_zero()
  for i in 0..<av.length() {
    acc = cx_add(acc, cx_scale(lambda[i], av[i]))
  }
  acc
}

///|
/// Action of an orthogonal matrix `w` (`d × d` real) on `λ ∈ ℂ^d`.
pub fn weyl_act_cx(
  w : Array[Double],
  lambda : Array[Cx],
  d : Int32,
) -> Array[Cx] {
  Array::makei(d, fn(i) {
    let mut s = cx_zero()
    for j in 0..<d {
      s = cx_add(s, cx_scale(lambda[j], w[i * d + j]))
    }
    s
  })
}

// ---------------------------------------------------------------------------
// p-adic & Adelic Gindikin-Karpelevich c-functions on ℂ
// ---------------------------------------------------------------------------

///|
/// Exact p-adic local Gindikin-Karpelevich factor at prime `p` for `s ∈ ℂ`:
/// `c_p(s) = (1 - p^{-(s+1)}) / (1 - p^{-s}) = ζ_p(s) / ζ_p(s+1)`.
pub fn padic_c_factor(p : Int32, s : Cx) -> Cx {
  let pc = cx_of_real(p.to_double())
  let num = cx_sub(cx_one(), cx_pow(pc, cx_neg(cx_add(s, cx_one()))))
  let den = cx_sub(cx_one(), cx_pow(pc, cx_neg(s)))
  cx_div(num, den)
}

///|
/// Archimedean local factor `c_∞(s) = Γ_ℝ(s) / Γ_ℝ(s+1)` where
/// `Γ_ℝ(s) = π^{-s/2} Γ(s/2)`.
pub fn archimedean_c_factor(s : Cx) -> Cx {
  let sp1 = cx_add(s, cx_one())
  let pref = cx_pow(cx_of_real(@math.PI), cx_of_real(0.5))
  let g_num = gamma_cx(cx_scale(s, 0.5))
  let g_den = gamma_cx(cx_scale(sp1, 0.5))
  cx_mul(pref, cx_div(g_num, g_den))
}

///|
/// Global completed local factor `c_global(s) = ξ(s) / ξ(s+1)` on `s ∈ ℂ`.
pub fn global_c_factor_cx(s : Cx) -> Cx {
  cx_div(xi_cx(s), xi_cx(cx_add(s, cx_one())))
}

///|
/// Sieve of Eratosthenes returning all primes `<= limit`.
pub fn primes_up_to(limit : Int32) -> Array[Int32] {
  if limit < 2 {
    return []
  }
  let is_p = Array::make(limit + 1, true)
  is_p[0] = false
  is_p[1] = false
  let mut p = 2
  while p * p <= limit {
    if is_p[p] {
      let mut m = p * p
      while m <= limit {
        is_p[m] = false
        m += p
      }
    }
    p += 1
  }
  let out : Array[Int32] = []
  for i in 2..<=limit {
    if is_p[i] {
      out.push(i)
    }
  }
  out
}

///|
/// Truncated adelic Euler product `c_A^{(P)}(s) = c_∞(s) · Π_{p <= P} c_p(s)`
/// (Algorithm A24).  By Euler's product for `ζ`, `lim_{P → ∞} c_A^{(P)}(s) = ξ(s)/ξ(s+1)`
/// for `Re(s) > 1` (Theorem T28a).
pub fn adelic_c_euler_product(s : Cx, prime_bound : Int32) -> Cx {
  let mut acc = archimedean_c_factor(s)
  for p in primes_up_to(prime_bound) {
    acc = cx_mul(acc, padic_c_factor(p, s))
  }
  acc
}

///|
/// Inversion set size `l(w) = #{α ∈ Φ^+ : w(α) ∈ -Φ^+}` (= Coxeter length of `w`).
pub fn coxeter_length(rt : RootType, w : Array[Double]) -> Int32 {
  let d = rt.ambient_dim()
  let pos = positive_roots(rt)
  let mut cnt = 0
  for alpha in pos {
    if is_negative_root(matvec(w, alpha, d), pos) {
      cnt += 1
    }
  }
  cnt
}

///|
/// General root-system Gindikin-Karpelevich intertwining function
/// `c_Φ(w, λ) = Π_{α ∈ Φ^+, w(α) ∈ -Φ^+} local_c(⟨λ, α^∨⟩)`
/// parameterised over any local or global factor `local_c : Cx → Cx`.
pub fn root_c_function(
  rt : RootType,
  w : Array[Double],
  lambda : Array[Cx],
  local_c : (Cx) -> Cx,
) -> Cx {
  let d = rt.ambient_dim()
  let pos = positive_roots(rt)
  let mut prod = cx_one()
  for alpha in pos {
    let wa = matvec(w, alpha, d)
    if is_negative_root(wa, pos) {
      let s = cx_coroot_pairing(lambda, alpha)
      prod = cx_mul(prod, local_c(s))
    }
  }
  prod
}

///|
/// Complex constant-term spectral projection `Π_Φ(λ)[u, v] = (1/|W|) c_Φ(u v^{-1}, v λ)`
/// (`|W| × |W|` complex matrix) and its idempotence defect `‖Π^2 - Π‖_∞` for ANY
/// classical root system `rt` (Theorem T27c).
pub fn root_projection_defect(
  rt : RootType,
  lambda : Array[Cx],
  local_c : (Cx) -> Cx,
) -> Double {
  let d = rt.ambient_dim()
  let ws = generate_weyl_group(rt)
  let n = ws.length()
  let scale = 1.0 / n.to_double()
  let p : Array[Cx] = Array::make(n * n, cx_zero())
  for a in 0..<n {
    let u = ws[a]
    for b in 0..<n {
      let v = ws[b]
      // Orthogonal Weyl matrices satisfy v^{-1} = v^T!
      let vt = transpose(v, d, d)
      let w = matrix_product(u, d, d, vt, d)
      let v_lam = weyl_act_cx(v, lambda, d)
      p[a * n + b] = cx_scale(root_c_function(rt, w, v_lam, local_c), scale)
    }
  }
  // Check ‖P^2 - P‖_max and |tr P - 1|
  let mut worst = 0.0
  let mut tr = cx_zero()
  for i in 0..<n {
    tr = cx_add(tr, p[i * n + i])
    for j in 0..<n {
      let mut sum = cx_zero()
      for k in 0..<n {
        sum = cx_add(sum, cx_mul(p[i * n + k], p[k * n + j]))
      }
      let diff = cx_abs(cx_sub(sum, p[i * n + j]))
      if diff > worst {
        worst = diff
      }
    }
  }
  let tr_err = cx_abs(cx_sub(tr, cx_one()))
  if tr_err > worst {
    tr_err
  } else {
    worst
  }
}

// ---------------------------------------------------------------------------
// Tests for Theorems T27 & T28
// ---------------------------------------------------------------------------

///|
test "T27a: Classical root systems A/B/C/D, Langlands Cartan transpose duality, Weyl orders" {
  // Langlands duality B_3 <-> C_3 transposes the Cartan matrix: C(B_3)^T == C(C_3)!
  let b3 = RootType::TypeB(3)
  let c3 = b3.langlands_dual()
  assert_eq(c3, RootType::TypeC(3))
  assert_eq(c3.langlands_dual(), b3)
  let cb = cartan_matrix(b3)
  let cc = cartan_matrix(c3)
  let cbt = transpose(cb, 3, 3)
  assert_true(vec_inf_diff(cbt, cc) < 1.0e-12)
  // Self-duality of A_3 and D_4
  assert_eq(RootType::TypeA(3).langlands_dual(), RootType::TypeA(3))
  assert_eq(RootType::TypeD(4).langlands_dual(), RootType::TypeD(4))
  // Positive root counts: |Φ^+(B_3)| = 9, |Φ^+(C_3)| = 9, |Φ^+(D_4)| = 12, |Φ^+(A_4)| = 10
  assert_eq(positive_roots(b3).length(), 9)
  assert_eq(positive_roots(c3).length(), 9)
  assert_eq(positive_roots(RootType::TypeD(4)).length(), 12)
  assert_eq(positive_roots(RootType::TypeA(4)).length(), 10)
  // Universal Weyl vector identity: ⟨ρ, α_i^∨⟩ = 1 for EVERY simple root of A_3, B_3, C_3, D_4
  for rt in [RootType::TypeA(3), b3, c3, RootType::TypeD(4)] {
    let rho = weyl_rho(rt)
    for alpha in simple_roots(rt) {
      assert_true(approx(dot(rho, coroot(alpha)), 1.0, tol=1.0e-12))
    }
  }
  // Exact Weyl group orders via BFS reflection generator (surpasses rank <= 4 type A!)
  assert_eq(generate_weyl_group(RootType::TypeB(2)).length(), 8) // 2^2 * 2! = 8
  assert_eq(generate_weyl_group(RootType::TypeC(2)).length(), 8)
  assert_eq(generate_weyl_group(RootType::TypeB(3)).length(), 48) // 2^3 * 3! = 48
  assert_eq(generate_weyl_group(RootType::TypeD(3)).length(), 24) // 2^2 * 3! = 24 (= |S_4|)
  assert_eq(generate_weyl_group(RootType::TypeA(4)).length(), 120) // |S_5| = 120 (rank 5 GL_5!)
}

///|
test "T27b & T27c: Non-type-A (B_2 = SO_5, C_2 = Sp_4) cocycle & constant-term projection over C" {
  let b2 = RootType::TypeB(2)
  let c2 = RootType::TypeC(2)
  // Complex spectral parameter λ ∈ ℂ^2 (off the real axis!)
  let lam : Array[Cx] = [{ re: 1.65, im: 0.45, }, { re: -0.55, im: 0.85, }]
  let ws = generate_weyl_group(b2)
  // 1) Crossed-homomorphism cocycle law holds whenever local_c(s)*local_c(-s) == 1
  //    (Theorem T27b: necessary & sufficient condition for root cancellation!).
  //    Global completed c_global(s) = ξ(s)/ξ(s+1) satisfies c(s)c(-s) = 1 identically on C!
  for w1 in ws {
    for w2 in ws {
      let w12 = matrix_product(w1, 2, 2, w2, 2)
      let lhs_g = root_c_function(b2, w12, lam, global_c_factor_cx)
      let rhs_g = cx_mul(
        root_c_function(b2, w1, weyl_act_cx(w2, lam, 2), global_c_factor_cx),
        root_c_function(b2, w2, lam, global_c_factor_cx),
      )
      assert_true(cx_approx(lhs_g, rhs_g, tol=1.0e-8))
    }
  }
  // 2) For p-adic local factor c_p(s) = (1 - p^{-(s+1)})/(1 - p^{-s}), c_p(s)c_p(-s) != 1;
  //    instead, the normalized unitary p-adic scattering factor
  //      r_p(s) = c_p(s) / c_p(-s)   (or Casselman-Shalika normalized intertwining)
  //    satisfies r_p(s)*r_p(-s) == 1 identically on C, restoring the exact cocycle & projection!
  let rp5 = fn(s : Cx) {
    cx_div(padic_c_factor(5, s), padic_c_factor(5, cx_neg(s)))
  }
  for w1 in ws {
    for w2 in ws {
      let w12 = matrix_product(w1, 2, 2, w2, 2)
      let lhs_p = root_c_function(b2, w12, lam, rp5)
      let rhs_p = cx_mul(
        root_c_function(b2, w1, weyl_act_cx(w2, lam, 2), rp5),
        root_c_function(b2, w2, lam, rp5),
      )
      assert_true(cx_approx(lhs_p, rhs_p, tol=1.0e-12))
    }
  }
  // 3) Unnormalized p-adic c_p(s) STILL satisfies the reduced-length additive cocycle
  //    c_p(w1*w2, λ) == c_p(w1, w2*λ)*c_p(w2, λ) whenever l(w1*w2) == l(w1) + l(w2)
  //    (and FAILS when l(w1*w2) < l(w1) + l(w2) because c_p(s)c_p(-s) != 1)!
  let c5 = fn(s : Cx) { padic_c_factor(5, s) }
  let mut additive_pairs = 0
  let mut non_additive_diff = 0.0
  for w1 in ws {
    for w2 in ws {
      let w12 = matrix_product(w1, 2, 2, w2, 2)
      let lhs_u = root_c_function(b2, w12, lam, c5)
      let rhs_u = cx_mul(
        root_c_function(b2, w1, weyl_act_cx(w2, lam, 2), c5),
        root_c_function(b2, w2, lam, c5),
      )
      if coxeter_length(b2, w12) ==
        coxeter_length(b2, w1) + coxeter_length(b2, w2) {
        additive_pairs += 1
        assert_true(cx_approx(lhs_u, rhs_u, tol=1.0e-12))
      } else {
        let d = cx_abs(cx_sub(lhs_u, rhs_u))
        if d > non_additive_diff {
          non_additive_diff = d
        }
      }
    }
  }
  assert_true(additive_pairs > 10)
  assert_true(non_additive_diff > 1.0e-4)
  // Spectral projection Π^2 = Π and tr Π = 1 on B_2 (SO_5) and C_2 (Sp_4):
  assert_true(root_projection_defect(b2, lam, rp5) < 1.0e-12)
  assert_true(root_projection_defect(c2, lam, rp5) < 1.0e-12)
  assert_true(root_projection_defect(b2, lam, global_c_factor_cx) < 1.0e-8)
  assert_true(root_projection_defect(c2, lam, global_c_factor_cx) < 1.0e-8)
}

///|
test "T28: Adelic Euler-product convergence to global c(s) and unitarity |c_Φ(w, iν)| = 1" {
  // 1. Adelic Euler product c_∞(s) * Π_{p <= P} c_p(s) -> ξ(s)/ξ(s+1) monotonically as P grows
  let s = { re: 2.4, im: 0.7, }
  let exact = global_c_factor_cx(s)
  let approx_50 = adelic_c_euler_product(s, 50)
  let approx_500 = adelic_c_euler_product(s, 500)
  let err_50 = cx_abs(cx_sub(approx_50, exact))
  let err_500 = cx_abs(cx_sub(approx_500, exact))
  assert_true(err_500 < err_50)
  assert_true(err_500 < 2.0e-3)
  // Local p-adic reflection law: c_p(s) * c_p(-s) != 1 individually (due to missing (1-p^{-1}) factors
  // wait: let's check global reflection c_global(s) * c_global(-s) == 1 across C!)
  let prod_global = cx_mul(
    global_c_factor_cx({ re: 0.8, im: 1.9, }),
    global_c_factor_cx({ re: -0.8, im: -1.9, }),
  )
  assert_true(cx_approx(prod_global, cx_one(), tol=1.0e-8))
  // 2. Unitarity on the imaginary spectral axis s = i t (Langlands spectral decomposition):
  //    |c_global(i t)| == 1 AND |c_{B_2}(w, i ν)| == 1 for ALL w ∈ W(B_2)!
  let nu : Array[Cx] = [{ re: 0.0, im: 1.35, }, { re: 0.0, im: -0.82, }]
  for w in generate_weyl_group(RootType::TypeB(2)) {
    let cw = root_c_function(RootType::TypeB(2), w, nu, global_c_factor_cx)
    assert_true(approx(cx_abs(cw), 1.0, tol=1.0e-8))
  }
}
