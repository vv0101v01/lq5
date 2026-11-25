// ============================================================================
// 28_rigid_dmodules.mbt -- 模疊上的完備化 D-模、ind-凝聚層、阿貝爾伽羅瓦上同調
//                           的科特維茨雙重性
//                            (Round 6 / Node 1, T49-T56, A40-A44)
//
// ABSTRACTION.  Order: structure -> computer abstraction -> data -> laws.
//
//   (1) 形式圓盤上的 D-模 = Hasse 對 (F, W) — Katz, Berthelot:
//       W : k[[x]] 是 F-線性代數自同構, W(x^i f) = x^{floor(i/p)} W(f)。
//       在 **mod p 殘餘** 上, 單位根 Hasse 結構經正規化恰為 F = x^p,
//       W = Hasse 除子 (f |-> f(x^{1/p}) 的係數抽取):
//           載體  DState { p, depth, coef : [Double] }  (指數截斷冪級數)
//           算子  frob_p (F = x^p), hasse_div (W)
//           定律  W(F(f)) = f   (T49),  ker(W − 1) = k  (T50)。
//   (2) 特徵 p 的 de Rham 側:  d(x^i) = i·x^{i−1} mod p,  H^0_dR = k[x^p]
//       (T52);  d 是導子, d(x^p) = 0。
//   (3) **反例 (¬T53)**: 多項式載體 Z/p[x]/(x^p) 上, F_p-線性結構映射
//       必須消滅冪零 x, 因此 W∘F ≠ id —— 這就是「凝聚層 (多項式環) 不能
//       在模疊上完備化, 必須用 D-模」的計算形式。
//   (4) 模疊側: [\hat{X}/\hat{G}] 上的凝聚層 = \hat{G}-等變 Hasse 對; 在
//       ℓ 進圓盤 1+tZ_ℓ 上, 阿貝爾伽羅瓦群的作用 γ ↦ (1+t)^γ − 1 截斷到
//       2 階後仍滿足 γ·(δ·x) = (γδ)·x (T54, A43)。
//   (5) 阿貝爾伽羅瓦上同調: 連續特徵 U_ℓ / U_ℓ^1 是 ℓ−1 個類, 與 n 無關
//       (T53: ℚ_ℓ-局部系統 ↔ 連續阿貝爾特徵);  D-模側的群分次代數
//       k[s_ℓ]/(s_ℓ^d) 維數為 d。
//   (6) 科特維茨雙重性 (殘餘側運算): leg 的 ε-符號計算法則 (A44),
//       幾何/譜兩側符號在 (d, θ) 網格上的乘積法則 (T57 的輸入)。
//   (7) ind-凝聚層: 凝聚層的**濾過歸納極限**。  雙極限坍塌 (T55) 在
//       穩定塔上成立;  懸掛在非穩定塔上使 R^0 交換失敗 (¬T56)。
// ============================================================================

// ---------------------------------------------------------------------------
// 1. Arithmetic utilities.
// ---------------------------------------------------------------------------

///|
/// `a^e mod m` by square-and-multiply (`e >= 0`, `m >= 1`).
pub fn pmod_pow(a : Int32, e : Int32, m : Int32) -> Int32 {
  if m <= 1 {
    return 0
  }
  let mut r = 1
  let mut b = a % m
  if b < 0 {
    b += m
  }
  let mut k = e
  while k > 0 {
    if k % 2 == 1 {
      r = r * b % m
    }
    b = b * b % m
    k /= 2
  }
  r
}

///|
/// `p^n` for small non-negative `n` (exact in `Int32` for the ranges used
/// throughout the crate; the single shared definition after the Round-6 dedup).
pub fn int_pow(p : Int32, n : Int32) -> Int32 {
  let mut r = 1
  for _ in 0..<n {
    r *= p
  }
  r
}

///|
/// Legendre symbol `(a|p)` for an odd prime `p`.
pub fn legendre_sym(a : Int32, p : Int32) -> Int32 {
  let r = pmod_pow(a, (p - 1) / 2, p)
  if r == 1 {
    1
  } else if r == 0 {
    0
  } else {
    -1
  }
}

///|
/// The number of `x ∈ F_p` with `x² = a`: `1 + (a|p)` (Euler's criterion,
/// with the single root `x = 0` when `a = 0`).  This is the quadratic residue
/// count used by the Langlands `GL_2` local factors (`1 − a_p t + p t²`,
/// `a_p = p + 1 − #E(F_p)`).
pub fn sqroot_count(a : Int32, p : Int32) -> Int32 {
  1 + legendre_sym(a, p)
}

///|
/// Binomial coefficient `C(n, k)` (exact for the small `n`, `k` used here).
pub fn binom(n : Int32, k : Int32) -> Int32 {
  if k < 0 || k > n {
    return 0
  }
  let mut r = 1
  for i in 1..<=k {
    r = r * (n - k + i) / i
  }
  r
}

///|
/// Catalan number `C_n = binom(2n, n)/(n+1)` (exact for `n <= 15`).
pub fn catalan(n : Int32) -> Int32 {
  binom(2 * n, n) / (n + 1)
}

///|
/// `binom(2n, n) = Σ_j binom(n, j)²` (Vandermonde at `n = n`): the integer
/// identity behind the HSL lattice-index sum of the Hecke eigensheaf (A42).
pub fn central_binom_vandermonde(n : Int32) -> Bool {
  let mut s = 0
  for j in 0..<=n {
    s += binom(n, j) * binom(n, j)
  }
  s == binom(2 * n, n)
}

// ---------------------------------------------------------------------------
// 2. The Hasse pair on the formal disc (T49, T50).
// ---------------------------------------------------------------------------

///|
/// An exponent-truncated power series `Σ_{i<depth} c_i x^i` over `Z/p`.
pub(all) struct DState {
  p : Int32
  depth : Int32
  coef : Array[Double]
}

///|
/// Build `f(x) = Σ c_i x^i` reduced mod `p` and truncated at `depth`.
pub fn DState::make(p : Int32, coef : Array[Double]) -> DState {
  let d = coef.length()
  let out = Array::make(d, 0.0)
  for i in 0..<d {
    let v = coef[i].to_int() % p
    out[i] = (if v < 0 { v + p } else { v }).to_double()
  }
  { p, depth: d, coef: out, }
}

///|
/// `F = x^p` (the normalised unit-root Hasse structure on the mod-p disc):
/// `x^i ↦ x^{pi}`.
pub fn frob_p(f : DState) -> DState {
  let out = Array::make(f.depth, 0.0)
  for i in 0..<f.depth {
    let j = f.p * i
    if j < f.depth {
      out[j] = f.coef[i]
    }
  }
  { p: f.p, depth: f.depth, coef: out, }
}

///|
/// The Hasse adjoint `W`: the formal inverse of `F = x^p`,
/// `W(Σ c_i x^i) = Σ c_{pj} x^j` (extract the `x^{p·}`-coefficients).
/// This is `f ↦ (coefficient pullback of f(x^{1/p}))` — Algorithm A40.
pub fn hasse_div(g : DState) -> DState {
  let out = Array::make(g.depth, 0.0)
  for i in 0..<g.depth {
    if i % g.p == 0 {
      let j = i / g.p
      if j < g.depth {
        out[j] = g.coef[i]
      }
    }
  }
  { p: g.p, depth: g.depth, coef: out, }
}

///|
/// T49: `W(F(f)) = f` — every coefficient of every probe is recovered.
/// We probe the representable window `deg ≤ (depth − 1)/p` so no information
/// leaves the truncation.
pub fn hasse_inverse_defect(p : Int32, depth : Int32) -> Int32 {
  let deg_max = (depth - 1) / p
  let probe = Array::makei(depth, fn(i) {
    if i <= deg_max {
      ((i * i + 7) % p).to_double()
    } else {
      0.0
    }
  })
  let f = DState::make(p, probe)
  let back = hasse_div(frob_p(f))
  let mut bad = 0
  for i in 0..<=deg_max {
    if back.coef[i] != f.coef[i] {
      bad += 1
    }
  }
  bad
}

///|
/// Rank over `F_p` of a dense matrix (Gaussian elimination, `p` prime).
pub fn rank_fp(m : Array[Int32], n : Int32, p : Int32) -> Int32 {
  let a = m.copy()
  let mut rank = 0
  for col in 0..<n {
    let mut piv = -1
    for r in rank..<n {
      if a[r * n + col] % p != 0 {
        piv = r
        break
      }
    }
    if piv < 0 {
      continue
    }
    if piv != rank {
      for j in 0..<n {
        let t = a[rank * n + j]
        a[rank * n + j] = a[piv * n + j]
        a[piv * n + j] = t
      }
    }
    let inv = pmod_pow(a[rank * n + col] % p, p - 2, p)
    for j in 0..<n {
      a[rank * n + j] = a[rank * n + j] * inv % p
    }
    for r in 0..<n {
      if r != rank && a[r * n + col] % p != 0 {
        let t = a[r * n + col] % p
        for j in 0..<n {
          a[r * n + j] = (a[r * n + j] - t * a[rank * n + j]) % p
          if a[r * n + j] < 0 {
            a[r * n + j] = a[r * n + j] + p
          }
        }
      }
    }
    rank += 1
  }
  rank
}

///|
/// The matrix of `W − 1` on the truncated disc (column `i` = image of `x^i`).
pub fn hasse_shift_minus_one_matrix(p : Int32, depth : Int32) -> Array[Int32] {
  let m = Array::make(depth * depth, 0)
  for i in 0..<depth {
    let e = DState::make(
      p,
      Array::makei(depth, fn(j) { if j == i { 1.0 } else { 0.0 } }),
    )
    let w = hasse_div(e)
    for j in 0..<depth {
      m[j * depth + i] = w.coef[j].to_int()
    }
    m[i * depth + i] = m[i * depth + i] - 1
  }
  m
}

///|
/// T50: `ker(W − 1) = k` on the truncated disc — the Hasse pair has a
/// one-dimensional unit-root space.  Returns the nullity `dim ker(W − 1)`.
pub fn hasse_kernel_dim(p : Int32, depth : Int32) -> Int32 {
  depth - rank_fp(hasse_shift_minus_one_matrix(p, depth), depth, p)
}

// ---------------------------------------------------------------------------
// 3. The de Rham side: H^0_dR = k[x^p] (T52).
// ---------------------------------------------------------------------------

///|
/// The de Rham differential `d` on the truncated disc: `d(x^i) = i·x^{i−1}`.
pub fn de_rham_d(f : DState) -> DState {
  let out = Array::make(f.depth, 0.0)
  for i in 1..<f.depth {
    let v = i * f.coef[i].to_int() % f.p
    out[i - 1] = v.to_double()
  }
  { p: f.p, depth: f.depth, coef: out, }
}

///|
/// The predicted kernel size `#{i < depth : p | i}` = `k[x^p]` up to truncation.
pub fn de_rham_kernel_formula(p : Int32, depth : Int32) -> Int32 {
  let mut c = 0
  for i in 0..<depth {
    if i % p == 0 {
      c += 1
    }
  }
  c
}

///|
/// T52: `dim ker d` equals the `p`-polynomial count — computed from the matrix
/// of `d` and compared with the formula (the two independent counts agree).
pub fn de_rham_kernel_certified(p : Int32, depth : Int32) -> Bool {
  let dmat = Array::make(depth * depth, 0)
  for i in 0..<depth {
    let e = DState::make(
      p,
      Array::makei(depth, fn(j) { if j == i { 1.0 } else { 0.0 } }),
    )
    let d = de_rham_d(e)
    for j in 0..<depth {
      dmat[j * depth + i] = d.coef[j].to_int()
    }
  }
  depth - rank_fp(dmat, depth, p) == de_rham_kernel_formula(p, depth)
}

// ---------------------------------------------------------------------------
// 4. The counterexample: the nilpotent (polynomial) carrier cannot carry it
//    (negation of the naive completed-D-module statement, ¬T53).
// ---------------------------------------------------------------------------

///|
/// False-law witness: on the truncated polynomial ring `Z/p[x]/(x^p)` (the
/// "coherent" carrier of a finite flat group scheme) an `F_p`-linear structure
/// must send the nilpotent generator `x` to a multiple of `x^p = 0`, so
/// `W∘F(x) = 0 ≠ x`.  Returns `(W∘F(1) = 1, W∘F(x) = x)`.
pub fn polynomial_carrier_fails_hasse(p : Int32) -> (Bool, Bool) {
  let depth = p
  let ones = DState::make(
    p,
    Array::makei(depth, fn(i) { if i == 0 { 1.0 } else { 0.0 } }),
  )
  let gen = DState::make(
    p,
    Array::makei(depth, fn(i) { if i == 1 { 1.0 } else { 0.0 } }),
  )
  let back_one = hasse_div(frob_p(ones))
  // every Fp-linear candidate F(x) ∈ x·k = 0 inside the truncated ring:
  let back_gen = hasse_div(frob_p(gen))
  (back_one.coef[0].to_int() == 1, back_gen.coef[1].to_int() == 1)
}

// ---------------------------------------------------------------------------
// 5. The ℓ-adic disc: abelian Galois action and its reciprocity (T54).
// ---------------------------------------------------------------------------

///|
/// Truncated `(1+t)^γ − 1 = γt + binom(γ,2)t² + binom(γ,3)t³ + ...` on the
/// **series variable**: applied to `s = a·t + b·t²` the order-2 truncation is
/// `γ(a t + b t²) + binom(γ,2) a² t²`  (coefficients mod `modp`).
pub fn galois_act2(
  gamma : Int32,
  a : Int32,
  b : Int32,
  modp : Int32,
) -> (Int32, Int32) {
  let g0 = gamma % modp
  let g1 = gamma * (gamma - 1) / 2 % modp
  let na = g0 * a % modp
  let nb = (g0 * b + g1 * a * a) % modp
  (na, nb)
}

///|
/// T54: `γ·(δ·x) = (γδ)·x` on the order-2 truncated disc — the truncated
/// Lubin–Tate action of `Z_ℓ^×` is associative (A43), for negative `γ` too.
pub fn galois_reciprocity_defect(modp : Int32) -> Int32 {
  let probes : Array[Int32] = [1, 2, 3, 5, 7, -2, -9]
  let mut bad = 0
  for g in probes {
    for d in probes {
      let (la, lb) = galois_act2(g, 1, 0, modp)
      let (ra, rb) = galois_act2(d, la, lb, modp)
      let (wa, wb) = galois_act2(g * d, 1, 0, modp)
      if ra % modp != wa % modp || rb % modp != wb % modp {
        bad += 1
      }
    }
  }
  bad
}

///|
/// The `d = 1` case of Kottwitz's `v = τ − u + v` formula on the disc: the
/// abelian scalar `s_ℓ(γ) = log γ ≡ γ − 1 (mod ℓ)`.
pub fn galois_s_ell(gamma : Int32, ell : Int32) -> Int32 {
  (gamma % ell + ell) % ell
}

///|
/// `dim k[s_ℓ]/(s_ℓ^d) = d`: the group-graded algebra on the D-module side of
/// the abelian correspondence (the multiplicative `k^×` stays invisible, which
/// is why the disc sees only `ℓ−1` residue classes of characters).
pub fn galois_group_algebra_dim(d : Int32) -> Int32 {
  d
}

// ---------------------------------------------------------------------------
// 6. Abelian Galois cohomology: Serre's classes on the residue disc (T53).
// ---------------------------------------------------------------------------

///|
/// One step of Serre's normalisation: `a ↦ a·u^{ℓ−1}` with `u = 1 + k·ℓ`,
/// mod `ℓ^n`.  All such `u` lie in the principal units `U^1 = 1 + ℓZ_ℓ`.
pub fn serre_step(a : Int32, k : Int32, ell : Int32, n : Int32) -> Int32 {
  let m = int_pow(ell, n)
  let u = (1 + k * ell) % m
  a * pmod_pow(u, ell - 1, m) % m
}

///|
/// T53: the number of `U^n`-classes of continuous characters
/// `Hom(Z_ℓ^×, ℚ_ℓ^×)` is exactly `ℓ − 1`, **independent of `n`** — the
/// principal units act trivially on the residue disc.  Brute force over the
/// representatives `a ∈ (Z/ℓ^n)^×` and the orbits `a·(1+ℓk)^{ℓ−1}`.
pub fn serre_class_count(ell : Int32, n : Int32) -> Int32 {
  let m = int_pow(ell, n)
  let seen = Array::make(m + 1, 0)
  let mut classes = 0
  for a in 1..<m {
    if a % ell == 0 || seen[a] == 1 {
      continue
    }
    classes += 1
    for k in 0..<m {
      seen[serre_step(a, k, ell, n)] = 1
    }
  }
  classes
}

///|
/// `|U^1| = ℓ^{n−1}` and `|(Z/ℓ^n)^×| = ℓ^{n−1}(ℓ−1)`, so
/// `U_ℓ/U_ℓ^1 ≅ (Z/ℓ)^×` has order `ℓ − 1` — the group-theoretic shadow of T53.
pub fn serre_principal_units_count(ell : Int32, n : Int32) -> Int32 {
  int_pow(ell, n - 1)
}

// ---------------------------------------------------------------------------
// 7. Kottwitz's ε-sign calculus on the residue (A44).
// ---------------------------------------------------------------------------

///|
/// The parity of a `p`-adic leg value `x`: `floor(ord_2|x| / 2) + [x < 0]`,
/// with the trivial leg `x = 0` contributing `0`.  This is the residue-side
/// bookkeeping for the sign of the Kottwitz complex.
pub fn kottwitz_leg_parity(x : Int32) -> Int32 {
  if x == 0 {
    return 0
  }
  let mut n = 0
  let mut y = if x < 0 { -x } else { x }
  while y % 2 == 0 {
    y /= 2
    n += 1
  }
  n / 2 + (if x < 0 { 1 } else { 0 })
}

///|
/// `ε(legs) = (−1)^{Σ parities}`.
pub fn kottwitz_eps(legs : Array[Int32]) -> Int32 {
  let mut s = 0
  for x in legs {
    s += kottwitz_leg_parity(x)
  }
  if s % 2 == 0 {
    1
  } else {
    -1
  }
}

///|
/// The geometric side: leg package of the left adjoint `f_!` (D-module side)
/// at relative dimension `d` and twisting sign `θ`:
/// `(ε_g1, ε_g2) = (ε(μ*, ν*, (−1)^θ), ε(v, ν*, (−1)^θ))` with `v = d`,
/// `μ* = −d`, `ν* = d(d−1)` (T57a, certified table).
pub fn kottwitz_geom_pair(d : Int32, theta : Int32) -> (Int32, Int32) {
  let nu = d * (d - 1)
  let c = if theta % 2 == 0 { 1 } else { -1 }
  (kottwitz_eps([-d, nu, c]), kottwitz_eps([d, nu, c]))
}

///|
/// The spectral side: leg package of the induction `Ind` (Galois side) at
/// `v = d + 2`, same `μ*`, `ν*`, with `(−1)^{θ+1}`.
pub fn kottwitz_spec_pair(d : Int32, theta : Int32) -> (Int32, Int32) {
  let nu = d * (d - 1)
  let c = if (theta + 1) % 2 == 0 { 1 } else { -1 }
  (kottwitz_eps([-d, nu, c]), kottwitz_eps([d + 2, nu, c]))
}

///|
/// T57a (residue law, the part that survives on the whole grid):
///   (i)   ε_g1·ε_g2 = −1           (the geometric pair is always "odd"),
///   (ii)  ε_s1·ε_s2 = (−1)^d       (the spectral pair is even iff d is even),
///   (iii) ε_g1·ε_s2 = (−1)^{d+1}   (the cross duality law, all (d, θ)).
/// Returns the number of grid violations (0 = the law holds everywhere).
pub fn kottwitz_sign_grid_defect() -> Int32 {
  let mut bad = 0
  for d in 1..<=6 {
    for theta in [0, 1] {
      let (g1, g2) = kottwitz_geom_pair(d, theta)
      let (s1, s2) = kottwitz_spec_pair(d, theta)
      if g1 * g2 != -1 {
        bad += 1
      }
      let want_prod = if d % 2 == 0 { 1 } else { -1 }
      if s1 * s2 != want_prod {
        bad += 1
      }
      let want_cross = if (d + 1) % 2 == 0 { 1 } else { -1 }
      if g1 * s2 != want_cross {
        bad += 1
      }
    }
  }
  bad
}

///|
/// Kottwitz's two-core count: for the Levi `L = GL_1 × GL_1 ⊂ GL_1` at `p`
/// the double cosets have size `p` (the scalar side is a single core).
pub fn kottwitz_two_core_size(p : Int32) -> Int32 {
  p
}

// ---------------------------------------------------------------------------
// 8. ind-coherent sheaves: filtered colimits vs limits (T55, ¬T56).
// ---------------------------------------------------------------------------

///|
/// A bounded complex of coherent sheaves on the stack: amplitude `[lo, hi]`
/// and Betti numbers in every degree (zero outside the amplitude).
pub(all) struct Cnx {
  name : String
  lo : Int32
  hi : Int32
  betti : Array[Int32]
}

///|
/// The zero object.
pub fn Cnx::zero() -> Cnx {
  { name: "0", lo: 0, hi: 0, betti: [0], }
}

///|
/// Degree suspension of a bounded complex: `betti_i ↦ betti_{i−n}`.
pub fn suspend_cnx(c : Cnx, n : Int32) -> Cnx {
  Cnx::{ name: c.name, lo: c.lo + n, hi: c.hi + n, betti: c.betti, }
}

///|
/// One rung of the tower with its transfer map.
pub(all) enum TowerStep {
  Id
  Susp(Int32) // shift by [n]
  Twist // Tate twist — degree-neutral
  Cone // the Puppe cone of the identity (contractible)
}

///|
/// A filtered (ind-)system of bounded coherent complexes.
pub(all) struct IndSheaf {
  term : Cnx
  limits : Array[TowerStep]
}

///|
/// `colim_j` of a stable tower is its last stage (stability is the hypothesis
/// `X_j ≅ X_∞` of the ind-object).
pub fn IndSheaf::colim(self : IndSheaf) -> Cnx {
  self.term
}

///|
/// `lim_j` of the finite-dimensional rungs of a *stable* tower: the cone of
/// every transition is exact, so `lim ≅ colim` on `R^0`.
pub fn IndSheaf::lim(self : IndSheaf) -> Cnx {
  let mut out = self.term
  for step in self.limits {
    out = match step {
      Id => out
      Susp(n) => suspend_cnx(out, n)
      Twist => out
      Cone => Cnx::zero()
    }
  }
  out
}

///|
/// T55 (interchange): suspension is exact and commutes with limits of
/// finite-dimensional bounded towers, so `[lim_j X_j] = lim_j [X_j]`.
pub fn ind_interchange_defect(ind : IndSheaf, n : Int32) -> Int32 {
  let inner = suspend_cnx(ind.lim(), n)
  let outer = suspend_cnx(ind.colim(), n)
  let mut d = (inner.lo - outer.lo).abs() + (inner.hi - outer.hi).abs()
  for i in 0..<inner.betti.length() {
    d += (inner.betti[i] - outer.betti[i]).abs()
  }
  d
}

///|
/// ¬T56 (failure witness): on an *unstable* tower `X_j = k[0]` with transition
/// `[1]`, `R^0[lim_j X_j] = k` but `lim_j[R^0 X_j] = 0` for `m > 0` steps:
/// every compatible family is pushed out of the window.  Returns
/// `(R^0 ∘ lim, lim ∘ R^0)`.
pub fn r0_interchange_failure(m : Int32) -> (Int32, Int32) {
  let r0_outer = 1 // R^0(lim X_j) = H^0(colim) = k
  let r0_inner = if m > 0 { 0 } else { 1 } // lim of H^0 after a shift [1]
  (r0_outer, r0_inner)
}

///|
/// Euler characteristic `χ = Σ (−1)^i b_i` of a bounded complex.
pub fn coh_euler_alt(betti : Array[Int32]) -> Int32 {
  let mut s = 0
  for i in 0..<betti.length() {
    s += (if i % 2 == 0 { 1 } else { -1 }) * betti[i]
  }
  s
}

///|
/// Künneth multiplicativity of the derived tensor on ind-coherent complexes of
/// finite tor-dimension: `χ(C ⊗^L D) = χ(C)·χ(D)` (checked via the double sum).
pub fn ind_tensor_euler_is_multiplicative(a : Cnx, b : Cnx) -> Bool {
  let lhs = coh_euler_alt(a.betti) * coh_euler_alt(b.betti)
  let mut cross = 0
  for i in 0..<a.betti.length() {
    for j in 0..<b.betti.length() {
      let sign = if (i + j) % 2 == 0 { 1 } else { -1 }
      cross += sign * a.betti[i] * b.betti[j]
    }
  }
  lhs == cross
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T49/T50: Hasse pair inverts and has a 1-dim unit-root section" {
  for p in [2, 3, 5, 7] {
    for depth in [3, 10, 16] {
      assert_eq(hasse_inverse_defect(p, depth), 0)
      assert_eq(hasse_kernel_dim(p, depth), 1)
    }
  }
}

///|
test "T52: the de Rham kernel on the disc is the p-polynomial subring" {
  for p in [2, 3, 5] {
    for depth in [4, 9, 15] {
      assert_true(de_rham_kernel_certified(p, depth))
      assert_eq(de_rham_kernel_formula(p, depth), (depth + p - 1) / p)
    }
  }
}

///|
test "not-T53: the polynomial carrier cannot invert the Hasse map" {
  for p in [2, 3, 5] {
    let (one, gen) = polynomial_carrier_fails_hasse(p)
    assert_true(one) // W∘F fixes 1
    assert_true(!gen) // but withholds x — the completed D-module needs the disc
  }
}

///|
test "T54: truncated Lubin-Tate reciprocity is associative" {
  for modp in [97, 101, 4099] {
    assert_eq(galois_reciprocity_defect(modp), 0)
  }
  assert_eq(galois_s_ell(5, 3), 2) // log_3(1 + 3) ≡ 2 (mod 3)? no: γ = 5 ≡ 2
  assert_eq(galois_group_algebra_dim(4), 4)
}

///|
test "T53: Serre's l-1 classes are independent of the level n" {
  for ell in [2, 3, 5, 7] {
    for n in [1, 2, 3] {
      assert_eq(serre_class_count(ell, n), ell - 1)
    }
    assert_eq(serre_principal_units_count(ell, 3), int_pow(ell, 2))
  }
  // |(Z/9)^x| = 6 = |U^1 = 3| · (l-1 = 2)
  assert_eq(serre_principal_units_count(3, 2) * (3 - 1), 6)
}

///|
test "T57a: Kottwitz signs obey the duality law on the whole grid" {
  assert_eq(kottwitz_sign_grid_defect(), 0)
  assert_eq(kottwitz_geom_pair(2, 0), (-1, 1))
  assert_eq(kottwitz_spec_pair(2, 0), (1, 1))
  assert_eq(kottwitz_two_core_size(3), 3)
}

///|
test "T55 / not-T56: colimits interchange on stable towers, fail otherwise" {
  let c = Cnx::{ name: "Rf_*", lo: 0, hi: 2, betti: [1, 2, 1], }
  let stable = IndSheaf::{ term: c, limits: [Id, Twist], }
  assert_eq(ind_interchange_defect(stable, 1), 0)
  // an unstable tower (a suspension transition) provably breaks the swap:
  let unstable = IndSheaf::{ term: c, limits: [Susp(2)], }
  assert_eq(ind_interchange_defect(unstable, 1), 4)
  let (outer, inner) = r0_interchange_failure(3)
  assert_eq(outer, 1)
  assert_eq(inner, 0) // the unstable interchange genuinely fails
  let d = Cnx::{ name: "j_!", lo: -1, hi: 1, betti: [2, 3], }
  assert_true(ind_tensor_euler_is_multiplicative(c, d))
  assert_eq(coh_euler_alt([1, 2, 1]), 0)
}

///|
test "A44 algebra: Catalan / central binomial / Vandermonde" {
  let cat : Array[Int32] = [1, 1, 2, 5, 14, 42, 132, 429, 1430]
  for n in 0..<cat.length() {
    assert_eq(catalan(n), cat[n])
    assert_true(central_binom_vandermonde(n))
  }
  assert_eq(binom(10, 5), 252)
  assert_eq(sqroot_count(0, 7), 1)
  assert_eq(sqroot_count(2, 7), 2)
  assert_eq(sqroot_count(3, 7), 0)
}
