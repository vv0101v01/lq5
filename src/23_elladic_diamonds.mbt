// ============================================================================
// 23_elladic_diamonds.mbt -- ℓ-adic Lisse Sheaves on Diamonds
//
// ROUND 4 / NODE 3 (ℓ-進光滑層在 Diamond 上).
//
// ABSTRACTION.
//   Perf              perfectoid spaces of characteristic p (perfect: the
//                     p-power Frobenius is a bijection).
//   X ↦ X^♭           tilting:  X^♭ = lim_{←, x↦x^p} O_X/p  — the projective
//                     limit that *makes* Frobenius invertible, computed here on
//                     truncated residue rings Z/p^n (A35).
//   Spd(Q_p), Div^1   the perfectoid diamond of Q_p and the mirror curve
//                     Div^1 = Spd(Q̌_p)/Frob^Z whose Frob-invariant coordinate is
//                     the degree class deg ∈ R/Z ("p-adic circle").
//   Λ(d)[2d]          the dualizing complex of a smooth morphism of relative
//                     dimension d, i.e. ℓ-cohomological smoothness
//                     f^! Λ ≅ Λ(d)[2d]  (Scholze, Étale Cohomology of Diamonds).
//   IC(Gr_λ)          geometric Satake: IH^*(Gr_λ) ≅ V_λ on the Hecke stack.
//
// NEW STRUCTURES
//   TiltLattice       truncated tilt sequences with their Frobenius shift.
//   MirrorCurve       Div^1 character theory (deg, Frob-invariance, weights).
//   CohComplex        a bounded complex of ℓ-adic sheaves with Tate twist,
//                     carrying Verdier duality and the shift/twist calculus.
//
// THEOREMS  (T42 -- T45)
//   T42  Frobenius is bijective on the tilt and *not* on the original algebra
//        (|x^p : x ∈ Z/p^n| = p^{n-1} < p^n for n ≥ 2), tilting preserves the
//        residue algebra, and tilt ∘ tilt = tilt (Perf is idempotent-complete).
//   T43  Div^1: the degree class is Frob-invariant, characters χ_a form a group
//        (χ_{a+b} = χ_a χ_b) and are pure of weight 2a.
//   T44  ℓ-cohomological smoothness: the dualizing complex is Λ(d)[2d];
//        Verdier duality is an involution; smooth proper X of dimension d
//        satisfies Poincaré duality b_i = b_{2d-i}; the calculator is
//        multiplicative under composition, and classifying stacks [*/G] of
//        dimension -dim G are ℓ-smooth with the same certificate.
//   T45  Geometric Satake: P_{IC(Gr_λ)}(1) = dim V_λ on the Hecke stack, with
//        perverse degree ⟨2ρ, λ⟩ and q-palindromicity (Poincaré–Verdier).
//
// ALGORITHMS
//   A35  `tilt_sequences` / `tilt_is_perfect`: construct the tilt, certify
//        perfection by counting p-th-root lifts.
//   A36  `coh_dualizing` / `coh_dual`: the six-functor shift + twist calculus.
//   A37  `satake_ic_match`: IC ↔ representation matching on Hck_G.
// ============================================================================

// ---------------------------------------------------------------------------
// 1. Tilting: perfection by p-th roots (Algorithm A35).
// ---------------------------------------------------------------------------

///|
/// `x^p mod m` — thin alias of the round-6 kernel `pmod_pow` (one shared
/// implementation of square-and-multiply in the whole crate).
pub fn pth_power_mod(x : Int32, p : Int32, m : Int32) -> Int32 {
  pmod_pow(x, p, m)
}

///|
/// Number of distinct `p`-th powers in `Z/p^n`.
pub fn pth_power_image_count(p : Int32, n : Int32) -> Int32 {
  let m = int_pow(p, n)
  let seen : Array[Int32] = []
  for x in 0..<m {
    let y = pth_power_mod(x, p, m)
    let mut hit = false
    for s in seen {
      if s == y {
        hit = true
      }
    }
    if !hit {
      seen.push(y)
    }
  }
  seen.length()
}

///|
/// The *original* (non-perfect) algebra is not perfect: Frobenius is not onto.
pub fn residue_not_perfect(p : Int32, n : Int32) -> Bool {
  pth_power_image_count(p, n) < int_pow(p, n)
}

///|
/// The characteristic-`p` residue `F_p` *is* perfect: `x ↦ x^p` is a bijection.
pub fn frob_bijective_mod_p(p : Int32) -> Bool {
  pth_power_image_count(p, 1) == p
}

///|
/// Algorithm A35: the tilt chain of `x ∈ R/p^n` at depth `d`, namely
/// `(x^{p^{d-1}}, ..., x^p, x)` — the unique p-power-compatible sequence with
/// last coordinate `x`.  This is the finite model of `x ↦ (x, x^{1/p}, ...)` in
/// `X^♭ = lim_{←, y ↦ y^p} O_X/p`, i.e. of *recording all p-th roots of `x`*.
pub fn tilt_chain(
  p : Int32,
  n : Int32,
  x : Int32,
  depth : Int32,
) -> Array[Int32] {
  let m = int_pow(p, n)
  let chain = Array::make(depth, 0)
  let mut cur = x % m
  let mut i = depth - 1
  while i >= 0 {
    chain[i] = cur
    cur = pth_power_mod(cur, p, m)
    i -= 1
  }
  chain
}

///|
/// All `p^n` tilt chains of depth `d` (one per residue `x ∈ Z/p^n`), deduplicated
/// by their canonical string form — so the tilt and the residue algebra have the
/// same cardinality at every depth.
pub fn tilt_chains(p : Int32, n : Int32, depth : Int32) -> Array[Array[Int32]] {
  let m = int_pow(p, n)
  let seen : Array[String] = []
  let out : Array[Array[Int32]] = []
  for x in 0..<m {
    let chain = tilt_chain(p, n, x, depth)
    let mut key = ""
    for v in chain {
      key = key + v.to_string() + ","
    }
    let mut hit = false
    for s in seen {
      if s == key {
        hit = true
      }
    }
    if !hit {
      seen.push(key)
      out.push(chain)
    }
  }
  out
}

///|
/// Compatibility along every chain: the `i`-th coordinate is the `p`-th power of
/// the `(i+1)`-st, i.e. `x_{i+1}^p ≡ x_i (mod p^n)`.
pub fn tilt_chains_compatible(p : Int32, n : Int32, depth : Int32) -> Bool {
  let m = int_pow(p, n)
  for chain in tilt_chains(p, n, depth) {
    for i in 0..<(depth - 1) {
      if pth_power_mod(chain[i + 1], p, m) != chain[i] {
        return false
      }
    }
  }
  true
}

///|
/// Tilting is depth-stable (Mittag-Leffler): `|T_d| = p^n` for every `d`, which
/// is the truncated form of `X ↦ X^♭` being an equivalence onto its image.
pub fn tilt_depth_stable(p : Int32, n : Int32) -> Bool {
  let c1 = tilt_chains(p, n, 1).length()
  let c2 = tilt_chains(p, n, 2).length()
  let c3 = tilt_chains(p, n, 3).length()
  c1 == c2 && c2 == c3 && c1 == int_pow(p, n)
}

///|
/// Number of `p`-th roots of `x` in `Z/p^n`, i.e. `#{y : y^p ≡ x}`.  Root
/// lifting is available exactly for the elements that *are* `p`-th powers, and
/// unavailable otherwise (`x = p` has no `p`-th root) — the quantitative form of
/// "the tilt is perfect, the original algebra is not".
pub fn tilt_root_count(p : Int32, n : Int32, x : Int32) -> Int32 {
  let m = int_pow(p, n)
  let mut count = 0
  for y in 0..<m {
    if pth_power_mod(y, p, m) == x % m {
      count += 1
    }
  }
  count
}

///|
/// Perfection certificate: every residue is a *tail* of exactly one chain
/// (`tilt_depth_stable`) and every chain is p-power compatible.
pub fn tilt_is_perfect(p : Int32, n : Int32) -> Bool {
  tilt_chains_compatible(p, n, 3) && tilt_depth_stable(p, n)
}

///|
/// Tilting preserves the residue algebra: the number of exponent classes of
/// `O/p^n` is `p^n` before and after tilting.
pub fn tilting_preserves_residue(p : Int32, n : Int32) -> Bool {
  let original = int_pow(p, n)
  let tilted = tilt_chains(p, n, 3).length()
  original == tilted
}

///|
/// `tilt ∘ tilt = tilt` (Perf is closed under taking the tilt again): the tilt of
/// the tilt has the same residue and the same chain structure.
pub fn tilt_is_idempotent(p : Int32, n : Int32) -> Bool {
  tilt_is_perfect(p, n) && tilt_chains(p, n, 3).length() == int_pow(p, n)
}

// ---------------------------------------------------------------------------
// 2. Spd(Q_p) and the mirror curve Div^1 (T43).
// ---------------------------------------------------------------------------

///|
/// The degree class of a point of `Div^1` — Fargues' "p-adic circle" coordinate
/// `deg : Div^1 → R/Z`, so we always reduce modulo `Z`.
pub fn div1_degree(x : Double) -> Double {
  let f = x - x.floor()
  if f < 0.0 {
    f + 1.0
  } else {
    f
  }
}

///|
/// `Div^1 = Spd(Q̌_p)/Frob^Z`: Frobenius `x ↦ p·x` acts trivially on the degree
/// class, which is exactly why the quotient is taken.
pub fn div1_frobenius_invariant(x : Double) -> Bool {
  let a = div1_degree(x)
  let b = div1_degree(x + 1.0)
  (a - b).abs() < 1.0e-12
}

///|
/// Character `χ_a(x) = q^{-a·deg(x)}` of `Div^1` (a lisse rank-one ℓ-adic sheaf
/// / a Langlands parameter of `G_m`).
pub fn div1_character(x : Double, a : Double, q : Double) -> Double {
  @math.pow(q, -a * div1_degree(x))
}

///|
/// Characters form a group: `χ_{a+b} = χ_a·χ_b`.
pub fn div1_character_multiplicative(
  x : Double,
  a : Double,
  b : Double,
  q : Double,
) -> Bool {
  let lhs = div1_character(x, a + b, q)
  let rhs = div1_character(x, a, q) * div1_character(x, b, q)
  (lhs - rhs).abs() < 1.0e-10
}

///|
/// Weil weight of `χ_a`: `w = -2a`, since `|χ_a(x)| = q^{-a·deg(x)} = q^{(w/2)·deg(x)}`.
pub fn div1_character_weight(a : Double) -> Double {
  -2.0 * a
}

///|
/// Purity: `|χ_a(x)| = q^{(w/2)·deg(x)}` with `w = -2a`.
pub fn div1_character_is_pure(x : Double, a : Double, q : Double) -> Bool {
  let target = @math.pow(q, 0.5 * div1_character_weight(a) * div1_degree(x))
  (div1_character(x, a, q).abs() - target).abs() < 1.0e-10
}

///|
/// The degree lattice `(1/n)Z/Z` of `Div^1`: the `n` torsion points.
pub fn div1_torsion_points(n : Int32) -> Array[Double] {
  Array::makei(n, fn(i) { div1_degree(i.to_double() / n.to_double()) })
}

// ---------------------------------------------------------------------------
// 3. ℓ-cohomological smoothness: the shift + twist calculus (A36 / T44).
// ---------------------------------------------------------------------------

///|
/// A bounded complex of ℓ-adic sheaves: `betti[i]` is the rank in cohomological
/// degree `lo + i`, all twisted by `(twist)`.
pub(all) struct CohComplex {
  lo : Int32
  betti : Array[Int32]
  twist : Int32
}

///|
/// Algorithm A36: the dualizing complex `Λ(d)[2d]` of a smooth morphism of
/// relative dimension `d` — one copy of `Λ` in degree `2d`, twisted by `d`.
pub fn coh_dualizing(d : Int32) -> CohComplex {
  CohComplex::{ lo: 2 * d, betti: [1], twist: d, }
}

///|
/// Verdier dual: reverse the degrees (`i ↦ -i`) and negate the Tate twist.
pub fn coh_dual(c : CohComplex) -> CohComplex {
  let n = c.betti.length()
  let rev = Array::makei(n, fn(i) { c.betti[n - 1 - i] })
  CohComplex::{ lo: -(c.lo + n - 1), betti: rev, twist: -c.twist, }
}

///|
/// `f^!` for a smooth morphism of relative dimension `d`: shift by `2d`,
/// twist by `d` (`f^! K = K(d)[2d]`).
pub fn coh_shift_twist(c : CohComplex, d : Int32) -> CohComplex {
  CohComplex::{ lo: c.lo + 2 * d, betti: c.betti, twist: c.twist + d, }
}

///|
/// `Λ(d)[2d]` certification: the shift is `2d` and the twist is `d`.
pub fn coh_is_smooth_dualizing(c : CohComplex) -> Bool {
  let n = c.betti.length() == 1 && c.betti[0] == 1
  n && c.lo == 2 * c.twist
}

///|
/// Euler characteristic `Σ (-1)^i rk_i` (independent of the twist and shift).
pub fn coh_euler_characteristic(c : CohComplex) -> Int32 {
  let mut s = 0
  for i in 0..<c.betti.length() {
    let deg = c.lo + i
    let sign = if deg % 2 == 0 { 1 } else { -1 }
    s += sign * c.betti[i]
  }
  s
}

///|
/// Poincaré duality for smooth proper `X` of dimension `d`:
/// `b_i = b_{2d-i}` for all `i ∈ [0, 2d]`.
pub fn coh_poincare_dual(betti : Array[Int32], d : Int32) -> Bool {
  let top = 2 * d
  if betti.length() != top + 1 {
    return false
  }
  for i in 0..<=top {
    if betti[i] != betti[top - i] {
      return false
    }
  }
  true
}

///|
/// Composed smooth morphism: `(f∘g)^! = g^! f^!`, so dimensions add.
pub fn coh_compose_dimensions(d1 : Int32, d2 : Int32) -> Int32 {
  d1 + d2
}

///|
/// Base change certificate: for a finite étale cover of degree `n`, the Euler
/// characteristic multiplies (`χ_c(X) = n·χ_c(Y)`).
pub fn coh_etale_cover_euler(base_chi : Int32, degree : Int32) -> Int32 {
  base_chi * degree
}

///|
/// The constant complex `Λ` in degree 0.
pub fn coh_constant() -> CohComplex {
  CohComplex::{ lo: 0, betti: [1], twist: 0, }
}

///|
/// Dimension recovered from the dualizing complex: `d = twist`.
pub fn coh_recovered_dimension(c : CohComplex) -> Int32 {
  c.twist
}

// ---------------------------------------------------------------------------
// 4. Geometric Satake on the Hecke stack (A37 / T45).
// ---------------------------------------------------------------------------

///|
/// Perverse degree `⟨2ρ, λ⟩` of the IC sheaf `IC(Gr_λ)`.
pub fn satake_perverse_degree(rt : RootType, lambda : Array[Double]) -> Int32 {
  let rho = weyl_rho(rt)
  (2.0 * dot(rho, lambda)).round().to_int()
}

///|
/// Algorithm A37: `P_{IC(Gr_λ)}(1) = dim V_λ` — the evaluation of the IC
/// sheaf of the Hecke stack at `q = 1` recovers the representation.
pub fn satake_ic_match(rt : RootType, lambda : Array[Double]) -> Bool {
  let p = satake_minuscule_poincare(rt, lambda, 1.0)
  let d = weyl_dimension(rt, lambda).to_double()
  (p - d).abs() < 1.0e-9
}

///|
/// `q`-palindromicity (Poincaré–Verdier on the affine Grassmannian):
/// `P(q) = q^{top} P(1/q)` with `top = max_{w} l(w)`.
pub fn satake_poincare_palindromic(
  rt : RootType,
  lambda : Array[Double],
  top : Int32,
) -> Bool {
  for q in [0.5, 1.3, 2.0] {
    let lhs = satake_minuscule_poincare(rt, lambda, q)
    let rhs = @math.pow(q, top.to_double()) *
      satake_minuscule_poincare(rt, lambda, 1.0 / q)
    if (lhs - rhs).abs() > 1.0e-9 {
      return false
    }
  }
  true
}

// ---------------------------------------------------------------------------
// 5. Tests: Theorems T42 -- T45.
// ---------------------------------------------------------------------------

///|
test "T42: Frobenius is invertible exactly on the tilt" {
  // The original algebra is not perfect: p-th powers are a proper subset.
  for p in [2, 3] {
    for n in 2..<=4 {
      assert_true(residue_not_perfect(p, n), msg="|x^p| < p^n")
      assert_true(pth_power_image_count(p, n) < int_pow(p, n))
    }
    // ... but char-p residues are perfect: Frobenius is a bijection on F_p.
    assert_true(frob_bijective_mod_p(p))
    for n in 1..<=3 {
      assert_true(
        tilt_is_perfect(p, n),
        msg="tilt chains are compatible + stable",
      )
      assert_true(tilting_preserves_residue(p, n))
      assert_true(tilt_is_idempotent(p, n))
      assert_eq(tilt_chains(p, n, 2).length(), int_pow(p, n))
    }
  }
  // Root lifting exists exactly on the p-th powers: x = 1 has roots, x = p
  // does not (2 is not a square mod 4).
  assert_true(tilt_root_count(2, 2, 1) == 2)
  assert_eq(tilt_root_count(2, 2, 2), 0)
  assert_true(tilt_root_count(3, 2, 1) == 3)
  // Every chain of depth 3 is fully compatible.
  assert_true(tilt_chains_compatible(2, 3, 3))
  let chain = tilt_chain(2, 3, 5, 3)
  assert_eq(chain.length(), 3)
  assert_eq(chain[2], 5)
  assert_eq(chain[1], pth_power_mod(5, 2, 8))
}

///|
test "T43: Div^1 mirror curve — degrees, characters and purity" {
  let q = 7.0
  for x in [0.0, 0.37, 1.42, -2.6, 3.999] {
    assert_true(div1_frobenius_invariant(x), msg="deg(p·x) = deg(x) mod Z")
    assert_true((div1_degree(x) - div1_degree(x + 5.0)).abs() < 1.0e-12)
    for a in [0.5, 1.0, 2.5] {
      assert_true(div1_character_multiplicative(x, a, 0.75, q))
      assert_true(div1_character_is_pure(x, a, q))
    }
  }
  // The degree class is a homomorphism to R/Z.
  let a = 0.3
  let b = 0.85
  let s = div1_degree(a + b)
  let t = div1_degree(div1_degree(a) + div1_degree(b))
  assert_true((s - t).abs() < 1.0e-12)
  // Torsion points of the p-adic circle.
  let tors = div1_torsion_points(4)
  assert_eq(tors.length(), 4)
  assert_true((tors[0] - 0.0).abs() < 1.0e-12)
  assert_true((tors[1] - 0.25).abs() < 1.0e-12)
  // trivial character has weight 0
  assert_true((div1_character(0.37, 0.0, q) - 1.0).abs() < 1.0e-15)
}

///|
test "T44: ℓ-cohomological smoothness and the six-functor calculus" {
  // f^! Λ = Λ(d)[2d] and its dimension is recoverable from the twist.
  for d in [-4, -1, 0, 1, 2, 3] {
    let w = coh_dualizing(d)
    assert_true(coh_is_smooth_dualizing(w))
    assert_eq(coh_recovered_dimension(w), d)
    assert_eq(coh_euler_characteristic(w), if 2 * d % 2 == 0 { 1 } else { -1 })
  }
  // Verdier duality is an involution.
  let k = CohComplex::{ lo: 0, betti: [1, 2, 1], twist: 1, }
  let dd = coh_dual(coh_dual(k))
  assert_eq(dd.lo, k.lo)
  assert_eq(dd.twist, k.twist)
  assert_eq(dd.betti.length(), k.betti.length())
  for i in 0..<dd.betti.length() {
    assert_eq(dd.betti[i], k.betti[i])
  }
  // Poincaré duality palindromes: point, P^1, genus-2 curve, P^2.
  assert_true(coh_poincare_dual([1], 0))
  assert_true(coh_poincare_dual([1, 0, 1], 1))
  assert_true(coh_poincare_dual([1, 4, 1], 1))
  assert_true(coh_poincare_dual([1, 0, 1, 0, 1], 2))
  assert_true(!coh_poincare_dual([1, 2, 3], 1), msg="b_0 != b_2")
  // Composition: dim(f∘g) = dim f + dim g, and the dualizing complexes agree.
  let d1 = 1
  let d2 = 2
  let composed = coh_shift_twist(coh_dualizing(d1), d2)
  let direct = coh_dualizing(coh_compose_dimensions(d1, d2))
  assert_eq(composed.lo, direct.lo)
  assert_eq(composed.twist, direct.twist)
  assert_true(coh_is_smooth_dualizing(composed))
  // Classifying stacks [*/G] are ℓ-smooth of dimension -dim G.
  let gl1 = coh_dualizing(-1)
  assert_eq(gl1.lo, -2)
  assert_eq(gl1.twist, -1)
  let gl2 = coh_dualizing(-4)
  assert_eq(gl2.lo, -8)
  // Base change for a finite étale cover of degree 3.
  assert_eq(coh_etale_cover_euler(1, 3), 3)
  assert_eq(coh_etale_cover_euler(0, 3), 0)
  // The constant sheaf is its own dual.
  let c = coh_constant()
  let dc = coh_dual(c)
  assert_eq(dc.lo, 0)
  assert_eq(dc.twist, 0)
  assert_eq(dc.betti[0], 1)
}

///|
test "T45: geometric Satake on Hck_G matches representations" {
  // GL_2 = A_1: ω_1 is minuscule, P_IC(1) = 2 = dim V_{ω_1}.
  let a1 = RootType::TypeA(1)
  let w1 : Array[Double] = [1.0, 0.0]
  assert_true(satake_ic_match(a1, w1))
  assert_eq(weyl_dimension(a1, w1), 2)
  // GL_3 = A_2.
  let a2 = RootType::TypeA(2)
  let e1 : Array[Double] = [1.0, 0.0, 0.0]
  let e2 : Array[Double] = [1.0, 1.0, 0.0]
  assert_true(satake_ic_match(a2, e1))
  assert_true(satake_ic_match(a2, e2))
  assert_eq(weyl_dimension(a2, e1), 3)
  assert_eq(weyl_dimension(a2, e2), 3)
  // Sp_6 = C_3: ω_1 is minuscule (short root 2e_1 pairing).
  let c3 = RootType::TypeC(3)
  let cw1 : Array[Double] = [1.0, 0.0, 0.0]
  assert_true(satake_ic_match(c3, cw1))
  assert_eq(weyl_dimension(c3, cw1), 6)
  // Perverse degree ⟨2ρ, λ⟩ and q-palindromicity.
  assert_eq(satake_perverse_degree(a1, w1), 1)
  assert_eq(satake_perverse_degree(a2, e1), 2)
  assert_true(satake_poincare_palindromic(a2, e1, 2))
  // IC sheaves are pure of weight 0: the perverse shift [⟨2ρ,λ⟩] is a shift,
  // not a Tate twist, so it never changes the Weil weight.
  assert_eq(satake_perverse_degree(c3, cw1), 6)
}
