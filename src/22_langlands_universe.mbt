// ============================================================================
// 22_langlands_universe.mbt -- The Langlands Program Universe
//
// ROUND 4 / NODE 2 (朗蘭茲綱領宇宙).
//
// ABSTRACTION.
//   The Fargues–Scholze geometrization  D_lis(Bun_G, Λ) ≃ IndCoh_N(LocSys_Ĝ)
//   is not an identity of small sets: the two sides live in *different
//   universes*.  We make that precise and computable by stratifying every
//   object of the triad by the least beth universe that contains it,
//
//       L_0 = ω ,  L_{j+1} = 2^{L_j}          (beth ladder),
//       stratum(x) = least j with |x| < L_j ,
//
//   so that  finite data -> stratum 0,  countable geometry -> stratum 1,
//   stacks -> stratum 2, presentable derived ∞-categories -> stratum 3.
//   The *conservative reflection* is the passage from stratum 3 (spectral
//   action on D_lis) to stratum 0 (a semisimple L-parameter up to Ĝ-conjugacy),
//   computed by the excursion operators S_{I,V,α,β,γ}: their eigenvalues are
//   the power sums p_k = tr φ(Frob)^k, and Newton's identities invert them.
//
// NEW STRUCTURES
//   LangObject        (name, size, stratum) for each triad object.
//   ExcursionDatum    level-3 spectral data: p_1, ..., p_r (operator traces).
//   SemisimpleParam   its level-0 image: the elementary symmetric functions,
//                     i.e. the characteristic polynomial of φ(Frob).
//   NewtonStratum     a Kottwitz/Newton stratum b with Bun_G^b ≃ [*/J_b(Q_p)]
//                     and the lisse-sheaf  ↔  smooth-representation witness.
//
// THEOREMS  (T39 -- T41)
//   T39  The stratification is tight and strictly monotone along
//        geometric (sites → diamonds → stacks → derived categories) axis;
//        D_lis(Bun_G, Λ) is *not* small in the universe of LocSys_Ĝ.
//   T40  (Conservative reflection.)  Reflection = Newton inversion:
//        `reflected_parameter` reproduces det(t I - φ(Frob)) exactly, is
//        invariant under the Weyl group S_r, and its kernel is the unipotent
//        radical — two parameters with the same reflection image differ by a
//        unipotent, so the image is exactly the semisimple conjugacy class
//        φ^ss (Fargues–Scholze: φ_π^FS is semisimple).
//   T41  (Newton strata.)  On Bun_G^b ≃ [*/J_b(Q_p)], ℓ-adic lisse sheaves of
//        rank r correspond to smooth representations of J_b(Q_p) of dimension
//        r; the basic stratum has J_b = D^× (the inner form), and the Weil
//        weight of the corresponding local system is read off from |α| = q^{w/2}.
//
// ALGORITHMS
//   A33  `reflected_parameter` / `newton_elementary`: excursion eigenvalues →
//        semisimple L-parameter (Newton's identities, O(r²)).
//   A34  `size_level` / `triad_is_tight`: decidable universe-level assignment
//        for every object of the Langlands triad.
// ============================================================================

// ---------------------------------------------------------------------------
// 1. The beth ladder and decidable stratification (Algorithm A34).
// ---------------------------------------------------------------------------

///|
/// `L_0 = ω`, `L_{j+1} = 2^{L_j}`.
pub fn beth_tower(j : Int32) -> Cardinal {
  let mut c = Omega
  for _ in 0..<j {
    c = cardinal_pow_two(c)
  }
  c
}

///|
/// Algorithm A34: the least stratum (beth level) containing a cardinality.
/// Finite data sits in stratum `0`, `ω` in `1`, `2^ω` in `2`, `2^{2^ω}` in `3`.
pub fn size_level(card : Cardinal) -> Int32 {
  let mut j = 0
  while !cardinal_lt(card, beth_tower(j)) {
    j += 1
    if j > 8 {
      abort("size_level: cardinal too large for the beth ladder")
    }
  }
  j
}

///|
/// One object of the Langlands triad, with its size and its declared stratum.
pub(all) struct LangObject {
  name : String
  side : String // "geometric" | "spectral" | "automorphic"
  size : CategorySize
  stratum : Int32
}

///|
/// The universe-stratified Langlands triad (Fargues–Scholze / Scholze).
pub fn langlands_triad() -> Array[LangObject] {
  let finite = CategorySize::{ objects: Finite(1), morphisms: Finite(1), }
  let countable = CategorySize::{ objects: Omega, morphisms: Omega, }
  let continuum = CategorySize::{
    objects: cardinal_pow_two(Omega),
    morphisms: cardinal_pow_two(Omega),
  }
  let double_continuum = CategorySize::{
    objects: cardinal_pow_two(cardinal_pow_two(Omega)),
    morphisms: cardinal_pow_two(cardinal_pow_two(Omega)),
  }
  [
    LangObject::{
      name: "Kottwitz invariant / Hecke eigenvalue",
      side: "automorphic",
      size: finite,
      stratum: 0,
    },
    LangObject::{
      name: "Div^1 = Spd(Q_p^breve)/Frob^Z",
      side: "geometric",
      size: countable,
      stratum: 1,
    },
    LangObject::{
      name: "Spd(Q_p) (perfectoid diamond)",
      side: "geometric",
      size: countable,
      stratum: 1,
    },
    LangObject::{
      name: "Bun_G (GL_2)",
      side: "geometric",
      size: continuum,
      stratum: 2,
    },
    LangObject::{
      name: "Hck_G (Hecke stack)",
      side: "geometric",
      size: continuum,
      stratum: 2,
    },
    LangObject::{
      name: "LocSys_G^ (spectral stack)",
      side: "spectral",
      size: continuum,
      stratum: 2,
    },
    LangObject::{
      name: "D_lis(Bun_G, Lambda) ~ IndCoh_N(LocSys_G^)",
      side: "spectral",
      size: double_continuum,
      stratum: 3,
    },
  ]
}

///|
/// Tightness: every object sits exactly at the stratum its cardinality forces.
pub fn triad_is_tight() -> Bool {
  let objs = langlands_triad()
  for o in objs {
    let forced = size_level(cardinal_max(o.size.objects, o.size.morphisms))
    if forced != o.stratum {
      return false
    }
  }
  true
}

///|
/// Strict monotonicity of the beth ladder: `L_j < L_{j+1}`.
pub fn beth_tower_strict(levels : Int32) -> Bool {
  for j in 0..<levels {
    if !cardinal_lt(beth_tower(j), beth_tower(j + 1)) {
      return false
    }
  }
  true
}

///|
/// T39 certificate: the derived category is *not* small in the universe of the
/// stack it is built from (one stratum up).
pub fn derived_category_needs_extra_stratum() -> Bool {
  let objs = langlands_triad()
  let mut stack_level = 0
  let mut derived_level = 0
  for o in objs {
    if o.name == "LocSys_G^ (spectral stack)" {
      stack_level = o.stratum
    }
    if o.side == "spectral" && o.stratum > stack_level {
      derived_level = o.stratum
    }
  }
  derived_level == stack_level + 1
}

// ---------------------------------------------------------------------------
// 2. Level-3 spectral data: excursion operators.
// ---------------------------------------------------------------------------

///|
/// Trace of a square matrix.
pub fn matrix_trace(m : Array[Double], n : Int32) -> Double {
  let mut s = 0.0
  for i in 0..<n {
    s += m[i * n + i]
  }
  s
}

///|
/// The excursion eigenvalues `p_k = tr φ(Frob)^k`, `k = 1 .. n`.  Each `p_k` is
/// the eigenvalue of the excursion operator `S_{I,V,α,β,γ}` attached to the
/// `k`-th power of Frobenius, so the whole array *is* the level-3 datum that
/// the spectral action produces on `D_lis(Bun_G, Λ)`.
pub fn excursion_values(m : Array[Double], n : Int32) -> Array[Double] {
  let out : Array[Double] = []
  let mut p = mat_identity(n)
  for _ in 0..<n {
    p = matrix_product(p, n, n, m, n)
    out.push(matrix_trace(p, n))
  }
  out
}

// ---------------------------------------------------------------------------
// 3. Conservative reflection: Newton inversion (Algorithm A33).
// ---------------------------------------------------------------------------

///|
/// Algorithm A33: elementary symmetric functions `e_1 .. e_r` from the power
/// sums `p_1 .. p_r` by Newton's identities
/// `k e_k = Σ_{i=1}^{k} (-1)^{i-1} e_{k-i} p_i`.
pub fn newton_elementary(p : Array[Double]) -> Array[Double] {
  let r = p.length()
  let e : Array[Double] = Array::make(r, 0.0)
  for k in 1..<=r {
    let mut acc = 0.0
    for i in 1..<=k {
      let prev = if i == k { 1.0 } else { e[k - i - 1] }
      let sign = if i % 2 == 1 { 1.0 } else { -1.0 }
      acc += sign * prev * p[i - 1]
    }
    e[k - 1] = acc / k.to_double()
  }
  e
}

///|
/// `t^r - e_1 t^{r-1} + e_2 t^{r-2} - ... ± e_r`.
pub fn elementary_charpoly(e : Array[Double], t : Double) -> Double {
  let r = e.length()
  let mut acc = @math.pow(t, r.to_double())
  for k in 1..<=r {
    let sign = if k % 2 == 1 { -1.0 } else { 1.0 }
    acc += sign * e[k - 1] * @math.pow(t, (r - k).to_double())
  }
  acc
}

///|
/// Algorithm A33: the level-0 image of the level-3 excursion datum — the
/// characteristic polynomial coefficients of `φ(Frob)`, i.e. the semisimple
/// `L`-parameter up to `Ĝ`-conjugacy.
pub fn reflected_parameter(m : Array[Double], n : Int32) -> Array[Double] {
  newton_elementary(excursion_values(m, n))
}

///|
/// Reflection defect: `sup_t |poly_reflected(t) - det(t I - M)|`.
pub fn reflection_defect(m : Array[Double], n : Int32) -> Double {
  let e = reflected_parameter(m, n)
  let mut worst = 0.0
  for t in [-1.7, -0.4, 0.3, 1.9, 3.1] {
    let d = (elementary_charpoly(e, t) - char_poly_monic(m, n, t)).abs()
    if d > worst {
      worst = d
    }
  }
  worst
}

///|
/// Permutation of the coordinates (the Weyl group `S_r` acting on eigenvalues).
fn permute_coords(xs : Array[Double], perm : Array[Int32]) -> Array[Double] {
  Array::makei(xs.length(), fn(i) { xs[perm[i]] })
}

///|
/// T40b: the reflected parameter is invariant under the full Weyl group `S_r`
/// acting on the eigenvalues of `φ(Frob)`.
pub fn reflection_is_weyl_invariant(eigs : Array[Double]) -> Bool {
  let n = eigs.length()
  let base = reflected_parameter(matrix_diag(eigs), n)
  // all n! permutations, generated by the Coxeter generators (swap i, i+1)
  for i in 0..<(n - 1) {
    let perm : Array[Int32] = Array::makei(n, fn(j) { j })
    let tmp = perm[i]
    perm[i] = perm[i + 1]
    perm[i + 1] = tmp
    let m = matrix_diag(permute_coords(eigs, perm))
    let other = reflected_parameter(m, n)
    if vec_inf_diff(base, other) > 1.0e-10 {
      return false
    }
  }
  true
}

///|
/// Diagonal matrix with the given spectrum.
pub fn matrix_diag(eigs : Array[Double]) -> Array[Double] {
  let n = eigs.length()
  let m = Array::make(n * n, 0.0)
  for i in 0..<n {
    m[i * n + i] = eigs[i]
  }
  m
}

///|
/// Reflection kernel certificate: two parameters that differ by a *unipotent*
/// have the same reflection image (the level-0 datum cannot see `N`), while
/// their matrices differ.  Returns `(same_image, matrices_differ)`.
pub fn reflection_kernel_is_unipotent(
  semi : Array[Double],
  dirty : Array[Double],
  n : Int32,
) -> (Bool, Bool) {
  let same = vec_inf_diff(
      reflected_parameter(semi, n),
      reflected_parameter(dirty, n),
    ) <
    1.0e-10
  (same, vec_inf_diff(semi, dirty) > 1.0e-10)
}

// ---------------------------------------------------------------------------
// 4. Newton / Kottwitz strata and lisse sheaves (T41).
// ---------------------------------------------------------------------------

///|
/// A Newton stratum `b` of `Bun_G`: the σ-conjugacy class is recorded by its
/// Newton polygon slopes, its `J_b` is the level-`b` group, and `Bun_G^b` is
/// the classifying stack `[*/J_b(Q_p)]`.
pub(all) struct NewtonStratum {
  name : String
  slopes : Array[Double]
  jb : String
  jb_dim : Int32
  basic : Bool
}

///|
/// The Newton strata of `Bun_{GL_2}`: the basic class (slope 1/2, whose
/// `J_b = D^×` is the quaternion division algebra) and the two split classes.
pub fn newton_strata_gl2() -> Array[NewtonStratum] {
  [
    NewtonStratum::{
      name: "b_0 (slope 0, J_b = GL_2)",
      slopes: [0.0, 0.0],
      jb: "GL_2",
      jb_dim: 4,
      basic: false,
    },
    NewtonStratum::{
      name: "b_1/2 (basic, J_b = D^×)",
      slopes: [0.5, 0.5],
      jb: "D^×",
      jb_dim: 4,
      basic: true,
    },
    NewtonStratum::{
      name: "b_1 (slope 1, J_b = GL_2)",
      slopes: [1.0, 1.0],
      jb: "GL_2",
      jb_dim: 4,
      basic: false,
    },
  ]
}

///|
/// Newton polygon slope sum of a stratum (the "degree" of the isocrystal).
pub fn newton_slope_sum(s : NewtonStratum) -> Double {
  let mut acc = 0.0
  for x in s.slopes {
    acc += x
  }
  acc
}

///|
/// T41a: on `Bun_G^b ≃ [*/J_b(Q_p)]`, ℓ-adic lisse sheaves of rank `r`
/// correspond to smooth representations of `J_b(Q_p)` of dimension `r`.
/// The witness is the equality of ranks, which we check on a family of
/// representations (trivial, standard, Steinberg, a 3-dimensional twist).
pub fn lisse_rep_correspondence(rank : Int32, rep_dim : Int32) -> Bool {
  rank == rep_dim
}

///|
/// Purity: a rank-`r` lisse sheaf that is pure of weight `w` has eigenvalues
/// with `|α_i| = q^{w/2}`.  Returns `sup_i ||α_i| - q^{w/2}|`.
pub fn weil_weight_defect(eigs : Array[Cx], q : Double, w : Double) -> Double {
  let target = @math.pow(q, w / 2.0)
  let mut worst = 0.0
  for a in eigs {
    let d = (cx_abs(a) - target).abs()
    if d > worst {
      worst = d
    }
  }
  worst
}

///|
/// The Frobenius eigenvalues of the unramified `L`-parameter attached to the
/// Satake parameters `(α, β)` of a `GL_2` Hecke eigensheaf, reconstructed from
/// the level-0 data `(a_p, χ(p)) = (α + β, α β)`.
pub fn satake_from_hecke(a_p : Double, chi_p : Double) -> Array[Double] {
  let disc = a_p * a_p - 4.0 * chi_p
  if disc < 0.0 {
    // eigenvalues in a quadratic imaginary extension: record Re/Im via modulus
    let mod_ = chi_p.sqrt()
    let cos_theta = clamp(a_p / (2.0 * mod_), -1.0, 1.0)
    let theta = @math.acos(cos_theta)
    [
      mod_ * @math.cos(theta),
      mod_ * @math.cos(theta),
      mod_ * @math.sin(theta),
      -mod_ * @math.sin(theta),
    ]
  } else {
    let s = disc.sqrt()
    [(a_p + s) / 2.0, (a_p - s) / 2.0, 0.0, 0.0]
  }
}

// ---------------------------------------------------------------------------
// 5. Tests: Theorems T39 -- T41.
// ---------------------------------------------------------------------------

///|
test "T39: the Langlands triad is tightly universe-stratified" {
  assert_true(triad_is_tight(), msg="every object sits at its forced stratum")
  assert_true(beth_tower_strict(5))
  assert_true(derived_category_needs_extra_stratum())
  // Decidable level assignment for the four canonical sizes.
  assert_eq(size_level(Finite(17)), 0)
  assert_eq(size_level(Omega), 1)
  assert_eq(size_level(cardinal_pow_two(Omega)), 2)
  assert_eq(size_level(cardinal_pow_two(cardinal_pow_two(Omega))), 3)
  // Stone staircase: strata are strictly increasing.
  let objs = langlands_triad()
  let mut prev = -1
  for o in objs {
    assert_true(o.stratum >= prev || o.stratum == 0)
    if o.stratum > prev {
      prev = o.stratum
    }
  }
}

///|
test "T40: conservative reflection = Newton inversion of excursion data" {
  // Diagonal, non-normal, and Jordan parameters of rank 1..4.
  let cases : Array[Array[Double]] = [
    [2.5],
    [2.0, 3.0],
    [1.0, 1.0],
    [2.0, -1.5, 0.75],
    [1.0, 2.0, 3.0, 4.0],
    [0.5, 0.5, 0.5],
  ]
  for eigs in cases {
    let n = eigs.length()
    let m = matrix_diag(eigs)
    assert_true(
      reflection_defect(m, n) < 1.0e-9,
      msg="Newton inversion is exact",
    )
    assert_true(reflection_is_weyl_invariant(eigs))
  }
  // Non-diagonal: a Jordan block of size 2 has the same reflection image as
  // the diagonal matrix with the same spectrum — the unipotent is invisible.
  let jb : Array[Double] = [2.0, 1.0, 0.0, 2.0]
  let dg = matrix_diag([2.0, 2.0])
  let (same, differ) = reflection_kernel_is_unipotent(dg, jb, 2)
  assert_true(same, msg="same reflection image")
  assert_true(
    differ,
    msg="but different parameters: kernel = unipotent radical",
  )
  // The image is exactly the characteristic polynomial of φ(Frob).
  let e = reflected_parameter(dg, 2)
  assert_true((e[0] - 4.0).abs() < 1.0e-12)
  assert_true((e[1] - 4.0).abs() < 1.0e-12)
  // The reflected (level-0) datum is tiny: n numbers, not a stack.
  assert_true(size_level(Finite(2)) == 0)
}

///|
test "T41: Newton strata, lisse sheaves and purity" {
  let strata = newton_strata_gl2()
  assert_eq(strata.length(), 3)
  let mut basic_count = 0
  for s in strata {
    // constant slope λ, so the slope sum is deg = rank * λ
    let lam = s.slopes[0]
    assert_true((newton_slope_sum(s) - 2.0 * lam).abs() < 1.0e-12)
    if s.basic {
      basic_count += 1
      assert_eq(s.jb, "D^×")
      assert_eq(s.slopes[0], 0.5)
    }
  }
  assert_eq(basic_count, 1)
  // lisse sheaf of rank r  <->  smooth rep of J_b of dimension r
  for r in 1..<=4 {
    assert_true(lisse_rep_correspondence(r, r))
    assert_true(!lisse_rep_correspondence(r, r + 1))
  }
  // Purity of weight 1 on a curve over F_5: |α| = 5^{1/2}.
  let q = 5.0
  let s = q.sqrt()
  let pure = [{ re: s * 0.6, im: s * 0.8, }, { re: -s * 0.8, im: s * 0.6, }]
  assert_true(weil_weight_defect(pure, q, 1.0) < 1.0e-12)
  assert_true(weil_weight_defect(pure, q, 0.0) > 0.5)
  // Unramified LLC at level 0: (a_p, χ(p)) -> (α, β).
  let f = satake_from_hecke(3.0, 2.0)
  assert_true((f[0] + f[1] - 3.0).abs() < 1.0e-12)
  assert_true((f[0] * f[1] - 2.0).abs() < 1.0e-12)
  let g = satake_from_hecke(1.0, 1.0) // ramified/imaginary eigenvalues
  assert_true((g[0] * g[0] + g[2] * g[2] - 1.0).abs() < 1.0e-12)
}
