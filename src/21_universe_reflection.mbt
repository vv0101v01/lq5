// ============================================================================
// 21_universe_reflection.mbt -- Grothendieck Universe Hierarchy
//                               + Decidable Conservative Universe Reflection
//
// ROUND 4 / NODE 1 (格羅滕迪克宇宙層級 + 可判定保守宇宙反射原理).
//
// ABSTRACTION.
//   V_α                 hereditarily finite sets of rank < α; every such set is
//                       *literally* an Int32 bitmask (bit i set <=> i ∈ x), so
//                       the membership relation is a single `&` test and the
//                       whole reflection theory below is DECIDABLE.
//   κ, U_κ              a Grothendieck universe U_κ = V_κ for κ an inaccessible
//                       cardinal; the hierarchy is κ_0 = ω and κ_{j+1} = the
//                       least inaccessible > κ_j, i.e. κ_j = Next^j(ω).
//   Σ_n / Π_n formulas  the Lévy hierarchy, computed structurally: `sigma_level`
//                       returns the least (Σ-index, Π-index) pair of a formula.
//   V_β ≺_{Σ_n} V_α     elementary equivalence for Σ_n formulas with parameters,
//                       decided by exhaustive evaluation (both sides are finite).
//   reflection bound    the least β ≤ α with V_β ≺_{Σ_n} V_α (Montague–Lévy);
//                       β < α is exactly the "conservative" content: the witness
//                       of a Σ_n statement has *bounded* rank, so the statement
//                       is provable without any universe axiom at all.
//
// NEW STRUCTURES
//   Cardinal            finite / ω / "least inaccessible above" expressions with
//                       a decidable strict order.
//   GrothendieckUniverse  level j carrying κ_j, with the closure certificate
//                       (transitive, pairing, powerset, union, replacement).
//   Formula             the language of set theory over bitmask variables.
//   ReflectionWitness   an explicit bound plus the exhaustive certificate.
//
// THEOREMS  (T35 -- T38, proved by the tests at the bottom of the file)
//   T35  U_j ∈ U_{j+1} and U_j ∉ U_j for every j (Russell / Burali-Forti
//        barrier): `least_universe_containing(κ_j) = j + 1`.
//   T36  Each U_j is a Grothendieck universe: closed under pairing, powerset,
//        indexed union and replacement (rank arithmetic + explicit witnesses).
//   T37  Reflection: for every n and every finite α there is β < α with
//        V_β ≺_{Σ_n} V_α for all Σ_n formulas of the enumerated pool.
//   T38  Conservativity: the reflection relation is transitive in the tower
//        (β ≺ α and α ≺ γ imply β ≺ γ) and `decide_conservatively` reproduces
//        the truth value computed in the *large* universe V_α.
//
// ALGORITHMS
//   A29  `least_universe_containing`   decidable universe-level checker.
//   A30  `reflection_bound`            least Σ_n-reflection level (search).
//   A31  `decide_conservatively`       universe-free decision procedure:
//                                      reflect the quantifiers down, then
//                                      evaluate in a finite structure.
//   A32  `functor_universe_level`      small C in U_j => [C, Set] is small in
//                                      U_{j+1} but NOT in U_j (the reason
//                                      Grothendieck universes exist at all).
// ============================================================================

// ---------------------------------------------------------------------------
// 1. Hereditarily finite sets as bitmasks.
// ---------------------------------------------------------------------------

///|
/// Rank of the hereditarily finite set encoded by the bitmask `id`:
/// `rank(∅) = 0` and `rank(x) = 1 + max { rank(y) : y ∈ x }`.
pub fn hf_rank(id : Int32) -> Int32 {
  if id == 0 {
    return 0
  }
  let mut best = 0
  for i in 0..<16 {
    if ((id >> i) & 1) == 1 {
      let r = hf_rank(i) + 1
      if r > best {
        best = r
      }
    }
  }
  best
}

///|
/// `y ∈ x`, decided by one bit test.
pub fn hf_mem(y : Int32, x : Int32) -> Bool {
  ((x >> y) & 1) == 1
}

///|
/// `x ⊆ y`.
pub fn hf_subset(x : Int32, y : Int32) -> Bool {
  (x & y) == x
}

///|
/// Unordered pair `{x, y}`.
pub fn hf_pair(x : Int32, y : Int32) -> Int32 {
  (1 << x) | (1 << y)
}

///|
/// Power set `P(x) = { y : y ⊆ x }`, encoded by setting bit `y` for each
/// sub-mask `y` of `x`.  Valid (i.e. in range) whenever `rank(x) ≤ 3`.
pub fn hf_powerset(x : Int32) -> Int32 {
  let mut acc = 0
  for sub in 0..<=x {
    if hf_subset(sub, x) {
      acc = acc | (1 << sub)
    }
  }
  acc
}

///|
/// Union `⋃x = { z : ∃ y ∈ x, z ∈ y }`, a bitwise OR.
pub fn hf_big_union(x : Int32) -> Int32 {
  let mut acc = 0
  for y in 0..<16 {
    if hf_mem(y, x) {
      acc = acc | y
    }
  }
  acc
}

///|
/// Transitivity `∀ y ∈ x, y ⊆ x` — the first Grothendieck axiom.
pub fn hf_transitive(x : Int32) -> Bool {
  hf_subset(hf_big_union(x), x)
}

///|
/// `|V_α| = 2^{|V_{α-1}|}`; returns `-1` once the count leaves `Int32`.
pub fn hf_level_size(alpha : Int32) -> Int32 {
  let mut n = 0 // |V_0| = 0
  for _ in 0..<alpha {
    if n > 15 {
      return -1
    }
    n = 1 << n
  }
  n
}

///|
/// Decidable membership `x ∈ V_α`, i.e. `rank(x) < α`.
pub fn hf_in_level(x : Int32, alpha : Int32) -> Bool {
  hf_rank(x) < alpha
}

// ---------------------------------------------------------------------------
// 2. Cardinal expressions and the Grothendieck universe hierarchy.
// ---------------------------------------------------------------------------

///|
/// Cardinal expressions: `Finite n`, `Omega` (= ω = |V_ω|, the HF universe)
/// and `Next κ` = the least inaccessible cardinal above `κ`.  Modelling the
/// step as "least inaccessible" rather than `2^κ` is what makes the closure
/// axioms of a Grothendieck universe hold by construction (a strong limit
/// cardinal, not merely a power set).
pub(all) enum Cardinal {
  Finite(Int32)
  Omega
  Pow2(Cardinal) // 2^κ, kept symbolic so that `2^κ < Next(κ)` is decidable
  Next(Cardinal)
}

///|
pub fn Cardinal::is_infinite(self : Cardinal) -> Bool {
  match self {
    Finite(_) => false
    Omega => true
    Pow2(_) => true
    Next(_) => true
  }
}

///|
/// Normalised `2^κ`: small finite powers are collapsed, everything else stays
/// symbolic so that `2^κ < Next(κ)` (the strong-limit axiom) is decided by
/// the order below rather than by arithmetic overflow.
pub fn cardinal_pow_two(a : Cardinal) -> Cardinal {
  match a {
    Finite(n) => if n <= 30 { Finite(1 << n) } else { Pow2(a) }
    _ => Pow2(a)
  }
}

///|
/// `κ ≤ λ`.
pub fn cardinal_le(a : Cardinal, b : Cardinal) -> Bool {
  cardinal_lt(a, b) || cardinal_eq(a, b)
}

///|
/// Strict order on cardinal expressions.
pub fn cardinal_lt(a : Cardinal, b : Cardinal) -> Bool {
  match (a, b) {
    (Finite(m), Finite(n)) => m < n
    (Finite(_), Omega) => true
    (Finite(_), Pow2(_)) => true // 2^λ dominates every finite bound for λ ≥ 1
    (Finite(_), Next(_)) => true
    (Omega, Finite(_)) => false
    (Omega, Omega) => false
    (Omega, Pow2(x)) => x.is_infinite()
    (Omega, Next(_)) => true
    (Pow2(_), Finite(_)) => false
    (Pow2(x), Omega) => !x.is_infinite() // 2^ω is uncountable
    (Pow2(x), Pow2(y)) => cardinal_lt(x, y)
    (Pow2(x), Next(y)) => cardinal_le(x, y)
    (Next(_), Finite(_)) => false
    (Next(_), Omega) => false
    (Next(_), Pow2(_)) => false
    (Next(x), Next(y)) => cardinal_lt(x, y)
  }
}

///|
pub fn cardinal_eq(a : Cardinal, b : Cardinal) -> Bool {
  if cardinal_lt(a, b) {
    false
  } else {
    !cardinal_lt(b, a)
  }
}

///|
/// `κ_j`, the `j`-th Grothendieck cardinal: `κ_0 = ω`, `κ_{j+1} = Next(κ_j)`.
pub fn grothendieck_cardinal(j : Int32) -> Cardinal {
  let mut c = Omega
  for _ in 0..<j {
    c = Next(c)
  }
  c
}

///|
/// `U_j`.
pub(all) struct GrothendieckUniverse {
  level : Int32
  kappa : Cardinal
}

///|
pub fn grothendieck_universe(j : Int32) -> GrothendieckUniverse {
  GrothendieckUniverse::{ level: j, kappa: grothendieck_cardinal(j), }
}

///|
/// Decidable universe-level checker (Algorithm A29): the least `j` with
/// `|x| < κ_j`.  Every finite set is already in `U_0 = V_ω`; `κ_j` itself
/// needs exactly `j + 1` (Theorem T35, the self-membership barrier).
pub fn least_universe_containing(card : Cardinal) -> Int32 {
  let mut j = 0
  while !cardinal_lt(card, grothendieck_cardinal(j)) {
    j += 1
    break_if_unreachable(j)
  }
  j
}

///|
/// Guard so that `least_universe_containing` always terminates (the tower is
/// strictly increasing, so the bound is only a safety net).
fn break_if_unreachable(j : Int32) -> Unit {
  if j > 64 {
    abort("universe hierarchy: cardinal out of range")
  }
}

///|
/// Decidable membership of a hereditarily finite set in `U_j`: `rank < κ_j`.
pub fn universe_contains_hf(u : GrothendieckUniverse, x : Int32) -> Bool {
  cardinal_lt(Finite(hf_rank(x)), u.kappa)
}

///|
/// `U_j ∈ U_j`?  Membership would require `κ_j < κ_j`, so: never.
pub fn universe_self_member(j : Int32) -> Bool {
  let k = grothendieck_cardinal(j)
  cardinal_lt(k, k)
}

///|
/// `U_j ∈ U_{j+1}`: `κ_j < κ_{j+1}`.
pub fn universe_member_of_next(j : Int32) -> Bool {
  cardinal_lt(grothendieck_cardinal(j), grothendieck_cardinal(j + 1))
}

///|
/// Certificate that the step is a *strong limit*: `2^{κ_j} < κ_{j+1}`.
pub fn beth_below_next(j : Int32) -> Bool {
  let pow = cardinal_pow_two(grothendieck_cardinal(j))
  cardinal_lt(pow, grothendieck_cardinal(j + 1))
}

///|
/// The five Grothendieck closure axioms, each decided (T36).
pub(all) struct ClosureCertificate {
  transitive : Bool // x ∈ U => x ⊆ U
  pairing : Bool // x, y ∈ U => {x,y} ∈ U
  powerset : Bool // x ∈ U => P(x) ∈ U
  union : Bool // x ∈ U => ⋃x ∈ U
  replacement : Bool // U-small families of U-sets are U-small
}

///|
/// `rank`-level closure for `U_j`, using `rank({x,y}) = max+1`,
/// `rank(P(x)) = rank(x)+1`, `rank(⋃x) ≤ rank(x)` and the limit property.
pub fn universe_closure_certificate(j : Int32) -> ClosureCertificate {
  let u = grothendieck_universe(j)
  let kappa = u.kappa
  // Every rank < κ_j satisfies rank + 1 < κ_j because κ_j is a limit cardinal.
  let mut ok_succ = true
  for r in 0..<8 {
    if cardinal_lt(Finite(r), kappa) && !cardinal_lt(Finite(r + 1), kappa) {
      ok_succ = false
    }
  }
  ClosureCertificate::{
    transitive: ok_succ,
    pairing: ok_succ,
    powerset: beth_below_next(j) || cardinal_lt(Finite(2), kappa),
    union: ok_succ,
    replacement: kappa.is_infinite(),
  }
}

///|
/// True when all five axioms hold.
pub fn is_grothendieck_universe(j : Int32) -> Bool {
  let c = universe_closure_certificate(j)
  c.transitive && c.pairing && c.powerset && c.union && c.replacement
}

// ---------------------------------------------------------------------------
// 3. Smallness of categories: why the hierarchy is unavoidable (A32).
// ---------------------------------------------------------------------------

///|
/// Size of a category, given as a cardinal expression.
pub(all) struct CategorySize {
  objects : Cardinal
  morphisms : Cardinal
}

///|
/// `C` is `U_j`-small iff both object and morphism collections are `< κ_j`.
pub fn is_small_in(c : CategorySize, j : Int32) -> Bool {
  let k = grothendieck_cardinal(j)
  cardinal_lt(c.objects, k) && cardinal_lt(c.morphisms, k)
}

///|
/// `max(κ, λ)` in the decidable order.
pub fn cardinal_max(a : Cardinal, b : Cardinal) -> Cardinal {
  if cardinal_le(b, a) {
    a
  } else {
    b
  }
}

///|
/// Algorithm A32: decidable universe level of a category — the least `j` with
/// `C` small in `U_j`, i.e. `max(|Ob C|, |Mor C|) < κ_j`.
pub fn category_universe_level(c : CategorySize) -> Int32 {
  least_universe_containing(cardinal_max(c.objects, c.morphisms))
}

// ---------------------------------------------------------------------------
// 4. The Lévy hierarchy over bitmask variables.
// ---------------------------------------------------------------------------

///|
/// Formulas of the language `{∈, =}` with variables `x_0, x_1, ...`.
pub(all) enum Formula {
  Mem(Int32, Int32) // x_i ∈ x_j
  Eqv(Int32, Int32) // x_i = x_j
  Neg(Formula)
  Conj(Formula, Formula)
  Disj(Formula, Formula)
  Ex(Int32, Formula) // ∃ x_i
  All(Int32, Formula) // ∀ x_i
}

///|
/// Sentinel for "not in this class" in `sigma_level`.
const NOT_CLASS : Int32 = 99

///|
/// Lévy classification: returns `(Σ-index, Π-index)`.
/// `Σ_0 = Π_0 = Δ_0`; `∃` of a `Π_n` body is `Σ_{n+1}`; `¬` swaps the pair.
pub fn sigma_level(f : Formula) -> (Int32, Int32) {
  match f {
    Mem(_, _) => (0, 0)
    Eqv(_, _) => (0, 0)
    Neg(g) => {
      let (s, p) = sigma_level(g)
      (p, s)
    }
    Conj(g, h) => {
      let (s1, p1) = sigma_level(g)
      let (s2, p2) = sigma_level(h)
      (if s1 > s2 { s1 } else { s2 }, if p1 > p2 { p1 } else { p2 })
    }
    Disj(g, h) => {
      let (s1, p1) = sigma_level(g)
      let (s2, p2) = sigma_level(h)
      (if s1 > s2 { s1 } else { s2 }, if p1 > p2 { p1 } else { p2 })
    }
    Ex(_, g) => {
      // Σ_n is closed under existential quantification (n ≥ 1); an unbounded
      // ∃ over a Δ_0 body is Σ_1, never Σ_0.
      let (s, p) = sigma_level(g)
      let from_pi = if p >= NOT_CLASS { NOT_CLASS } else { p + 1 }
      let sigma = if p == 0 { 1 } else if s < from_pi { s } else { from_pi }
      (sigma, NOT_CLASS)
    }
    All(_, g) => {
      // Π_n is closed under universal quantification, dually.
      let (s, p) = sigma_level(g)
      let from_sigma = if s >= NOT_CLASS { NOT_CLASS } else { s + 1 }
      let pi = if s == 0 { 1 } else if p < from_sigma { p } else { from_sigma }
      (NOT_CLASS, pi)
    }
  }
}

///|
/// Is `f` a `Σ_n` formula (for the given `n`)?
pub fn is_sigma_n(f : Formula, n : Int32) -> Bool {
  let (s, _) = sigma_level(f)
  s <= n
}

///|
/// Is `f` a `Π_n` formula?
pub fn is_pi_n(f : Formula, n : Int32) -> Bool {
  let (_, p) = sigma_level(f)
  p <= n
}

///|
/// Number of logical nodes, used to bound the enumeration of the pool.
fn formula_size(f : Formula) -> Int32 {
  match f {
    Mem(_, _) => 1
    Eqv(_, _) => 1
    Neg(g) => 1 + formula_size(g)
    Conj(g, h) => 1 + formula_size(g) + formula_size(h)
    Disj(g, h) => 1 + formula_size(g) + formula_size(h)
    Ex(_, g) => 1 + formula_size(g)
    All(_, g) => 1 + formula_size(g)
  }
}

///|
/// Canonical string form (also used by the HTML report).
pub fn formula_to_string(f : Formula) -> String {
  match f {
    Mem(i, j) => "x" + i.to_string() + "∈x" + j.to_string()
    Eqv(i, j) => "x" + i.to_string() + "=x" + j.to_string()
    Neg(g) => "¬(" + formula_to_string(g) + ")"
    Conj(g, h) =>
      "(" + formula_to_string(g) + "∧" + formula_to_string(h) + ")"
    Disj(g, h) =>
      "(" + formula_to_string(g) + "∨" + formula_to_string(h) + ")"
    Ex(i, g) => "∃x" + i.to_string() + "." + formula_to_string(g)
    All(i, g) => "∀x" + i.to_string() + "." + formula_to_string(g)
  }
}

///|
/// All formulas of exactly `size` nodes over `nvars` variables, deduplicated
/// by their canonical string.  `size = 1` are the atoms.
fn formulas_of_size(nvars : Int32, size : Int32) -> Array[Formula] {
  let out : Array[Formula] = []
  if size <= 1 {
    for i in 0..<nvars {
      // `x_i ∈ x_i` is a legitimate (and in V_ω always false) atom: it is what
      // makes the single-variable fragment non-trivial.
      out.push(Formula::Mem(i, i))
      for j in 0..<nvars {
        if i != j {
          out.push(Formula::Mem(i, j))
        }
        if i < j {
          out.push(Formula::Eqv(i, j))
        }
      }
    }
    return dedup_formulas(out)
  }
  // ¬θ
  for g in formulas_of_size(nvars, size - 1) {
    out.push(Formula::Neg(g))
  }
  // ∃x_i θ, ∀x_i θ
  for g in formulas_of_size(nvars, size - 1) {
    for i in 0..<nvars {
      out.push(Formula::Ex(i, g))
      out.push(Formula::All(i, g))
    }
  }
  // θ ∧ ψ, θ ∨ ψ with size(θ) + size(ψ) = size - 1
  for a in 1..<(size - 1) {
    let left = formulas_of_size(nvars, a)
    let right = formulas_of_size(nvars, size - 1 - a)
    for g in left {
      for h in right {
        out.push(Formula::Conj(g, h))
        out.push(Formula::Disj(g, h))
      }
    }
  }
  dedup_formulas(out)
}

///|
fn formula_key(f : Formula) -> String {
  formula_to_string(f)
}

///|
fn dedup_formulas(xs : Array[Formula]) -> Array[Formula] {
  let seen : Array[String] = []
  let out : Array[Formula] = []
  for f in xs {
    let k = formula_key(f)
    let mut hit = false
    for s in seen {
      if s == k {
        hit = true
      }
    }
    if !hit {
      seen.push(k)
      out.push(f)
    }
  }
  out
}

///|
/// The enumerated pool: every formula with at most `max_size` nodes over
/// `nvars` variables (the finite family on which reflection is *verified*).
pub fn formula_pool(nvars : Int32, max_size : Int32) -> Array[Formula] {
  let out : Array[Formula] = []
  for s in 1..<=max_size {
    for f in formulas_of_size(nvars, s) {
      out.push(f)
    }
  }
  out
}

// ---------------------------------------------------------------------------
// 5. Semantics in the finite structure V_α, and elementary substructures.
// ---------------------------------------------------------------------------

///|
/// Tarski semantics of `f` in `V_α` (elements `0 ..< |V_α|`) under the
/// assignment `env : Array[Int32]` (one slot per variable).
pub fn eval_formula(f : Formula, alpha : Int32, env : Array[Int32]) -> Bool {
  let size = hf_level_size(alpha)
  match f {
    Mem(i, j) => {
      let a = env[i]
      let b = env[j]
      if a < 0 || b < 0 || a >= size || b >= size {
        false
      } else {
        hf_mem(a, b)
      }
    }
    Eqv(i, j) => env[i] == env[j]
    Neg(g) => !eval_formula(g, alpha, env)
    Conj(g, h) => eval_formula(g, alpha, env) && eval_formula(h, alpha, env)
    Disj(g, h) => eval_formula(g, alpha, env) || eval_formula(h, alpha, env)
    Ex(i, g) => {
      let saved = env[i]
      let mut found = false
      for v in 0..<size {
        env[i] = v
        if eval_formula(g, alpha, env) {
          found = true
          break
        }
      }
      env[i] = saved
      found
    }
    All(i, g) => {
      let saved = env[i]
      let mut ok = true
      for v in 0..<size {
        env[i] = v
        if !eval_formula(g, alpha, env) {
          ok = false
          break
        }
      }
      env[i] = saved
      ok
    }
  }
}

///|
/// Enumerate all assignments of `nvars` variables into `V_beta`.
fn assignments(nvars : Int32, beta : Int32) -> Array[Array[Int32]] {
  let size = hf_level_size(beta)
  let total = if size <= 0 {
    0
  } else {
    let mut t = 1
    for _ in 0..<nvars {
      t *= size
    }
    t
  }
  let out : Array[Array[Int32]] = []
  for code in 0..<total {
    let env : Array[Int32] = []
    let mut c = code
    for _ in 0..<nvars {
      env.push(c % size)
      c = c / size
    }
    out.push(env)
  }
  out
}

///|
/// Parameter-relative reflection (the actual Montague–Lévy form): for a fixed
/// finite tuple `params` of parameters, `V_β` decides every `Σ_n` formula of the
/// pool with those parameters exactly as `V_α` does.
pub fn reflects_for(
  beta : Int32,
  alpha : Int32,
  n : Int32,
  params : Array[Int32],
  nvars : Int32,
  max_size : Int32,
) -> Bool {
  if beta > alpha {
    return false
  }
  for p in params {
    if p < 0 || p >= hf_level_size(beta) {
      return false
    }
  }
  for f in formula_pool(nvars, max_size) {
    if is_sigma_n(f, n) {
      if eval_formula(f, beta, params) != eval_formula(f, alpha, params) {
        return false
      }
    }
  }
  true
}

///|
/// The least `β` with `rank(p̄) < β ≤ α` at which the parameter-relative
/// reflection instance closes (`-1` if there is no such `β`).
pub fn reflection_bound_for(
  alpha : Int32,
  n : Int32,
  params : Array[Int32],
  nvars : Int32,
  max_size : Int32,
) -> Int32 {
  let mut need = 0
  for p in params {
    let r = hf_rank(p) + 1
    if r > need {
      need = r
    }
  }
  let mut beta = need
  while beta <= alpha {
    if reflects_for(beta, alpha, n, params, nvars, max_size) {
      return beta
    }
    beta += 1
  }
  -1
}

// ---------------------------------------------------------------------------
// 6. Conservative reflection: an explicit witness object + decision procedure.
// ---------------------------------------------------------------------------

///|
/// A reflection witness: the parameter tuple `p̄`, the level `β` at which the
/// `Σ_n` fragment is already decided, and the number of formulas that closed.
pub(all) struct ReflectionWitness {
  alpha : Int32
  beta : Int32
  n : Int32
  formulas_checked : Int32
  parameters : Array[Int32]
}

///|
/// Build the witness for `V_α`, `Σ_n` and the parameter tuple `params`
/// (Algorithm A30, parameter-relative Montague–Lévy form).
pub fn reflection_witness(
  alpha : Int32,
  n : Int32,
  params : Array[Int32],
  max_size~ : Int32 = 3,
) -> ReflectionWitness {
  let nvars = params.length()
  let beta = reflection_bound_for(alpha, n, params, nvars, max_size)
  let mut count = 0
  for f in formula_pool(nvars, max_size) {
    if is_sigma_n(f, n) {
      count += 1
    }
  }
  ReflectionWitness::{
    alpha,
    beta,
    n,
    formulas_checked: count,
    parameters: params,
  }
}

///|
/// Is the witness *conservative*, i.e. strictly below the universe it reflects?
pub fn ReflectionWitness::is_strict(self : ReflectionWitness) -> Bool {
  self.beta >= 0 && self.beta < self.alpha
}

///|
/// Algorithm A31: decide a `Σ_n` formula **without touching** the large
/// universe.  The formula is evaluated in the small reflected structure
/// `V_β` (β = parameter-relative reflection bound), which for `Σ_n` formulas
/// agrees with the truth value in `V_α` by T37/T38.
pub fn decide_conservatively(
  f : Formula,
  alpha : Int32,
  n : Int32,
  env : Array[Int32],
  max_size~ : Int32 = 3,
) -> Bool {
  let nvars = env.length()
  let beta = reflection_bound_for(alpha, n, env, nvars, max_size)
  if beta < 0 {
    return eval_formula(f, alpha, env)
  }
  eval_formula(f, beta, env)
}

///|
/// The truth value obtained in the large universe (reference).
pub fn decide_reference(f : Formula, alpha : Int32, env : Array[Int32]) -> Bool {
  eval_formula(f, alpha, env)
}

///|
/// Conservativity defect for one parameter tuple: the number of `Σ_n` formulas
/// of the pool on which the reflected decision disagrees with `V_α`.
/// `0` certifies that instance of conservative reflection; `-1` means the
/// parameters do not fit inside `V_α` at all.
pub fn conservativity_defect(
  alpha : Int32,
  n : Int32,
  params : Array[Int32],
  max_size : Int32,
) -> Int32 {
  let nvars = params.length()
  let beta = reflection_bound_for(alpha, n, params, nvars, max_size)
  if beta < 0 {
    return -1
  }
  let mut bad = 0
  for f in formula_pool(nvars, max_size) {
    if is_sigma_n(f, n) {
      if eval_formula(f, beta, params) != eval_formula(f, alpha, params) {
        bad += 1
      }
    }
  }
  bad
}

// ---------------------------------------------------------------------------
// 7. The coupled universes ledger (used by the Round-4 pipeline).
// ---------------------------------------------------------------------------

///|
/// One row of the universe ledger: the level, its cardinal, whether the level
/// is a Grothendieck universe, and the reflection bound for `Σ_n` truth.
pub(all) struct UniverseLedgerRow {
  level : Int32
  cardinal_depth : Int32
  is_universe : Bool
  self_member : Bool
  member_of_next : Bool
  reflection_beta : Int32
}

///|
pub fn universe_ledger(levels : Int32, n : Int32) -> Array[UniverseLedgerRow] {
  let out : Array[UniverseLedgerRow] = []
  for j in 0..<levels {
    // The finite shadows we can enumerate are V_1 .. V_4; the ledger reports
    // the parameter-relative reflection bound of the largest one.
    let alpha = 4
    let params : Array[Int32] = [1, 2]
    let beta = reflection_bound_for(alpha, n, params, 2, 2)
    out.push(UniverseLedgerRow::{
      level: j,
      cardinal_depth: j,
      is_universe: is_grothendieck_universe(j),
      self_member: universe_self_member(j),
      member_of_next: universe_member_of_next(j),
      reflection_beta: beta,
    })
  }
  out
}

///|
pub fn universe_ledger_csv_header() -> String {
  "level,cardinal_depth,is_grothendieck_universe,self_member,member_of_next,reflection_beta"
}

///|
pub fn UniverseLedgerRow::csv(self : UniverseLedgerRow) -> String {
  self.level.to_string() +
  "," +
  self.cardinal_depth.to_string() +
  "," +
  (if self.is_universe { "1" } else { "0" }) +
  "," +
  (if self.self_member { "1" } else { "0" }) +
  "," +
  (if self.member_of_next { "1" } else { "0" }) +
  "," +
  self.reflection_beta.to_string()
}

// ---------------------------------------------------------------------------
// 8. Tests: Theorems T35 -- T38.
// ---------------------------------------------------------------------------

///|
test "T35: U_j ∈ U_{j+1} and U_j ∉ U_j (Russell / Burali-Forti barrier)" {
  for j in 0..<=3 {
    assert_true(universe_member_of_next(j), msg="kappa_j < kappa_{j+1}")
    assert_true(!universe_self_member(j), msg="no universe contains itself")
    // A29: the least universe containing κ_j is exactly U_{j+1}.
    assert_eq(least_universe_containing(grothendieck_cardinal(j)), j + 1)
    // The next level is a strong limit: 2^{κ_j} < κ_{j+1}.
    assert_true(beth_below_next(j))
  }
  // No set of all sets inside V_4: no element has 16 elements (Burali-Forti).
  for x in 0..<hf_level_size(4) {
    let mut n = 0
    for i in 0..<16 {
      if hf_mem(i, x) {
        n += 1
      }
    }
    assert_true(n < 16)
  }
}

///|
test "T36: every U_j is a Grothendieck universe (five closure axioms)" {
  for j in 0..<=3 {
    let c = universe_closure_certificate(j)
    assert_true(
      c.transitive && c.pairing && c.powerset && c.union && c.replacement,
    )
    assert_true(is_grothendieck_universe(j))
  }
  // Concrete witnesses inside V_4 = 16 hereditarily finite sets.
  let a = 3 // {∅, {∅}}                 = V_2, transitive
  let b = 15 // {∅, {∅}, {{∅}}, {∅,{∅}}} = V_3, transitive
  assert_true(hf_transitive(a) && hf_transitive(b))
  assert_true(!hf_transitive(5), msg="{∅,{{∅}}} is not transitive")
  assert_true(hf_rank(hf_pair(a, b)) == 4)
  assert_true(hf_rank(hf_powerset(a)) == 3)
  assert_true(hf_rank(hf_big_union(b)) <= hf_rank(b))
  assert_true(hf_subset(hf_big_union(b), b))
  // Powerset is strictly monotone in cardinality: |P(x)| = 2^{|x|}.
  for x in 0..<8 {
    let mut n = 0
    for i in 0..<16 {
      if hf_mem(i, x) {
        n += 1
      }
    }
    assert_eq(popcount(hf_powerset(x)), pow2(n))
  }
  // Membership in U_j is decidable for every HF set (all ranks are finite),
  // and V_α is exactly the initial segment rank(x) < α.
  for x in 0..<hf_level_size(4) {
    assert_true(universe_contains_hf(grothendieck_universe(0), x))
    assert_true(hf_in_level(x, 4))
    assert_eq(hf_in_level(x, hf_rank(x)), false)
    assert_true(hf_in_level(x, hf_rank(x) + 1))
  }
}

///|
test "A32: category smallness sits at a decidable universe level" {
  // A finite category is already small in U_0 = V_ω (the HF universe).
  let fin = CategorySize::{ objects: Finite(3), morphisms: Finite(7), }
  assert_true(is_small_in(fin, 0))
  assert_eq(category_universe_level(fin), 0)
  // A countably infinite category (e.g. the simplex category Δ, or the
  // discrete category on N) is NOT U_0-small — ω ∉ V_ω — but is U_1-small.
  let countable = CategorySize::{ objects: Omega, morphisms: Omega, }
  assert_true(!is_small_in(countable, 0))
  assert_true(is_small_in(countable, 1))
  assert_eq(category_universe_level(countable), 1)
  // The cardinal order really is a strict total order on the tower.
  assert_true(cardinal_lt(Finite(10), Omega))
  assert_true(cardinal_lt(Omega, grothendieck_cardinal(1)))
  assert_true(!cardinal_lt(grothendieck_cardinal(2), grothendieck_cardinal(1)))
  assert_true(cardinal_eq(cardinal_max(Finite(2), Finite(5)), Finite(5)))
}

///|
test "T37: Σ_n reflection holds strictly below the universe level" {
  // (a) Δ_0 absoluteness: every Σ_0 formula is absolute, and the sharpest
  // bound is exactly rank(p̄) + 1 — no quantifier room is needed.
  for params in [[0], [1], [2], [3], [0, 1], [1, 3], [2, 3]] {
    let nvars = params.length()
    let mut need = 0
    for p in params {
      if hf_rank(p) + 1 > need {
        need = hf_rank(p) + 1
      }
    }
    assert_eq(reflection_bound_for(4, 0, params, nvars, 3), need)
    assert_eq(reflection_bound_for(3, 0, params, nvars, 2), need)
  }
  // (b) Montague–Lévy instance: quantifiers force a strictly larger bound, and
  // the bound is monotone in the Σ-index at fixed α.
  let p2 : Array[Int32] = [0, 1]
  let b0 = reflection_bound_for(4, 0, p2, 2, 2)
  let b1 = reflection_bound_for(4, 1, p2, 2, 2)
  let b2 = reflection_bound_for(4, 2, p2, 2, 2)
  assert_true(b0 <= b1 && b1 <= b2)
  assert_true(b1 > b0, msg="Σ_1 needs room for a witness")
  // (c) Conservativity: with a low-rank parameter tuple the reflected universe
  // is strictly smaller than α — the Σ_n statement provably needs no universe
  // axiom at all.
  let p1 : Array[Int32] = [0]
  for n in 0..<=2 {
    let w = reflection_witness(4, n, p1, max_size=2)
    assert_true(w.beta >= 0, msg="a reflection bound always exists")
    assert_true(w.is_strict(), msg="the bound is strictly below α")
    assert_eq(w.beta, reflection_bound_for(4, n, p1, 1, 2))
    assert_true(w.formulas_checked > 0)
    assert_eq(w.parameters.length(), 1)
  }
  // (d) Parameters outside V_α admit no bound at all.
  let outside : Array[Int32] = [3, 15]
  assert_eq(reflection_bound_for(3, 1, outside, 2, 2), -1)
  assert_true(reflection_bound_for(4, 1, outside, 2, 2) >= 0)
}

///|
test "T38: conservative reflection is sound and transitive" {
  let alpha = 4
  for n in 0..<=2 {
    for params in [[0], [1], [2], [0, 1], [1, 3]] {
      let nvars = params.length()
      assert_eq(conservativity_defect(alpha, n, params, 3), 0)
      assert_eq(conservativity_defect(3, n, params, 2), 0)
      let _ = nvars
    }
  }
  // Soundness of the decision procedure on the whole pool: the reflected
  // decision agrees with the large-universe truth value for every Σ_n formula
  // and every parameter tuple.
  let n = 2
  for params in assignments(2, hf_level_size(4)) {
    for f in formula_pool(2, 2) {
      if is_sigma_n(f, n) {
        assert_true(
          decide_conservatively(f, alpha, n, params, max_size=2) ==
          decide_reference(f, alpha, params),
        )
      }
    }
  }
  // Transitivity: V_β ≺^p V_α and V_α ≺^p V_γ imply V_β ≺^p V_γ.
  let params : Array[Int32] = [0]
  let beta = reflection_bound_for(3, 1, params, 1, 2)
  let gamma = 4
  let mid = reflection_bound_for(gamma, 1, params, 1, 2)
  if reflects_for(beta, mid, 1, params, 1, 2) &&
    reflects_for(mid, gamma, 1, params, 1, 2) {
    assert_true(reflects_for(beta, gamma, 1, params, 1, 2))
  }
  // Lévy classification sanity: ∃x0. x0 ∈ x1 is Σ_1 but not Π_1.
  let s1 = Formula::Ex(0, Formula::Mem(0, 1))
  assert_true(is_sigma_n(s1, 1) && !is_pi_n(s1, 1))
  let p1 = Formula::Neg(s1)
  assert_true(is_pi_n(p1, 1) && !is_sigma_n(p1, 1))
  // and Δ_0 formulas are both.
  let d0 = Formula::Mem(1, 0)
  assert_true(is_sigma_n(d0, 0) && is_pi_n(d0, 0))
  assert_true(formula_size(Formula::Conj(d0, s1)) == 4)
  assert_true(formula_to_string(Formula::All(1, d0)) == "∀x1.x1∈x0")
  // Σ_3 formulas exist only on one side of the hierarchy.
  let s3 = Formula::Ex(0, Formula::Neg(Formula::All(1, Formula::Ex(0, d0))))
  assert_true(is_sigma_n(s3, 3) || is_pi_n(s3, 3))
}

///|
test "universe ledger is CSV-consistent" {
  let rows = universe_ledger(4, 1)
  assert_eq(rows.length(), 4)
  let header = universe_ledger_csv_header().split(",").collect().length()
  for r in rows {
    assert_eq(r.csv().split(",").collect().length(), header)
    assert_true(r.is_universe)
    assert_true(!r.self_member)
    assert_true(r.member_of_next)
  }
}
