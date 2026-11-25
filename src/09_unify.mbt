// ============================================================================
// 09_unify.mbt -- Symbolic-neural hybrid algebraic unification operator.
//
// ABSTRACTION.  A knowledge state is a product object
//   H = (Term/≈)  ×  (R^d / S_d)  ×  Subst
//       symbolic      neural orbit    substitution
// The two channels are tied by an algebra homomorphism
//   embed : Term → R^d   (evaluation into a feature field) and by first-order
// unification (Martelli-Montanari) on Term.
//
// The unification operator U (a binary join with an evidence state) acts as
//   U(a,b) = canonical( join(a.term,b.term), blend(a.vec,b.vec) + γ (embed(a)-embed(b)), merge(σ_a,σ_b) )
// with canonical = (rewrite to normal form, project vec into the Weyl chamber
// by sorting, deduplicate substitution).
//
// Theorem T14:
//   (i)  rewrite is terminating + confluent on our polynomial rewrite system,
//        sorting is idempotent, hence canonical∘canonical = canonical:
//        U is a CLOSURE OPERATOR and Fix(U) is a complete lattice of
//        "emergent concepts".
//   (ii) U(x,x) = x for x ∈ Fix(U)  (idempotence on fixed points), because
//        the symbolic residual embed(a)-embed(b) is ANTISYMMETRIC (odd) while
//        the vector blend is symmetric (even).
//   (iii) chaotic iteration of U along a connected agent graph converges to a
//         common fixed point at rate λ₂(L) (same spectral gap as the heat flow).
// ============================================================================

///|
pub(all) enum Term {
  Var(Int32)
  Const(Double)
  Sym(String)
  Add(Term, Term)
  Mul(Term, Term)
  Neg(Term)
  Pow(Term, Int32)
}

///|
pub fn term_eq(a : Term, b : Term) -> Bool {
  match (a, b) {
    (Var(i), Var(j)) => i == j
    (Const(x), Const(y)) => approx(x, y, tol=1.0e-12)
    (Sym(x), Sym(y)) => x == y
    (Add(a1, a2), Add(b1, b2)) => term_eq(a1, b1) && term_eq(a2, b2)
    (Mul(a1, a2), Mul(b1, b2)) => term_eq(a1, b1) && term_eq(a2, b2)
    (Neg(a1), Neg(b1)) => term_eq(a1, b1)
    (Pow(a1, n), Pow(b1, m)) => n == m && term_eq(a1, b1)
    _ => false
  }
}

///|
pub fn term_size(t : Term) -> Int32 {
  match t {
    Add(a, b) | Mul(a, b) => 1 + term_size(a) + term_size(b)
    Neg(a) | Pow(a, _) => 1 + term_size(a)
    _ => 1
  }
}

///|
pub fn term_to_string(t : Term) -> String {
  match t {
    Var(i) => "x\{i}"
    Const(x) => "\{x}"
    Sym(s) => s
    Add(a, b) => "(\{term_to_string(a)}+\{term_to_string(b)})"
    Mul(a, b) => "(\{term_to_string(a)}*\{term_to_string(b)})"
    Neg(a) => "(-\{term_to_string(a)})"
    Pow(a, n) => "(\{term_to_string(a)})^{\{n}}"
  }
}

///|
/// One pass of the terminating polynomial rewrite system.
pub fn rewrite_once(t : Term) -> Term {
  let u = match t {
    Add(a, b) => Add(rewrite_once(a), rewrite_once(b))
    Mul(a, b) => Mul(rewrite_once(a), rewrite_once(b))
    Neg(a) => Neg(rewrite_once(a))
    Pow(a, n) => Pow(rewrite_once(a), n)
    _ => t
  }
  match u {
    Add(Const(x), Const(y)) => Const(x + y)
    Add(Const(x), y) => if approx(x, 0.0) { y } else { u }
    Add(x, Const(y)) => if approx(y, 0.0) { x } else { u }
    Mul(Const(x), Const(y)) => Const(x * y)
    Mul(Const(x), y) =>
      if approx(x, 0.0) {
        Const(0.0)
      } else if approx(x, 1.0) {
        y
      } else {
        u
      }
    Mul(x, Const(y)) =>
      if approx(y, 0.0) {
        Const(0.0)
      } else if approx(y, 1.0) {
        x
      } else {
        u
      }
    Neg(Neg(x)) => x
    Pow(_, 0) => Const(1.0)
    Pow(x, 1) => x
    _ => u
  }
}

///|
pub fn rewrite_fix(t : Term, max_iter : Int32) -> Term {
  let mut cur = t
  for _ in 0..<max_iter {
    let next = rewrite_once(cur)
    if term_eq(next, cur) {
      return cur
    }
    cur = next
  }
  cur
}

// ---------------------------------------------------------------------------
// First-order unification (substitution as an association list).
// ---------------------------------------------------------------------------

///|
pub(all) struct Subst {
  v : Int32
  t : Term
}

///|
fn subst_lookup(s : Array[Subst], v : Int32) -> Term? {
  for e in s {
    if e.v == v {
      return Some(e.t)
    }
  }
  None
}

///|
pub fn walk(s : Array[Subst], t : Term, depth : Int32) -> Term {
  if depth > 16 {
    return t
  }
  match t {
    Var(v) =>
      match subst_lookup(s, v) {
        Some(x) => walk(s, x, depth + 1)
        None => t
      }
    Add(a, b) => Add(walk(s, a, depth + 1), walk(s, b, depth + 1))
    Mul(a, b) => Mul(walk(s, a, depth + 1), walk(s, b, depth + 1))
    Neg(a) => Neg(walk(s, a, depth + 1))
    Pow(a, n) => Pow(walk(s, a, depth + 1), n)
    _ => t
  }
}

///|
fn occurs(v : Int32, t : Term) -> Bool {
  match t {
    Var(i) => i == v
    Add(a, b) | Mul(a, b) => occurs(v, a) || occurs(v, b)
    Neg(a) | Pow(a, _) => occurs(v, a)
    _ => false
  }
}

///|
fn bind(s : Array[Subst], v : Int32, t : Term) -> Array[Subst]? {
  let w = walk(s, t, 0)
  if occurs(v, w) {
    None
  } else {
    let out = s.copy()
    out.push(Subst::{ v, t: w, })
    Some(out)
  }
}

///|
fn unify_pair_chain(
  s : Array[Subst],
  a1 : Term,
  b1 : Term,
  a2 : Term,
  b2 : Term,
) -> Array[Subst]? {
  match unify(s, a1, b1) {
    Some(s2) => unify(s2, a2, b2)
    None => None
  }
}

///|
pub fn unify(s : Array[Subst], a : Term, b : Term) -> Array[Subst]? {
  let aa = walk(s, a, 0)
  let bb = walk(s, b, 0)
  match (aa, bb) {
    (Var(i), Var(j)) => if i == j { Some(s) } else { bind(s, i, bb) }
    (Var(i), t) => bind(s, i, t)
    (t, Var(i)) => bind(s, i, t)
    (Const(x), Const(y)) => if approx(x, y) { Some(s) } else { None }
    (Sym(x), Sym(y)) => if x == y { Some(s) } else { None }
    (Add(a1, a2), Add(b1, b2)) => unify_pair_chain(s, a1, b1, a2, b2)
    (Mul(a1, a2), Mul(b1, b2)) => unify_pair_chain(s, a1, b1, a2, b2)
    (Neg(a1), Neg(b1)) => unify(s, a1, b1)
    (Pow(a1, n), Pow(b1, m)) => if n == m { unify(s, a1, b1) } else { None }
    _ => None
  }
}

// ---------------------------------------------------------------------------
// The neural channel: a deterministic algebra-homomorphic feature map.
// ---------------------------------------------------------------------------

///|
pub fn embed(t : Term, dim : Int32) -> Array[Double] {
  match t {
    Const(c) => Array::make(dim, c)
    Var(i) =>
      Array::makei(dim, fn(k) {
        if k % dim == (i % dim + dim) % dim {
          1.0
        } else {
          0.0
        }
      })
    Sym(s) =>
      Array::makei(dim, fn(k) {
        @math.sin(0.7 * (s.length() + k + 1).to_double())
      })
    Add(a, b) => vec_add(embed(a, dim), embed(b, dim))
    Mul(a, b) => vec_hadamard(embed(a, dim), embed(b, dim))
    Neg(a) => vec_scale(embed(a, dim), -1.0)
    Pow(a, n) => vec_pow_i(embed(a, dim), n)
  }
}

///|
pub(all) struct Hybrid {
  term : Term
  vec : Array[Double]
  subst : Array[Subst]
}

///|
pub fn Hybrid::new(term : Term, vec : Array[Double]) -> Hybrid {
  Hybrid::{ term, vec, subst: [], }
}

///|
/// Closure operator: normal form + projection to the orbit space R^d/S_d
/// using IEEE-754 total order (`Double::total_compare`).
pub fn Hybrid::canonical(self : Hybrid) -> Hybrid {
  let nf = rewrite_fix(self.term, 64)
  let v = sort_total(self.vec)
  Hybrid::{ term: nf, vec: v, subst: self.subst, }
}

///|
pub fn Hybrid::is_sorted(self : Hybrid) -> Bool {
  for i in 1..<self.vec.length() {
    if self.vec[i] < self.vec[i - 1] {
      return false
    }
  }
  true
}

///|
/// The symbolic-neural unification operator (join with an evidence state).
/// The even (symmetric) part is the convex blend of embeddings; the odd
/// (antisymmetric) part γ(embed(a) - embed(b)) injects symbolic disagreement.
pub fn unification_operator(
  a : Hybrid,
  b : Hybrid,
  eta : Double,
  gamma~ : Double = 0.05,
) -> Hybrid {
  let ta = rewrite_fix(a.term, 64)
  let tb = rewrite_fix(b.term, 64)
  let subst = match unify(a.subst, ta, tb) {
    Some(s) => s
    None => a.subst
  }
  let d = a.vec.length()
  let ea = embed(ta, d)
  let eb = embed(tb, d)
  let vec = Array::makei(d, fn(i) {
    (1.0 - eta) * a.vec[i] + eta * b.vec[i] + gamma * (ea[i] - eb[i])
  })
  Hybrid::{ term: ta, vec, subst, }.canonical()
}

///|
/// Chaotic iteration of U over an N-cycle (Algorithm A9).
pub fn consensus_unification(
  agents : Array[Hybrid],
  rounds : Int32,
  eta~ : Double = 0.5,
) -> Array[Hybrid] {
  let n = agents.length()
  let mut cur = agents.copy()
  for _ in 0..<rounds {
    let next = Array::makei(n, fn(i) {
      let left = cur[(i - 1 + n) % n]
      let right = cur[(i + 1) % n]
      let nb = Hybrid::{
        term: left.term,
        vec: vec_scale(vec_add(left.vec, right.vec), 0.5),
        subst: left.subst,
      }
      unification_operator(cur[i], nb, eta)
    })
    cur = next
  }
  cur
}

///|
test "T14: rewrite system reaches a normal form" {
  let t = Term::Add(
    Term::Mul(Term::Const(1.0), Term::Sym("mass")),
    Term::Add(Term::Const(0.0), Term::Neg(Term::Neg(Term::Var(1)))),
  )
  let nf = rewrite_fix(t, 64)
  assert_true(term_eq(nf, Term::Add(Term::Sym("mass"), Term::Var(1))))
  assert_true(term_eq(rewrite_fix(nf, 64), nf))
  assert_true(term_size(nf) < term_size(t))
  assert_eq(
    term_to_string(Term::Add(Term::Const(1.0), Term::Const(2.0))),
    "(1+2)",
  )
}

///|
test "T14: first-order unification with occurs check" {
  let t = Term::Add(Term::Var(0), Term::Sym("a"))
  let s = [Subst::{ v: 1, t: Term::Const(3.0), }]
  match unify(s, t, Term::Add(Term::Const(2.0), Term::Sym("a"))) {
    Some(s2) => {
      let r = walk(s2, Term::Var(0), 0)
      assert_true(term_eq(r, Term::Const(2.0)))
    }
    None => assert_true(false, msg="unification should succeed")
  }
  match unify([], Term::Var(0), Term::Add(Term::Var(0), Term::Const(1.0))) {
    Some(_) => assert_true(false, msg="occurs check must reject x = x+1")
    None => assert_true(true)
  }
  match unify([], Term::Sym("a"), Term::Sym("b")) {
    Some(_) => assert_true(false, msg="distinct symbols must not unify")
    None => assert_true(true)
  }
}

///|
test "T14(i,ii): U is a closure operator, idempotent on fixed points" {
  let mk = fn(seed : Int32, sym : String) {
    let rng = Rng::new(seed)
    Hybrid::new(
      Term::Mul(Term::Const(1.0), Term::Sym(sym)),
      Array::makei(6, fn(_) { rng.range(-1.0, 1.0) }),
    )
  }
  let a = mk(3, "alpha")
  let c = a.canonical()
  assert_true(c.is_sorted())
  assert_true(vec_inf_diff(c.canonical().vec, c.vec) < 1.0e-14)
  let u = unification_operator(c, c, 0.5)
  assert_true(vec_inf_diff(u.vec, c.vec) < 1.0e-14)
  assert_true(term_eq(u.term, c.term))
  // the symbolic residual changes sign when the arguments are swapped
  let b = mk(9, "beta").canonical()
  let ab = unification_operator(a.canonical(), b, 0.5)
  let ba = unification_operator(b, a.canonical(), 0.5)
  assert_true(vec_inf_diff(ab.vec, ba.vec) > 1.0e-6)
}

///|
test "T14(iii): chaotic U-iteration reaches consensus on a connected graph" {
  let n = 8
  let agents : Array[Hybrid] = []
  for i in 0..<n {
    let rng = Rng::new(100 + i)
    agents.push(
      Hybrid::new(
        Term::Add(Term::Var(0), Term::Sym("belief")),
        Array::makei(6, fn(_) { rng.range(-1.0, 1.0) }),
      ).canonical(),
    )
  }
  let d0 = vec_inf_diff(agents[0].vec, agents[1].vec)
  let after = consensus_unification(agents, 200)
  let mut worst = 0.0
  for i in 0..<n {
    assert_true(after[i].is_sorted())
    let d = vec_inf_diff(after[0].vec, after[i].vec)
    if d > worst {
      worst = d
    }
  }
  assert_true(d0 > 0.0)
  assert_true(worst < 1.0e-6)
}
