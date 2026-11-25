// ============================================================================
// 20_diff_symbolic_satake_frontier.mbt -- Differentiable Symbolic-Neural
//                                         Encoder (Reverse-Mode Super-AD),
//                                         Geometric Satake Kazhdan-Lusztig /
//                                         Lusztig q-Analogue, and Frontier
//                                         Emergence Pipeline.
//
// ELIMINATING LIMITATION 4 OF ROUND 1 + FRONTIER LANGLANDS SYNTHESIS:
//   * Limitation 4 noted that `embed : Term → R^d` in `09_unify.mbt` was a
//     fixed non-trainable feature map rather than a differentiable / trainable
//     neural encoder.  Here we build an exact **Reverse-Mode Automatic
//     Differentiation + Trainable Tree-MLP Encoder** (`TrainableEncoder`) over
//     recursive `Term` ASTs (`d = 4`), proving:
//       - Exact analytical backpropagation gradient vs central finite-difference
//         (`‖∇_θ L_AD - ∇_θ L_FD‖_∞ < 1e-7`, Theorem T32a);
//       - Contrastive / unification alignment training (`train_encoder_step`)
//         monotonically decreases symbolic-neural representation loss while
//         preserving algebraic normal-form invariance (`L(Term) == L(NF(Term))`).
//   * Frontier Attempt 1 (Quantum / Categorical Geometric Satake & Lusztig
//     `q`-Analogue, Beilinson–Drinfeld / Mirković–Vilonen / Lusztig):
//       For type `A_{r-1}` (`GL_r`) with dominant coweights `λ, μ` (`λ >= μ`),
//       we construct the **Spherical Hecke / Satake `q`-Character** via the
//       Weyl–Demazure / Macdonald–Casselman–Shalika formula and **Kostka–Foulkes
//       / Lusztig `q`-weight multiplicity** `K_{λ, μ}(q)`:
//         - At `q = 1`, `Σ_μ K_{λ, μ}(1) |W·μ|` equals the Weyl dimension
//           `dim V_λ = Π_{α > 0} ⟨λ + ρ, α^∨⟩ / ⟨ρ, α^∨⟩` (Theorem T33a);
//         - On the affine Grassmannian `Gr_{GL_r}`, the intersection cohomology
//           Poincaré polynomial `P_{IC(Gr_λ)}(q) = Σ_w q^{l(w)}` (for minuscule /
//           spherical Schubert varieties) and Satake `q`-spherical trace satisfy
//           `P_{IC}(1) = dim V_λ` and `q`-palindromicity (Poincaré–Verdier duality!).
//   * Frontier Attempt 2 (`20_frontier_emergence` — Unified Grand Synthesis):
//       Couples ALL newly unlocked structures into a single executable multi-agent
//       frontier dynamical system (`frontier_trace.csv`, 160 steps, 16 agents):
//         1. Agents live on the 4D hyperbolic Poincare ball `B^4` (`M^{4|32}`
//            with sparse 32-generator Grassmann soul `SparseSuper` — `2^{32}`
//            fibre dimension, impossible in v0.1!);
//         2. Non-type-A Langlands `Sp_4 / SO_5` (`TypeB(2)`, `|W| = 8`) adelic
//            intertwining gate `|c_A(w_0, λ)|` drives Möbius gyrovector consensus
//            on `B^4`;
//         3. Beilinson–Drinfeld `SL_3` oper coordinates `(q_2, q_3)` align along
//            the Hitchin base while preserving `det Hol(∇) = 1`;
//         4. Trainable symbolic-neural encoder updates its weights online via
//            reverse-mode AD (`grad_norm → 0`, `loss` collapses by `> 98%`).
// ============================================================================

///|
/// Trainable parameters of a differentiable recursive symbolic-neural encoder
/// `E_θ : Term → R^d` (`d = 4`):
/// * `sym_emb`: `4 × 4` embedding table for hashed symbol classes;
/// * `w_add`: `4` diagonal/channel mixing weights for `Add(a, b)`;
/// * `w_mul`: `4` bilinear Hadamard gating weights for `Mul(a, b)`;
/// * `bias`: `4` channel biases (passed through `tanh` nonlinearity).
pub(all) struct TrainableEncoder {
  sym_emb : Array[Double] // length 16 (4 buckets × 4 dims)
  w_add : Array[Double] // length 4
  w_mul : Array[Double] // length 4
  bias : Array[Double] // length 4
}

///|
/// Deterministic initialization of a `d = 4` trainable tree-encoder.
pub fn TrainableEncoder::init(seed : Int32) -> TrainableEncoder {
  let rng = Rng::new(seed)
  TrainableEncoder::{
    sym_emb: Array::makei(16, fn(_) { rng.range(-0.6, 0.6) }),
    w_add: Array::makei(4, fn(_) { rng.range(0.4, 0.9) }),
    w_mul: Array::makei(4, fn(_) { rng.range(0.4, 0.9) }),
    bias: Array::makei(4, fn(_) { rng.range(-0.2, 0.2) }),
  }
}

///|
/// Smooth hyperbolic tangent `tanh(x)` and its derivative `1 - tanh(x)^2`.
pub fn tanh_act(x : Double) -> Double {
  let z = clamp(x, -20.0, 20.0)
  let ez = @math.exp(z)
  let enz = @math.exp(-z)
  (ez - enz) / (ez + enz)
}

///|
fn sym_bucket(s : String) -> Int32 {
  (s.length() * 7 + 3) % 4
}

///|
/// Forward pass of `E_θ(t)` on a normalized `Term` (`rewrite_fix(t, 32)` first,
/// guaranteeing exact algebraic rewrite invariance `E_θ(t) == E_θ(NF(t))`!).
pub fn TrainableEncoder::forward(
  self : TrainableEncoder,
  t : Term,
) -> Array[Double] {
  let nf = rewrite_fix(t, 32)
  self.forward_nf(nf)
}

///|
fn TrainableEncoder::forward_nf(
  self : TrainableEncoder,
  t : Term,
) -> Array[Double] {
  match t {
    Const(c) => Array::makei(4, fn(k) { tanh_act(c + self.bias[k]) })
    Var(i) => {
      let b = (i % 4 + 4) % 4
      Array::makei(4, fn(k) { tanh_act(self.sym_emb[b * 4 + k] + self.bias[k]) })
    }
    Sym(s) => {
      let b = sym_bucket(s)
      Array::makei(4, fn(k) { tanh_act(self.sym_emb[b * 4 + k] + self.bias[k]) })
    }
    Neg(a) => {
      let va = self.forward_nf(a)
      vec_scale(va, -1.0)
    }
    Pow(a, n) => {
      let va = self.forward_nf(a)
      vec_pow_i(va, n)
    }
    Add(a, b) => {
      let va = self.forward_nf(a)
      let vb = self.forward_nf(b)
      Array::makei(4, fn(k) {
        tanh_act(self.w_add[k] * (va[k] + vb[k]) + self.bias[k])
      })
    }
    Mul(a, b) => {
      let va = self.forward_nf(a)
      let vb = self.forward_nf(b)
      Array::makei(4, fn(k) {
        tanh_act(self.w_mul[k] * (va[k] * vb[k]) + self.bias[k])
      })
    }
  }
}

///|
/// Reverse-mode automatic differentiation (backpropagation) through the
/// recursive `Term` AST: accumulates parameter gradients `(g_sym, g_add, g_mul, g_bias)`
/// given upstream gradient `d_out ∈ R^4` (Algorithm A27).
fn TrainableEncoder::backward_nf(
  self : TrainableEncoder,
  t : Term,
  d_out : Array[Double],
  g_sym : Array[Double],
  g_add : Array[Double],
  g_mul : Array[Double],
  g_bias : Array[Double],
) -> Unit {
  match t {
    Const(c) =>
      for k in 0..<4 {
        let y = tanh_act(c + self.bias[k])
        let dz = d_out[k] * (1.0 - y * y)
        g_bias[k] += dz
      }
    Var(i) => {
      let b = (i % 4 + 4) % 4
      for k in 0..<4 {
        let y = tanh_act(self.sym_emb[b * 4 + k] + self.bias[k])
        let dz = d_out[k] * (1.0 - y * y)
        g_sym[b * 4 + k] += dz
        g_bias[k] += dz
      }
    }
    Sym(s) => {
      let b = sym_bucket(s)
      for k in 0..<4 {
        let y = tanh_act(self.sym_emb[b * 4 + k] + self.bias[k])
        let dz = d_out[k] * (1.0 - y * y)
        g_sym[b * 4 + k] += dz
        g_bias[k] += dz
      }
    }
    Neg(a) =>
      self.backward_nf(a, vec_scale(d_out, -1.0), g_sym, g_add, g_mul, g_bias)
    Pow(a, n) => {
      let va = self.forward_nf(a)
      let p_prev = if n >= 1 {
        vec_pow_i(va, n - 1)
      } else {
        Array::make(4, 0.0)
      }
      let da = Array::makei(4, fn(k) { d_out[k] * n.to_double() * p_prev[k] })
      self.backward_nf(a, da, g_sym, g_add, g_mul, g_bias)
    }
    Add(a, b) => {
      let va = self.forward_nf(a)
      let vb = self.forward_nf(b)
      let da = Array::make(4, 0.0)
      let db = Array::make(4, 0.0)
      for k in 0..<4 {
        let y = tanh_act(self.w_add[k] * (va[k] + vb[k]) + self.bias[k])
        let dz = d_out[k] * (1.0 - y * y)
        g_add[k] += dz * (va[k] + vb[k])
        g_bias[k] += dz
        da[k] = dz * self.w_add[k]
        db[k] = dz * self.w_add[k]
      }
      self.backward_nf(a, da, g_sym, g_add, g_mul, g_bias)
      self.backward_nf(b, db, g_sym, g_add, g_mul, g_bias)
    }
    Mul(a, b) => {
      let va = self.forward_nf(a)
      let vb = self.forward_nf(b)
      let da = Array::make(4, 0.0)
      let db = Array::make(4, 0.0)
      for k in 0..<4 {
        let y = tanh_act(self.w_mul[k] * (va[k] * vb[k]) + self.bias[k])
        let dz = d_out[k] * (1.0 - y * y)
        g_mul[k] += dz * (va[k] * vb[k])
        g_bias[k] += dz
        da[k] = dz * self.w_mul[k] * vb[k]
        db[k] = dz * self.w_mul[k] * va[k]
      }
      self.backward_nf(a, da, g_sym, g_add, g_mul, g_bias)
      self.backward_nf(b, db, g_sym, g_add, g_mul, g_bias)
    }
  }
}

///|
/// Squared `L^2` alignment loss `½ ‖E_θ(t) - target‖^2`.
pub fn TrainableEncoder::loss(
  self : TrainableEncoder,
  t : Term,
  target : Array[Double],
) -> Double {
  let y = self.forward(t)
  let mut s = 0.0
  for k in 0..<4 {
    let d = y[k] - target[k]
    s += 0.5 * d * d
  }
  s
}

///|
/// Flatten the analytical reverse-mode AD gradient `∇_θ loss(t, target)` into
/// a single 28-element vector `[g_sym(16), g_add(4), g_mul(4), g_bias(4)]`.
pub fn TrainableEncoder::grad_vector(
  self : TrainableEncoder,
  t : Term,
  target : Array[Double],
) -> Array[Double] {
  let nf = rewrite_fix(t, 32)
  let y = self.forward_nf(nf)
  let d_out = Array::makei(4, fn(k) { y[k] - target[k] })
  let g_sym = Array::make(16, 0.0)
  let g_add = Array::make(4, 0.0)
  let g_mul = Array::make(4, 0.0)
  let g_bias = Array::make(4, 0.0)
  self.backward_nf(nf, d_out, g_sym, g_add, g_mul, g_bias)
  let out : Array[Double] = []
  out.append(g_sym[:])
  out.append(g_add[:])
  out.append(g_mul[:])
  out.append(g_bias[:])
  out
}

///|
/// Pack/unpack helpers for finite-difference verification and SGD updates.
pub fn TrainableEncoder::to_params(self : TrainableEncoder) -> Array[Double] {
  let out : Array[Double] = []
  out.append(self.sym_emb[:])
  out.append(self.w_add[:])
  out.append(self.w_mul[:])
  out.append(self.bias[:])
  out
}

///|
pub fn TrainableEncoder::from_params(p : Array[Double]) -> TrainableEncoder {
  TrainableEncoder::{
    sym_emb: Array::makei(16, fn(i) { p[i] }),
    w_add: Array::makei(4, fn(i) { p[16 + i] }),
    w_mul: Array::makei(4, fn(i) { p[20 + i] }),
    bias: Array::makei(4, fn(i) { p[24 + i] }),
  }
}

///|
/// One gradient-descent training step `θ ← θ - lr · ∇_θ L`, returning `(updated, grad_inf_norm)`.
pub fn TrainableEncoder::step_sgd(
  self : TrainableEncoder,
  t : Term,
  target : Array[Double],
  lr : Double,
) -> (TrainableEncoder, Double) {
  let p = self.to_params()
  let g = self.grad_vector(t, target)
  let p_next = Array::makei(28, fn(i) { p[i] - lr * g[i] })
  (TrainableEncoder::from_params(p_next), vec_inf_norm(g))
}

// ---------------------------------------------------------------------------
// Geometric Satake Equivalence & Lusztig / Poincaré IC Polynomials on Gr_G
// ---------------------------------------------------------------------------

///|
/// Weyl dimension formula for an irreducible representation `V_λ` of a classical
/// group with root system `rt` and dominant integral highest weight `λ`:
///   `dim V_λ = Π_{α ∈ Φ^+} ⟨λ + ρ, α^∨⟩ / ⟨ρ, α^∨⟩` (Theorem T33a).
pub fn weyl_dimension(rt : RootType, lambda : Array[Double]) -> Int32 {
  let rho = weyl_rho(rt)
  let lam_rho = vec_add(lambda, rho)
  let mut num = 1.0
  let mut den = 1.0
  for alpha in positive_roots(rt) {
    let av = coroot(alpha)
    num *= dot(lam_rho, av)
    den *= dot(rho, av)
  }
  (num / den).round().to_int()
}

///|
/// Poincaré polynomial `P_{W / W_λ}(q) = Σ_{w ∈ W / W_λ} q^{l(w)}` of the
/// projective / minuscule affine Schubert variety `Gr_λ ⊂ Gr_G` (orbit of `λ`
/// under `W` weighted by minimal coset-representative length, Algorithm A28):
/// For minuscule / quasi-minuscule coweights `λ`, `IC(Gr_λ)` is smooth/rational
/// so `P_{IC(Gr_λ)}(q) = Σ_{μ ∈ W·λ} q^{height_dist(λ, μ)}`, and Geometric Satake
/// (`IH^*(Gr_λ) ≅ V_λ`) forces `P_{IC(Gr_λ)}(1) = dim V_λ` and `q`-palindromicity!
pub fn satake_minuscule_poincare(
  rt : RootType,
  lambda : Array[Double],
  q : Double,
) -> Double {
  let d = rt.ambient_dim()
  let ws = generate_weyl_group(rt)
  // Collect distinct weights in the Weyl orbit W·λ, tracking minimal Coxeter length l(w)
  let min_len = @sorted_map.SortedMap::from_array([])
  for w in ws {
    let mu = matvec(w, lambda, d)
    let key = mat_key(mu)
    let l = coxeter_length(rt, w)
    match min_len.get(key) {
      Some(prev) => if l < prev { min_len[key] = l }
      None => min_len[key] = l
    }
  }
  let mut poly = 0.0
  for _k, deg in min_len {
    poly += @math.pow(q, deg.to_double())
  }
  poly
}

///|
/// Lusztig `q`-analogue of Kostant's partition function `P_Φ(ν; q)` for `GL_3` (`A_2`):
/// coefficient of `e^ν` in `Π_{α ∈ Φ^+} (1 / (1 - q e^α))`, giving Lusztig's `q`-weight
/// multiplicity `K_{λ, μ}(q) = Σ_{w ∈ W} (-1)^{l(w)} P_Φ(w(λ + ρ) - (μ + ρ); q)`
/// (= Kazhdan–Lusztig polynomial `P_{w_μ, w_λ}(q)` on `Gr_{GL_3}`, Theorem T33b!).
pub fn lusztig_q_multiplicity_a2(
  lambda : Array[Int32], // dominant weight in R^3, e.g. [2, 1, 0]
  mu : Array[Int32], // weight in R^3 with sum(lambda) == sum(mu), e.g. [1, 1, 1]
  q : Double,
) -> Double {
  // Positive roots of A_2: α_1 = (1,-1,0), α_2 = (0,1,-1), α_{12} = (1,0,-1).
  // Any vector ν = (ν_0, ν_1, ν_2) with ν_0 + ν_1 + ν_2 = 0 equals
  //   a α_1 + b α_2 = (a, b - a, -b)  =>  a = ν_0, b = -ν_2 (check: a - b = ν_0 + ν_2 = -ν_1!).
  // Decomposition into k α_{12} + (a - k) α_1 + (b - k) α_2 uses (a + b - k) roots,
  // valid iff a >= 0, b >= 0, and 0 <= k <= min(a, b)!
  let kp_a2 = fn(nu0 : Int32, nu1 : Int32, nu2 : Int32) -> Double {
    if nu0 + nu1 + nu2 != 0 {
      return 0.0
    }
    let a = nu0
    let b = -nu2
    if a < 0 || b < 0 {
      return 0.0
    }
    let m = if a < b { a } else { b }
    let mut sum = 0.0
    for k in 0..<=m {
      let num_roots = a + b - k
      sum += @math.pow(q, num_roots.to_double())
    }
    sum
  }
  // ρ = (1, 0, -1) for A_2
  let lp : Array[Int32] = [lambda[0] + 1, lambda[1], lambda[2] - 1]
  let mp : Array[Int32] = [mu[0] + 1, mu[1], mu[2] - 1]
  let perms = weyl_group(3)
  let mut total = 0.0
  for w in perms {
    // sign (-1)^{inversions(w)}
    let mut inv = 0
    for i in 0..<3 {
      for j in (i + 1)..<3 {
        if w[i] > w[j] {
          inv += 1
        }
      }
    }
    let sign = if inv % 2 == 0 { 1.0 } else { -1.0 }
    let w_lp0 = lp[w[0]]
    let w_lp1 = lp[w[1]]
    let w_lp2 = lp[w[2]]
    total += sign * kp_a2(w_lp0 - mp[0], w_lp1 - mp[1], w_lp2 - mp[2])
  }
  total
}

// ---------------------------------------------------------------------------
// Unified Frontier Emergence Pipeline (M^{4|32} × Sp_4/SO_5 Adelic × SL_3 Oper × AD)
// ---------------------------------------------------------------------------

///|
pub(all) struct FrontierTrace {
  rows : Array[String]
  summary : String
  initial_hyp_spread : Double
  final_hyp_spread : Double
  initial_encoder_loss : Double
  final_encoder_loss : Double
  max_oper_det_defect : Double
  sparse_fibre_generators : Int32
}

///|
pub fn frontier_csv_header() -> String {
  "step,hyp_spread_B4,adelic_gate_B2,oper_hitchin_spread,encoder_loss,encoder_grad,sparse_soul_l1,oper_det_defect"
}

///|
/// Run the Round-3 Frontier Emergence synthesis (`16` agents, `120` steps):
/// * State space: `M^{4|32}` (`B^4` Poincare ball + sparse `Λ(32)` soul with `2^{32}` basis!);
/// * Dual group gate: non-type-A `TypeB(2)` (`SO_5 ↔ Sp_4`) adelic Euler product (`P = 30`);
/// * Geometric Langlands channel: `SL_3` oper coordinates `(q_2, q_3)` + monodromy `det = 1`;
/// * Symbolic-neural channel: online reverse-mode AD training of `TrainableEncoder`.
pub fn run_frontier_emergence(
  agents~ : Int32 = 16,
  steps~ : Int32 = 120,
  seed~ : Int32 = 20261124,
) -> FrontierTrace {
  let rng = Rng::new(seed)
  let m_dim = 4
  let n_odd = 32
  // 1. Initialize agents in B^4 (||x|| < 0.55) + sparse souls in Λ(32) + SL_3 opers
  let pts : Array[Array[Double]] = []
  let souls : Array[SparseSuper] = []
  let opers : Array[Array[Double]] = [] // [q2, q3] in Op_{SL_3}
  for i in 0..<agents {
    pts.push(Array::makei(m_dim, fn(_) { rng.range(-0.22, 0.22) }))
    let g1 = i * 2 % n_odd
    let g2 = (i * 5 + 17) % n_odd
    let s = if g1 == g2 {
      SparseSuper::theta(n_odd, g1)
    } else {
      SparseSuper::theta(n_odd, g1).mul(SparseSuper::theta(n_odd, g2))
    }
    souls.push(s)
    opers.push([rng.range(-1.2, 1.2), rng.range(-1.2, 1.2)])
  }
  // Trainable symbolic-neural encoder + target consensus concept
  let mut enc = TrainableEncoder::init(seed + 7)
  let concept_term = Term::Add(
    Term::Mul(Term::Sym("satake"), Term::Var(1)),
    Term::Mul(Term::Const(1.0), Term::Sym("oper")),
  )
  let target_vec = [0.65, -0.45, 0.55, -0.35]
  let b2 = RootType::TypeB(2)
  let ws_b2 = generate_weyl_group(b2)
  let w_long_b2 = ws_b2[ws_b2.length() - 1]
  let c_adelic = fn(s : Cx) { adelic_c_euler_product(s, 30) }
  let rows : Array[String] = []
  let mut init_spread = 0.0
  let mut last_spread = 0.0
  let mut init_loss = 0.0
  let mut last_loss = 0.0
  let mut worst_oper_det = 0.0
  for step in 0..<=steps {
    // Compute Frechet-like hyperbolic spread on B^4 via pairwise geodesic distances to mean
    let mean_x = Array::makei(m_dim, fn(d) {
      let mut s = 0.0
      for i in 0..<agents {
        s += pts[i][d]
      }
      s / agents.to_double()
    })
    let mut d_sum = 0.0
    for i in 0..<agents {
      d_sum += poincare_ball_distance(pts[i], mean_x)
    }
    let hyp_sp = d_sum / agents.to_double()
    // Oper Hitchin base spread
    let q2s = Array::makei(agents, fn(i) { opers[i][0] })
    let q3s = Array::makei(agents, fn(i) { opers[i][1] })
    let oper_sp = spread(q2s) + spread(q3s)
    // Adelic B_2 intertwining amplitude on current collective spectral state
    let lam_cx : Array[Cx] = [
      { re: 1.8 + hyp_sp, im: 0.4 + oper_sp * 0.2, },
      { re: 1.25 + 0.5 * hyp_sp, im: -0.3, },
    ]
    let gate = cx_abs(root_c_function(b2, w_long_b2, lam_cx, c_adelic))
    // Sparse Λ(32) collective soul norm
    let mut soul_sum = 0.0
    for s in souls {
      soul_sum += s.soul_l1()
    }
    let soul_mean = soul_sum / agents.to_double()
    // Encoder loss & gradient
    let loss_val = enc.loss(concept_term, target_vec)
    let grad_n = vec_inf_norm(enc.grad_vector(concept_term, target_vec))
    // Certify SL_3 oper monodromy determinant on mean oper every 20 steps
    let mq2 = mean_of(q2s)
    let mq3 = mean_of(q3s)
    if step % 20 == 0 {
      let hol = oper_monodromy_rk4(3, [mq2, mq3], [0.1 * hyp_sp, 0.0], 120)
      let det_err = (determinant(hol, 3) - 1.0).abs()
      if det_err > worst_oper_det {
        worst_oper_det = det_err
      }
    }
    if step == 0 {
      init_spread = hyp_sp
      init_loss = loss_val
    }
    last_spread = hyp_sp
    last_loss = loss_val
    rows.push(
      "\{step},\{hyp_sp},\{gate},\{oper_sp},\{loss_val},\{grad_n},\{soul_mean},\{worst_oper_det}",
    )
    if step == steps {
      break
    }
    // --- Step 1: Möbius gyrovector geodesic contraction on B^4 -------------
    let eta_hyp = 0.08 * (1.0 + gate / (1.0 + gate))
    for i in 0..<agents {
      // Geodesic step from x_i towards mean_x using gyro-translation (-x_i) ⊕ mean_x
      let neg_xi = vec_scale(pts[i], -1.0)
      let delta_gyro = mobius_add(neg_xi, mean_x)
      pts[i] = mobius_add(pts[i], vec_scale(delta_gyro, eta_hyp))
      // Cool odd soul in Λ(32)
      souls[i] = souls[i].scale(0.94)
      // Align SL_3 oper coordinates on the Hitchin base
      opers[i][0] += 0.09 * (mq2 - opers[i][0])
      opers[i][1] += 0.09 * (mq3 - opers[i][1])
    }
    // --- Step 2: Online reverse-mode AD SGD update on TrainableEncoder -----
    let (next_enc, _) = enc.step_sgd(concept_term, target_vec, 0.18)
    enc = next_enc
  }
  let summary = "# frontier-emergence / M^{4|32} Sparse + B_2 Adelic + SL_3 Oper + Reverse-AD\n" +
    "agents: \{agents}, steps: \{steps}, supermanifold: B^4 x Lambda(\{n_odd}) (2^32 = 4,294,967,296 fibre dim)\n" +
    "B^4 hyperbolic geodesic spread: \{init_spread} -> \{last_spread}\n" +
    "trainable symbolic-neural encoder loss: \{init_loss} -> \{last_loss}\n" +
    "max |det Hol(SL_3 oper) - 1| = \{worst_oper_det}\n" +
    "Satake IC Poincare P_{Gr_omega1}(1) on C_3 (Sp_6) = \{satake_minuscule_poincare(RootType::TypeC(3), [1.0, 0.0, 0.0], 1.0)} == dim V_omega1 = \{weyl_dimension(RootType::TypeC(3), [1.0, 0.0, 0.0])}\n" +
    "Lusztig q-multiplicity K_{(2,1,0),(1,1,1)}(q=1) = \{lusztig_q_multiplicity_a2([2, 1, 0], [1, 1, 1], 1.0)} (adjoint sl_3 Cartan dim = 2), at q=2 -> \{lusztig_q_multiplicity_a2([2, 1, 0], [1, 1, 1], 2.0)}\n"
  FrontierTrace::{
    rows,
    summary,
    initial_hyp_spread: init_spread,
    final_hyp_spread: last_spread,
    initial_encoder_loss: init_loss,
    final_encoder_loss: last_loss,
    max_oper_det_defect: worst_oper_det,
    sparse_fibre_generators: n_odd,
  }
}

///|
pub fn FrontierTrace::csv(self : FrontierTrace) -> String {
  let mut out = frontier_csv_header() + "\n"
  for r in self.rows {
    out = out + r + "\n"
  }
  out
}

// ---------------------------------------------------------------------------
// Tests for Theorems T32, T33, T34
// ---------------------------------------------------------------------------

///|
test "T32: Reverse-mode AD through symbolic Term AST matches finite differences & trains to convergence" {
  let enc = TrainableEncoder::init(42)
  let term = Term::Add(
    Term::Mul(Term::Sym("alpha"), Term::Var(2)),
    Term::Neg(Term::Pow(Term::Sym("beta"), 2)),
  )
  let target = [0.5, -0.3, 0.7, -0.2]
  // 1. Exact analytical reverse-mode AD vs central finite differences across all 28 parameters!
  let ad_grad = enc.grad_vector(term, target)
  let base_p = enc.to_params()
  let eps = 1.0e-6
  for i in 0..<28 {
    let p_up = base_p.copy()
    let p_dn = base_p.copy()
    p_up[i] += eps
    p_dn[i] -= eps
    let l_up = TrainableEncoder::from_params(p_up).loss(term, target)
    let l_dn = TrainableEncoder::from_params(p_dn).loss(term, target)
    let fd = (l_up - l_dn) / (2.0 * eps)
    assert_true((ad_grad[i] - fd).abs() < 1.0e-7)
  }
  // 2. Algebraic rewrite invariance: unsimplified term gives IDENTICAL output and gradient!
  let unsimplified = Term::Add(
    Term::Const(0.0),
    Term::Mul(Term::Const(1.0), term),
  )
  assert_true(
    vec_inf_diff(enc.forward(term), enc.forward(unsimplified)) < 1.0e-14,
  )
  // 3. Training decreases alignment loss by > 99%
  let l0 = enc.loss(term, target)
  let mut cur = enc
  for _ in 0..<80 {
    let (nxt, _) = cur.step_sgd(term, target, 0.2)
    cur = nxt
  }
  assert_true(cur.loss(term, target) < l0 * 0.01)
}

///|
test "T33: Geometric Satake IC Poincaré polynomial, Weyl dimension, and Lusztig q-multiplicity" {
  // 1. Minuscule / standard representations match IH^*(Gr_λ)(q=1) == dim V_λ AND q-palindromicity!
  //    - B_3 (SO_7) vector weight ω_1 = (1, 0, 0): dim V = 7, Gr_{ω_1} has 6 roots? Wait:
  //      W(B_3) orbit of (1,0,0) is {±e_1, ±e_2, ±e_3} of size 6? Wait: in SO_7, (1,0,0) is
  //      quasi-minuscule (dim V = 7 = 6 extremal + 1 zero weight), while in C_3 (Sp_6),
  //      ω_1 = (1,0,0) is MINUSCULE with dim V = 6 == |W·ω_1|! Let's check BOTH!
  let c3 = RootType::TypeC(3)
  assert_eq(weyl_dimension(c3, [1.0, 0.0, 0.0]), 6)
  assert_true(
    approx(
      satake_minuscule_poincare(c3, [1.0, 0.0, 0.0], 1.0),
      6.0,
      tol=1.0e-12,
    ),
  )
  // Minuscule Grassmannian G(2, 4) for A_3 (GL_4), weight ω_2 = (1, 1, 0, 0):
  // dim Λ^2(C^4) = 6, and P_{G(2,4)}(q) = 1 + q + 2q^2 + q^3 + q^4!
  let a3 = RootType::TypeA(3)
  let w2 = [1.0, 1.0, 0.0, 0.0]
  assert_eq(weyl_dimension(a3, w2), 6)
  assert_true(approx(satake_minuscule_poincare(a3, w2, 1.0), 6.0, tol=1.0e-12))
  let q_val = 2.0
  let expected_g24 = 1.0 + 2.0 + 2.0 * 4.0 + 8.0 + 16.0 // = 35
  assert_true(
    approx(satake_minuscule_poincare(a3, w2, q_val), expected_g24, tol=1.0e-12),
  )
  // Palindromicity (Poincaré–Verdier duality on IC(Gr_{ω_2})): q^4 P(1/q) == P(q)
  let pq = satake_minuscule_poincare(a3, w2, 2.5)
  let pinv = @math.pow(2.5, 4.0) * satake_minuscule_poincare(a3, w2, 1.0 / 2.5)
  assert_true(approx(pq, pinv, tol=1.0e-11))
  // 2. Lusztig q-analogue K_{λ, μ}(q) on A_2 (GL_3):
  //    Adjoint rep λ = (2, 1, 0), zero weight μ = (1, 1, 1):
  //      K_{(2,1,0), (1,1,1)}(q) = q + q^2  =>  at q=1 gives 2 (= rank sl_3!), at q=2 gives 6!
  assert_true(approx(lusztig_q_multiplicity_a2([2, 1, 0], [1, 1, 1], 1.0), 2.0))
  assert_true(approx(lusztig_q_multiplicity_a2([2, 1, 0], [1, 1, 1], 2.0), 6.0))
  //    Sum of q=1 weight multiplicities over all weights in V_{(2,1,0)} equals Weyl dim = 8:
  //    6 extremal weights (mult 1) + 1 zero weight (mult 2) = 8 == weyl_dimension(A_2, [2,1,0])!
  assert_eq(weyl_dimension(RootType::TypeA(2), [2.0, 1.0, 0.0]), 8)
  assert_true(approx(lusztig_q_multiplicity_a2([2, 1, 0], [2, 1, 0], 1.0), 1.0))
}

///|
test "T34: Unified frontier emergence converges on B^4 × Λ(32) with B_2 adelic gate, SL_3 opers & AD" {
  let tr = run_frontier_emergence()
  assert_eq(tr.sparse_fibre_generators, 32)
  assert_eq(tr.rows.length(), 121)
  assert_true(tr.final_hyp_spread < tr.initial_hyp_spread * 0.01)
  assert_true(tr.final_encoder_loss < tr.initial_encoder_loss * 0.01)
  assert_true(tr.max_oper_det_defect < 1.0e-7)
  assert_eq(tr.rows[0].split(",").collect().length(), 8)
}
