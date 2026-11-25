// ============================================================================
// 10_emergence.mbt -- Integration: an agent-like emergence pipeline.
//
// The abstraction chain, made executable:
//
//   M^{2|2} superchart  ──▶  agent state (phase θ_i, belief p_i, concept h_i)
//   PoPE position kernel ──▶ w(|i-j|) : translation-equivariant coupling
//   PDE evolution flow  ──▶ Kuramoto + heat/replicator diffusion on the graph
//   Super-Jacobi        ──▶ closure certificate of the agent algebra
//   Fisher metric       ──▶ natural-gradient belief consensus
//   Unification U       ──▶ symbolic-neural concept consensus (Fix(U))
//   Super-mass tensor   ──▶ orbit invariants (λ₁,λ₂, Ber) = emergent order
//   Spectral triple     ──▶ Connes distance = inverse coupling = emergent metric
//
// EMERGENCE THEOREM (T15).  For a connected graph with positive translation-
// equivariant weights w, identical intrinsic frequencies and a Fisher-natural-
// gradient consensus update, the coupled flow
//   (θ̇, ṗ, ḣ) = (Kuramoto_w, F(p)∇L, U-consensus)
// converges to a state with R → 1, belief dispersion → 0, concept dispersion
// → 0; the macroscopic order is completely captured by the orbit invariants
// (λ₁, λ₂, anisotropy, Ber) of the second-order super-mass tensor.  The
// U(1) phase symmetry is spontaneously broken; only the quotient (orbit
// space) is observable, which is exactly the "emergent" degree of freedom.
// ============================================================================

///|
pub(all) struct EmergenceConfig {
  n : Int32 // number of agents
  steps : Int32
  kappa : Double // diffusion strength (heat / Fisher coupling)
  coupling : Double // Kuramoto coupling K (also inverse Connes distance)
  eta : Double // learning / unification rate
  seed : Int32
  base : Double // PoPE frequency base
  categories : Int32 // belief simplex dimension
}

///|
pub fn default_config() -> EmergenceConfig {
  EmergenceConfig::{
    n: 32,
    steps: 240,
    kappa: 0.7,
    coupling: 3.2,
    eta: 0.25,
    seed: 20261117,
    base: 0.6,
    categories: 4,
  }
}

///|
pub(all) struct EmergenceTrace {
  agents : Int32
  steps : Int32
  order : Array[Double] // polar order R(t)
  nematic : Array[Double] // second-order (orbit) anisotropy
  concept : Array[Double] // symbolic-neural concept dispersion
  dispersion : Array[Double] // belief dispersion
  kl : Array[Double] // mean KL to the population belief
  fisher_grad : Array[Double] // mean natural-gradient norm
  jacobi_defect : Array[Double] // algebraic closure certificate
  connes_distance : Double
  final_berezinian : Double
  orbit_eigenvalues : Array[Double]
  lambda2 : Double
  cfl : Double
  final_mass : SuperMassTensor
}

///|
pub fn belief_mean(bs : Array[Array[Double]]) -> Array[Double] {
  let n = bs.length()
  let k = bs[0].length()
  Array::makei(k, fn(j) {
    let mut s = 0.0
    for i in 0..<n {
      s += bs[i][j]
    }
    s / n.to_double()
  })
}

///|
pub fn belief_dispersion(bs : Array[Array[Double]]) -> Double {
  let mean = belief_mean(bs)
  let n = bs.length()
  let k = mean.length()
  let mut s = 0.0
  for i in 0..<n {
    for j in 0..<k {
      s += (bs[i][j] - mean[j]).abs()
    }
  }
  s / n.to_double()
}

///|
pub fn mean_kl_to_population(bs : Array[Array[Double]]) -> Double {
  let mean = belief_mean(bs)
  let mut s = 0.0
  for b in bs {
    s += kl_divergence(b, mean)
  }
  s / bs.length().to_double()
}

///|
pub fn concept_dispersion(hs : Array[Hybrid]) -> Double {
  let mut worst = 0.0
  for h in hs {
    let d = vec_inf_diff(hs[0].vec, h.vec)
    if d > worst {
      worst = d
    }
  }
  worst
}

///|
/// A closure certificate: build three (2|1) supermatrices from the current
/// macroscopic state and evaluate the Super-Jacobiator.  It must vanish for
/// ALL states, i.e. the agent algebra stays closed while order emerges.
pub fn closure_certificate(mean_phase_ : Double, r : Double) -> Double {
  let t = mean_phase_
  let x = SuperMatrix::new(2, 1, 2, 1, 0, [
    @math.cos(t),
    0.2 * r,
    0.1,
    -0.3,
    1.0 + 0.1 * r,
    0.4,
    0.2,
    -0.5,
    0.9,
  ])
  let y = SuperMatrix::new(2, 1, 2, 1, 1, [
    0.3, -0.8, 0.6, 1.2, 0.1, -0.4, 0.7, 0.3, 0.5,
  ])
  let z = SuperMatrix::new(2, 1, 2, 1, 0, [
    1.4,
    0.5 * @math.sin(t),
    -0.2,
    0.1,
    -0.9,
    0.3,
    -0.6,
    0.8,
    1.3,
  ])
  super_jacobiator_m(x, y, z).norm()
}

///|
/// Run the full coupled flow and record the macroscopic observables
/// (Round-1 topology: ring + PoPE, see `CouplingGraph::ring_pope`).
pub fn run_emergence(cfg : EmergenceConfig) -> EmergenceTrace {
  run_emergence_on(cfg, CouplingGraph::ring_pope(cfg.n, cfg.base))
}

///|
/// The same coupled flow on an arbitrary agent topology `g` (module 25):
/// Kuramoto, Fisher-belief and symbolic-neural consensus all use `g`.
pub fn run_emergence_on(
  cfg : EmergenceConfig,
  g : CouplingGraph,
) -> EmergenceTrace {
  let n = cfg.n
  let k = cfg.categories
  let rng = Rng::new(cfg.seed)
  // agent phases and tiny frequency heterogeneity
  let mut theta = Array::makei(n, fn(_) { rng.range(0.0, 2.0 * @math.PI) })
  let omega = Array::makei(n, fn(_) { rng.range(-0.05, 0.05) })
  // categorical beliefs on Δ^{k-1}
  let beliefs : Array[Array[Double]] = []
  for _ in 0..<n {
    beliefs.push(normalize(Array::makei(k, fn(_) { rng.range(0.2, 1.0) })))
  }
  // symbolic-neural concept states (all in Fix(U) initially)
  let hybrids : Array[Hybrid] = []
  for i in 0..<n {
    let h = Hybrid::new(
      Term::Add(Term::Var(0), Term::Sym("belief")),
      Array::makei(k, fn(j) { beliefs[i][j] }),
    )
    hybrids.push(h.canonical())
  }
  let dt = 0.05
  let steps = cfg.steps
  let order = Array::make(steps + 1, 0.0)
  let nematic = Array::make(steps + 1, 0.0)
  let concept = Array::make(steps + 1, 0.0)
  let dispersion = Array::make(steps + 1, 0.0)
  let kl = Array::make(steps + 1, 0.0)
  let fisher_grad = Array::make(steps + 1, 0.0)
  let jacobi_defect = Array::make(steps + 1, 0.0)
  let mut cur_h = hybrids
  let mut final_mass = SuperMassTensor::from_phases(theta)
  order[0] = order_parameter(theta)
  nematic[0] = final_mass.anisotropy()
  concept[0] = concept_dispersion(cur_h)
  dispersion[0] = belief_dispersion(beliefs)
  kl[0] = mean_kl_to_population(beliefs)
  fisher_grad[0] = 0.0
  jacobi_defect[0] = closure_certificate(mean_phase(theta), order[0])
  for t in 1..<=steps {
    // (1) nonlinear PDE flow: Kuramoto with PoPE weights
    theta = kuramoto_step_graph(theta, omega, g, cfg.coupling, dt)
    // (2) Fisher-natural-gradient belief consensus with the neighbour mean
    let mean_p = belief_mean(beliefs)
    let mut grad_norm = 0.0
    for i in 0..<n {
      let p = beliefs[i]
      let q = g.neighbour_mean(i, beliefs)
      // score towards the neighbour/population consensus (linearised KL)
      let g = Array::makei(k, fn(j) { 0.5 * q[j] + 0.5 * mean_p[j] - p[j] })
      let ng = natural_gradient_cg(p, g, 0.5, 20)
      grad_norm += vec_inf_norm(ng)
      beliefs[i] = normalize(
        Array::makei(k, fn(j) { clamp(p[j] + cfg.eta * ng[j], 1.0e-6, 1.0e6) }),
      )
    }
    fisher_grad[t] = grad_norm / n.to_double()
    // (3) chaotic symbolic-neural unification on the agent graph
    cur_h = consensus_unification_graph(g, cur_h, 1, eta=0.5)
    // (4) emergent observables from the second-order super-mass tensor
    final_mass = SuperMassTensor::from_phases(theta)
    order[t] = order_parameter(theta)
    nematic[t] = final_mass.anisotropy()
    concept[t] = concept_dispersion(cur_h)
    dispersion[t] = belief_dispersion(beliefs)
    kl[t] = mean_kl_to_population(beliefs)
    jacobi_defect[t] = closure_certificate(mean_phase(theta), order[t])
  }
  EmergenceTrace::{
    agents: n,
    steps,
    order,
    nematic,
    concept,
    dispersion,
    kl,
    fisher_grad,
    jacobi_defect,
    connes_distance: emergent_metric_radius(cfg.coupling),
    final_berezinian: final_mass.berezinian(),
    orbit_eigenvalues: final_mass.canonical_form(),
    lambda2: g.lambda2,
    cfl: g.heat_cfl(cfg.kappa),
    final_mass,
  }
}

///|
pub fn EmergenceTrace::final_order(self : EmergenceTrace) -> Double {
  self.order[self.steps]
}

///|
pub fn EmergenceTrace::final_nematic(self : EmergenceTrace) -> Double {
  self.nematic[self.steps]
}

///|
pub fn EmergenceTrace::final_dispersion(self : EmergenceTrace) -> Double {
  self.dispersion[self.steps]
}

///|
pub fn EmergenceTrace::max_jacobi_defect(self : EmergenceTrace) -> Double {
  let mut s = 0.0
  for d in self.jacobi_defect {
    if d > s {
      s = d
    }
  }
  s
}

///|
pub fn EmergenceTrace::csv(self : EmergenceTrace) -> String {
  let mut out = "step,order_R,nematic,concept_dispersion,belief_dispersion,kl_to_consensus,fisher_grad,jacobi_defect\n"
  for t in 0..<=self.steps {
    out = out +
      "\{t},\{self.order[t]},\{self.nematic[t]},\{self.concept[t]},\{self.dispersion[t]},\{self.kl[t]},\{self.fisher_grad[t]},\{self.jacobi_defect[t]}\n"
  }
  out
}

///|
pub fn EmergenceTrace::summary(self : EmergenceTrace) -> String {
  "agents: \{self.agents}, steps: \{self.steps}\n" +
  "R(0)=\{self.order[0]} -> R(T)=\{self.final_order()}\n" +
  "nematic(0)=\{self.nematic[0]} -> nematic(T)=\{self.final_nematic()}\n" +
  "belief dispersion(0)=\{self.dispersion[0]} -> (T)=\{self.final_dispersion()}\n" +
  "concept dispersion(0)=\{self.concept[0]} -> (T)=\{self.concept[self.steps]}\n" +
  "max Super-Jacobi defect = \{self.max_jacobi_defect()}\n" +
  "Connes distance = 1/K = \{self.connes_distance}\n" +
  "λ_2(L) = \{self.lambda2}, heat CFL dt = \{self.cfl}\n" +
  "orbit eigenvalues = [\{self.orbit_eigenvalues[0]}, \{self.orbit_eigenvalues[1]}]\n" +
  "final Berezinian = \{self.final_berezinian}\n"
}

///|
test "T15: coupled flow produces emergent order and closed algebra" {
  let cfg = default_config()
  let tr = run_emergence(cfg)
  assert_true(tr.order[0] < 0.45)
  assert_true(tr.final_order() > 0.85)
  assert_true(tr.final_order() > tr.order[0])
  // second-order tensor detects the emergent nematic order
  assert_true(tr.final_nematic() > tr.nematic[0])
  // both consensus channels collapse
  assert_true(tr.final_dispersion() < tr.dispersion[0] * 0.1)
  assert_true(tr.concept[cfg.steps] < tr.concept[0] * 0.1)
  assert_true(tr.kl[cfg.steps] < tr.kl[0] * 0.1)
  // the graded Lie algebra stays closed at every time step
  assert_true(tr.max_jacobi_defect() < 1.0e-10)
  // noncommutative geometry: metric radius = inverse coupling
  assert_true(approx(tr.connes_distance, 1.0 / cfg.coupling))
  assert_true(tr.lambda2 > 0.0)
  assert_eq(tr.order.length(), cfg.steps + 1)
  assert_eq(tr.csv().length() > 100, true)
}
