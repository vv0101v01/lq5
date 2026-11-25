// ============================================================================
// 16_langlands_emergence.mbt -- the coupled emergence pipeline.
//
// The four new mathematical objects are run as *one* dynamical system on a
// population of agents.  Each agent carries
//
//   * a point of the Betti moduli space  (x,y,z) ∈ M_B(Σ_{1,1}, SL_2),
//     i.e. the trace coordinates of a `π_1`-representation (what the agent
//     "represents" about its world),
//   * a non-semisimple L-parameter encoded by its unipotent defect, i.e. a
//     Jordan type; the defect measures how far the parameter is from
//     semisimple ("how locked" the agent is),
//   * an effective spectral parameter λ, from which the c-function gate and
//     the ^G-operator grading are read off.
//
// One step performs, in this order:
//
//   (1) c-function gate      g_i = |c(w0, λ_i)|, the intertwining amplitude of
//                            the agent's Eisenstein series (module 12);
//   (2) unipotent cooling    defect ← defect·(1 - decay·(1 - 1/(1+g)))
//                            (modules 13 and 14: the Jordan type relaxes
//                            towards the semisimple one);
//   (3) Hamiltonian alignment on each symplectic leaf of the character
//                            variety:  p ← p + dt·{p, H}  with the Poisson
//                            bracket of module 15.  This *preserves* the
//                            Casimir κ, so the measured drift is a certificate
//                            for Theorems T23b;
//   (4) Eisenstein projection  p ← p + η·(p̄ - p):  the rank-one spectral
//                            projection of Theorem T20 is what allows agents
//                            to leave their leaf and agree (consensus).
//
// Recorded per step: order parameter R, κ-spread, trace-spread, Jacobi
// defect, mean unipotent defect, Eisenstein fraction, Casimir drift and the
// ^G-operator fixed-space ratio.
// ============================================================================

///|
/// Configuration of the coupled run.
pub(all) struct LanglandsConfig {
  agents : Int32
  steps : Int32
  genus : Int32
  punctures : Int32
  eta : Double // Eisenstein projection (consensus) strength
  coupling : Double // Hamiltonian flow strength
  decay : Double // unipotent cooling rate
  dt : Double // flow step
  seed : Int32
}

///|
/// Default configuration (documented in README.md).
pub fn langlands_default_config() -> LanglandsConfig {
  {
    agents: 24,
    steps: 200,
    genus: 1,
    punctures: 1,
    eta: 0.55,
    coupling: 1.2,
    decay: 0.12,
    dt: 0.05,
    seed: 20261117,
  }
}

///|
/// Result of a run: CSV rows plus the final certificate bundle.
pub(all) struct LanglandsTrace {
  rows : Array[String]
  summary : String
  final_order : Double
  final_kappa_spread : Double
  final_defect : Double
  max_casimir_drift : Double
  max_jacobi : Double
}

///|
/// CSV header of the trace.
pub fn langlands_csv_header() -> String {
  "step,order_R,kappa_spread,trace_spread,jacobi_defect,unipotent_defect,eisenstein_fraction,casimir_drift,gcheck_fixed_ratio"
}

///|
/// Spread (standard deviation) of an array.
pub fn spread(x : Array[Double]) -> Double {
  let n = x.length()
  if n == 0 {
    return 0.0
  }
  let m = mean_of(x)
  let mut s = 0.0
  for v in x {
    s += (v - m) * (v - m)
  }
  (s / n.to_double()).sqrt()
}

///|
/// Arithmetic mean of an array.
pub fn mean_of(x : Array[Double]) -> Double {
  let mut s = 0.0
  for v in x {
    s += v
  }
  s / x.length().to_double()
}

///|
/// The full coupled run.  Deterministic given the seed.
pub fn run_langlands(cfg : LanglandsConfig) -> LanglandsTrace {
  let rng = Rng::new(cfg.seed)
  let n = cfg.agents
  let xs : Array[Double] = []
  let ys : Array[Double] = []
  let zs : Array[Double] = []
  let defects : Array[Double] = []
  for i in 0..<n {
    // A random point of the one-holed-torus character variety: the trace
    // coordinates of a random SL_2-representation.
    let a = random_sl2(rng)
    let b = random_sl2(rng)
    xs.push(mat2_tr(a))
    ys.push(mat2_tr(b))
    zs.push(mat2_tr(mat2_mul(a, b)))
    // Non-semisimple start: Jordan types of GL_3, with the defect of the type.
    let part : Array[Int32] = if i % 3 == 0 {
      [3]
    } else if i % 3 == 1 {
      [2, 1]
    } else {
      [1, 1, 1]
    }
    defects.push((unipotent_defect(part) + 1).to_double() / 3.0)
  }
  let rows : Array[String] = []
  // initial width of the order parameter
  let k0 : Array[Double] = []
  for i in 0..<n {
    k0.push(fricke_kappa(xs[i], ys[i], zs[i]))
  }
  let kappa0 = spread(k0)
  let mut worst_casimir_drift = 0.0
  let mut worst_jacobi = 0.0
  let mut first_def = 0.0
  let mut first_order = 0.0
  let mut first_ks = 0.0
  let mut last_order = 0.0
  let mut last_ks = 0.0
  let mut last_def = 0.0
  for step in 0..<=cfg.steps {
    // ---- measurements -----------------------------------------------------
    let kappas : Array[Double] = []
    for i in 0..<n {
      kappas.push(fricke_kappa(xs[i], ys[i], zs[i]))
    }
    let ks = spread(kappas)
    let txs : Array[Double] = []
    let tzs : Array[Double] = []
    for i in 0..<n {
      txs.push(xs[i])
      tzs.push(zs[i])
    }
    let trace_sp = spread(txs) + spread(tzs)
    let order = 1.0 - ks / (kappa0 + 1.0e-12)
    let mp : Array[Double] = [mean_of(xs), mean_of(ys), mean_of(zs)]
    let jd = jacobi_defect(
      fn(q : Array[Double]) { q[0] },
      fn(q : Array[Double]) { q[1] * q[1] },
      fn(q : Array[Double]) { fricke_kappa(q[0], q[1], q[2]) },
      mp,
    )
    let mut def_sum = 0.0
    for d in defects {
      def_sum += d
    }
    let def_mean = def_sum / n.to_double()
    let mut eis_cnt = 0.0
    for d in defects {
      if d < 0.05 {
        eis_cnt += 1.0
      }
    }
    let eis = eis_cnt / n.to_double()
    // ^G-operator "grading activation": def² · Σ_{ij}(w_i - w_j)² / n² for the
    // reference weights w = (1,-1,0) of type (2,1); it vanishes exactly when
    // the parameters become semisimple (h = 0).
    let fixed_ratio = def_mean * def_mean * 12.0 / 9.0
    if jd > worst_jacobi {
      worst_jacobi = jd
    }
    if step == 0 {
      first_order = order
      first_ks = ks
      first_def = def_mean
    }
    last_order = order
    last_ks = ks
    last_def = def_mean
    let row = step.to_string() +
      "," +
      order.to_string() +
      "," +
      ks.to_string() +
      "," +
      trace_sp.to_string() +
      "," +
      jd.to_string() +
      "," +
      def_mean.to_string() +
      "," +
      eis.to_string() +
      "," +
      worst_casimir_drift.to_string() +
      "," +
      fixed_ratio.to_string()
    rows.push(row)
    if step == cfg.steps {
      break
    }
    // ---- (1)+(2) c-function gate and unipotent cooling --------------------
    for i in 0..<n {
      let lambda : Array[Double] = [
        defects[i] + 0.5,
        1.0 - defects[i],
        zs[i] * 0.25,
      ]
      let gate = c_function(weyl_longest(3), lambda).abs()
      let factor = 1.0 - cfg.decay * (1.0 - 1.0 / (1.0 + gate))
      let nd = defects[i] * factor
      defects[i] = if nd < 1.0e-4 { 0.0 } else { nd }
    }
    // ---- (3) Hamiltonian alignment inside each symplectic leaf ------------
    let mx = mean_of(xs)
    let my = mean_of(ys)
    let mz = mean_of(zs)
    for i in 0..<n {
      let x0 = xs[i]
      let y0 = ys[i]
      let z0 = zs[i]
      let half = -0.5 * cfg.coupling
      let ham = fn(q : Array[Double]) {
        let dx = q[0] - mx
        let dy = q[1] - my
        let dz = q[2] - mz
        half * (dx * dx + dy * dy + dz * dz)
      }
      let p0 : Array[Double] = [x0, y0, z0]
      let px = poisson_bracket(fn(q : Array[Double]) { q[0] }, ham, p0)
      let py = poisson_bracket(fn(q : Array[Double]) { q[1] }, ham, p0)
      let pz = poisson_bracket(fn(q : Array[Double]) { q[2] }, ham, p0)
      // Structural certificate: dκ/dt = {κ, H} must vanish (Theorem T23b).
      let kap = fn(q : Array[Double]) { fricke_kappa(q[0], q[1], q[2]) }
      let casimir_bracket = poisson_bracket(kap, ham, p0).abs()
      if casimir_bracket > worst_casimir_drift {
        worst_casimir_drift = casimir_bracket
      }
      // Euler step followed by exact leaf projection (κ = const).
      let tentative : Array[Double] = [
        x0 + cfg.dt * px,
        y0 + cfg.dt * py,
        z0 + cfg.dt * pz,
      ]
      let proj = project_to_leaf(p0, tentative)
      xs[i] = proj[0]
      ys[i] = proj[1]
      zs[i] = proj[2]
    }
    // ---- (4) Eisenstein projection (consensus, leaves can be crossed) -----
    let mx2 = mean_of(xs)
    let my2 = mean_of(ys)
    let mz2 = mean_of(zs)
    for i in 0..<n {
      xs[i] = xs[i] + cfg.eta * (mx2 - xs[i])
      ys[i] = ys[i] + cfg.eta * (my2 - ys[i])
      zs[i] = zs[i] + cfg.eta * (mz2 - zs[i])
    }
  }
  let first = rows[0]
  let last_row = rows[rows.length() - 1]
  let summary = "# langlands emergence / Eisenstein + Betti + ^G / MoonBit\n" +
    "agents: \{cfg.agents}, steps: \{cfg.steps}, surface: genus \{cfg.genus} with \{cfg.punctures} puncture(s)\n" +
    "dim M_B(SL_2) = \{betti_dim(cfg.genus, cfg.punctures, 2)} (Poisson rank 2, one Casimir kappa)\n" +
    "first row  : \{first}\n" +
    "last row   : \{last_row}\n" +
    "R(0) = \{first_order} -> R(T) = \{last_order}, kappa spread \{first_ks} -> \{last_ks}\n" +
    "unipotent defect \{first_def} -> \{last_def}\n" +
    "max |{kappa, H}| (Casimir certificate) = \{worst_casimir_drift}\n" +
    "max Jacobi defect over the whole run = \{worst_jacobi}\n" +
    "Eisenstein projection defect (rank 2) = \{eisenstein_datum([1.4, -0.9]).projection_defect()}\n" +
    "Eisenstein projection defect (rank 3) = \{eisenstein_datum([1.7, -0.4, 2.3]).projection_defect()}\n" +
    "c-function functional-equation defect = \{c_functional_eq_defect([1.7, -0.4, 2.3])}\n" +
    "CT-Eis defect (GL_4, type (2,2)) = \{ct_eis_defect([2, 2], [1.5, 1.5, -0.5, -0.5])}\n" +
    "unipotent orbit dim, type (2,2) = \{ct_fiber_dim([2, 2])}, ^G fixed algebra = \{gcheck_kernel_dim([2, 2])}, total = 16\n" +
    "dim ker Theta (3,1) = \{gcheck_kernel_dim([3, 1])}, dim Cent(N) = \{centralizer_dim_formula([3, 1])}, Ext^1 = \{ext1_dim([3, 1])}\n" +
    "nilpotent orbits in GL_6 = \{nilpotent_orbit_count(6)}\n"
  {
    rows,
    summary,
    final_order: last_order,
    final_kappa_spread: last_ks,
    final_defect: last_def,
    max_casimir_drift: worst_casimir_drift,
    max_jacobi: worst_jacobi,
  }
}

///|
/// CSV text of a trace (header + rows).
pub fn LanglandsTrace::csv(self : LanglandsTrace) -> String {
  let mut out = langlands_csv_header() + "\n"
  for r in self.rows {
    out = out + r + "\n"
  }
  out
}

// ---------------------------------------------------------------------------
// Tests: the coupled pipeline must reproduce the certificates.
// ---------------------------------------------------------------------------

///|
test "T24: coupled run converges and preserves every certificate" {
  let cfg = langlands_default_config()
  let tr = run_langlands(cfg)
  assert_true(tr.rows.length() == cfg.steps.to_int() + 1)
  // order parameter goes from ~0 to ~1 and the kappa spread collapses
  assert_true(tr.final_order > 0.99)
  assert_true(tr.final_kappa_spread < 1.0e-3)
  // unipotent defects cool down to (almost) zero: semisimplification
  assert_true(tr.final_defect < 1.0e-3)
  // the Hamiltonian step never leaks the Casimir (Theorem T23b certificate)
  assert_true(tr.max_casimir_drift < 1.0e-5)
  // the Jacobi identity survives the whole run
  assert_true(tr.max_jacobi < 1.0e-4)
  // the trace carries the documented header
  assert_true(tr.rows[0].split(",").collect().length() == 9)
}
