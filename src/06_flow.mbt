// ============================================================================
// 06_flow.mbt -- Partial differential evolution flows on the agent graph.
//
// ABSTRACTION.  A population of N agents is a discrete manifold: the graph
// Laplacian L replaces Δ, and a PDE u_t = κ Δ u becomes the ODE
//   du/dt = -κ L u.
// Theorem T11: Σ u (mass) is conserved, E = ½||u||² decreases monotonically
// because dE/dt = -κ uᵀ L u <= 0, and explicit Euler is stable iff
// dt <= 2/(κ λ_max(L))  (CFL condition, Algorithm A7).
// The nonlinear flow is Kuramoto: θ̇_i = ω_i + (K/N) Σ_j w(|i-j|) sin(θ_j-θ_i),
// i.e. a reaction-diffusion PDE on the graph; w comes from the PoPE kernel and
// depends only on relative position (translation equivariance).
// ============================================================================

///|
/// Discrete Laplacian of the N-cycle: L = 2I - S - Sᵀ.
pub fn ring_laplacian(n : Int32) -> Array[Double] {
  Array::makei(n * n, fn(t) {
    let i = t / n
    let j = t % n
    if i == j {
      2.0
    } else if (i - j + n) % n == 1 || (j - i + n) % n == 1 {
      -1.0
    } else {
      0.0
    }
  })
}

///|
/// λ_k(L) = 2 - 2 cos(2πk/N);  λ_max is the largest.
pub fn ring_lambda_max(n : Int32) -> Double {
  let mut s = 0.0
  for k in 0..<n {
    let v = 2.0 -
      2.0 * @math.cos(2.0 * @math.PI * k.to_double() / n.to_double())
    if v > s {
      s = v
    }
  }
  s
}

///|
/// λ_2(L), the algebraic connectivity: convergence rate to consensus.
pub fn algebraic_connectivity(n : Int32) -> Double {
  2.0 - 2.0 * @math.cos(2.0 * @math.PI / n.to_double())
}

///|
/// CFL bound for the explicit heat integrator: dt <= 2/(κ λ_max).
pub fn heat_cfl(kappa : Double, n : Int32) -> Double {
  2.0 / (kappa * ring_lambda_max(n))
}

///|
/// Apply the (dense) graph Laplacian to `u`; a named alias of `matvec` kept
/// for readability inside the PDE layer.
pub fn laplacian_apply(
  l : Array[Double],
  u : Array[Double],
  n : Int32,
) -> Array[Double] {
  matvec(l, u, n)
}

///|
pub fn heat_step(
  u : Array[Double],
  l : Array[Double],
  kappa : Double,
  dt : Double,
  n : Int32,
) -> Array[Double] {
  let lu = laplacian_apply(l, u, n)
  Array::makei(n, fn(i) { u[i] - kappa * dt * lu[i] })
}

///|
pub fn energy(u : Array[Double]) -> Double {
  0.5 * dot(u, u)
}

///|
pub fn mass(u : Array[Double]) -> Double {
  let mut s = 0.0
  for v in u {
    s += v
  }
  s
}

///|
pub fn centered_energy(u : Array[Double]) -> Double {
  let mean = mass(u) / u.length().to_double()
  let mut s = 0.0
  for v in u {
    s += 0.5 * (v - mean) * (v - mean)
  }
  s
}

///|
pub(all) struct HeatTrace {
  energies : Array[Double]
  masses : Array[Double]
  state : Array[Double]
  dt : Double
}

///|
/// Explicit Euler heat flow with a CFL-safe adaptive step.
pub fn heat_flow(
  u0 : Array[Double],
  l : Array[Double],
  kappa : Double,
  steps : Int32,
  n : Int32,
  safety~ : Double = 0.9,
) -> HeatTrace {
  let dt = safety * heat_cfl(kappa, n)
  let energies = Array::make(steps + 1, 0.0)
  let masses = Array::make(steps + 1, 0.0)
  let mut u = u0.copy()
  energies[0] = energy(u)
  masses[0] = mass(u)
  for t in 1..<=steps {
    u = heat_step(u, l, kappa, dt, n)
    energies[t] = energy(u)
    masses[t] = mass(u)
  }
  HeatTrace::{ energies, masses, state: u, dt, }
}

///|
/// PoPE-derived coupling weights; depends only on the relative offset.
pub fn pope_coupling_weights(k : PolarPoPE, n : Int32) -> Array[Double] {
  Array::makei(n / 2 + 1, fn(off) { 1.0 + 0.5 * k.eval(off) })
}

///|
pub fn kuramoto_derivative(
  theta : Array[Double],
  omega : Array[Double],
  weights : Array[Double],
  coupling : Double,
) -> Array[Double] {
  let n = theta.length()
  Array::makei(n, fn(i) {
    let mut s = omega[i]
    for j in 0..<n {
      let d = (j - i + n) % n
      let off = if d <= n - d { d } else { n - d }
      s += coupling /
        n.to_double() *
        weights[off] *
        @math.sin(theta[j] - theta[i])
    }
    s
  })
}

///|
/// Heun (RK2) step for the nonlinear flow.
pub fn kuramoto_step(
  theta : Array[Double],
  omega : Array[Double],
  weights : Array[Double],
  coupling : Double,
  dt : Double,
) -> Array[Double] {
  let k1 = kuramoto_derivative(theta, omega, weights, coupling)
  let mid = Array::makei(theta.length(), fn(i) { theta[i] + dt * k1[i] })
  let k2 = kuramoto_derivative(mid, omega, weights, coupling)
  Array::makei(theta.length(), fn(i) { theta[i] + 0.5 * dt * (k1[i] + k2[i]) })
}

///|
/// Kuramoto-type order parameter of harmonic `k`:
/// `r_k = |Σ_j e^{i k θ_j}| / n`.  `k = 1` is the usual phase coherence,
/// `k = 2` the nematic (head-tail symmetric) coherence of module 07.
pub fn order_harmonic(theta : Array[Double], harmonic : Int32) -> Double {
  let mut c = 0.0
  let mut s = 0.0
  for x in theta {
    c += @math.cos(harmonic.to_double() * x)
    s += @math.sin(harmonic.to_double() * x)
  }
  let n = theta.length().to_double()
  (c * c + s * s).sqrt() / n
}

///|
/// Phase coherence `r_1` (special case of `order_harmonic`).
pub fn order_parameter(theta : Array[Double]) -> Double {
  order_harmonic(theta, 1)
}

///|
pub fn mean_phase(theta : Array[Double]) -> Double {
  let mut c = 0.0
  let mut s = 0.0
  for x in theta {
    c += @math.cos(x)
    s += @math.sin(x)
  }
  @math.atan2(s, c)
}

///|
pub(all) struct SyncTrace {
  order : Array[Double]
  theta : Array[Double]
}

///|
pub fn kuramoto_flow(
  theta0 : Array[Double],
  omega : Array[Double],
  weights : Array[Double],
  coupling : Double,
  dt : Double,
  steps : Int32,
) -> SyncTrace {
  let order = Array::make(steps + 1, 0.0)
  let mut theta = theta0.copy()
  order[0] = order_parameter(theta)
  for t in 1..<=steps {
    theta = kuramoto_step(theta, omega, weights, coupling, dt)
    order[t] = order_parameter(theta)
  }
  SyncTrace::{ order, theta, }
}

///|
test "T11: heat flow conserves mass, dissipates energy, obeys CFL" {
  let n = 16
  let l = ring_laplacian(n)
  let rng = Rng::new(20261117)
  let u0 = Array::makei(n, fn(_) { rng.range(-1.0, 1.0) })
  let kappa = 0.7
  let tr = heat_flow(u0, l, kappa, 400, n)
  assert_true(tr.dt <= heat_cfl(kappa, n))
  assert_true(approx(mass(tr.state), mass(u0), tol=1.0e-10))
  for t in 0..<tr.masses.length() {
    assert_true(approx(tr.masses[t], mass(u0), tol=1.0e-10))
    if t > 0 {
      assert_true(tr.energies[t] <= tr.energies[t - 1] + 1.0e-12)
    }
  }
  assert_true(tr.energies[400] < tr.energies[0])
  // λ_max(L) verified by power iteration
  let spectral = largest_eigenvalue_power(l, n, 300)
  assert_true(approx(spectral, ring_lambda_max(n), tol=1.0e-8))
  // consensus: the centered energy decays to zero
  assert_true(centered_energy(tr.state) < centered_energy(u0) * 1.0e-3)
}

///|
test "T11: algebraic connectivity controls convergence rate" {
  let n = 8
  let l = ring_laplacian(n)
  let rng = Rng::new(7)
  let u0 = Array::makei(n, fn(_) { rng.range(0.0, 1.0) })
  let cfl = heat_cfl(1.0, n)
  let dt = 0.9 * cfl
  let mut u = u0.copy()
  let e0 = centered_energy(u)
  let mut e1 = e0
  for _ in 0..<100 {
    u = heat_step(u, l, 1.0, dt, n)
    e1 = centered_energy(u)
  }
  // theoretical envelope exp(-κ λ_2 t) with t = 100 dt
  let t = 100.0 * dt
  let envelope = @math.exp(-algebraic_connectivity(n) * t)
  assert_true(e1 < e0)
  assert_true(e1 / e0 <= envelope + 1.0e-6)
}

///|
test "T11: Kuramoto flow breaks U(1) symmetry (emergence)" {
  let n = 24
  let rng = Rng::new(99)
  let theta0 = Array::makei(n, fn(_) { rng.range(0.0, 2.0 * @math.PI) })
  let omega = Array::make(n, 0.0) // identical agents
  let weights = pope_coupling_weights(PolarPoPE::new(8, 0.6), n)
  let tr = kuramoto_flow(theta0, omega, weights, 3.0, 0.05, 600)
  assert_true(tr.order[0] < 0.5)
  assert_true(tr.order[600] > tr.order[0])
  assert_true(tr.order[600] > 0.7)
  // the order parameter is invariant under a global phase shift
  let shifted = Array::makei(n, fn(i) { tr.theta[i] + 1.234 })
  assert_true(
    approx(order_parameter(shifted), order_parameter(tr.theta), tol=1.0e-10),
  )
}
