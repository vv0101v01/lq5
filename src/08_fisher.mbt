// ============================================================================
// 08_fisher.mbt -- Fisher information metric and natural gradient.
//
// ABSTRACTION.  Each agent holds a categorical belief p ∈ Δ^{K-1}.  The
// Fisher-Rao metric on the simplex is g_p(u,v) = Σ_i u_i v_i / p_i and the
// associated COMETRIC (Fisher information matrix on score functions) is
//   F(p) = diag(p) - p pᵀ.
// Theorem T13:
//   (a) F is PSD of rank K-1 and is the inverse metric restricted to the
//       tangent space; KL(p+εu || p) = ½ ε² g_p(u,u) + O(ε³).
//   (b) The natural gradient of mean fitness L(p) = ⟨p,f⟩ is
//       F f = p ⊙ (f - ⟨p,f⟩), i.e. exactly the REPLICATOR equation.
//   (c) Natural gradient is equivariant under any reparameterisation
//       (here: a permutation of categories).
// Algorithm A8: damped CG solves (F + λI) x = ∇L in O(K · iters) time without
// materialising a dense inverse (F is singular on the simplex).
// ============================================================================

///|
/// Fisher information matrix F(p) = diag(p) - p pᵀ.
pub fn fisher_matrix(p : Array[Double]) -> Array[Double] {
  let k = p.length()
  Array::makei(k * k, fn(t) {
    let i = t / k
    let j = t % k
    (if i == j { p[i] } else { 0.0 }) - p[i] * p[j]
  })
}

///|
/// Matrix-free cometric action: F v = p ⊙ (v - ⟨p,v⟩).
pub fn fisher_apply(p : Array[Double], v : Array[Double]) -> Array[Double] {
  let pv = dot(p, v)
  Array::makei(p.length(), fn(i) { p[i] * (v[i] - pv) })
}

///|
pub fn fisher_damped_apply(
  p : Array[Double],
  v : Array[Double],
  lambda : Double,
) -> Array[Double] {
  let fv = fisher_apply(p, v)
  Array::makei(p.length(), fn(i) { fv[i] + lambda * v[i] })
}

///|
/// Fisher-Rao quadratic form g_p(u,u) = Σ u_i²/p_i.
pub fn fisher_rao_quadratic(p : Array[Double], u : Array[Double]) -> Double {
  let mut s = 0.0
  for i in 0..<p.length() {
    s += u[i] * u[i] / (p[i] + 1.0e-15)
  }
  s
}

///|
pub fn kl_divergence(p : Array[Double], q : Array[Double]) -> Double {
  let mut s = 0.0
  for i in 0..<p.length() {
    if p[i] > 0.0 {
      s += p[i] * @math.ln(p[i] / (q[i] + 1.0e-15))
    }
  }
  s
}

///|
pub fn entropy(p : Array[Double]) -> Double {
  let mut s = 0.0
  for v in p {
    if v > 0.0 {
      s -= v * @math.ln(v)
    }
  }
  s
}

///|
/// Damped natural gradient x = (F + λI)^{-1} grad by conjugate gradients.
pub fn natural_gradient_cg(
  p : Array[Double],
  grad : Array[Double],
  lambda : Double,
  iters : Int32,
) -> Array[Double] {
  let k = p.length()
  let mut x = Array::make(k, 0.0)
  let mut r = grad.copy()
  let mut dir = r.copy()
  let mut rs = dot(r, r)
  for _ in 0..<iters {
    let ad = fisher_damped_apply(p, dir, lambda)
    let denom = dot(dir, ad)
    if denom.abs() < 1.0e-18 {
      return x
    }
    let alpha = rs / denom
    let x_new = Array::makei(k, fn(i) { x[i] + alpha * dir[i] })
    let r_new = Array::makei(k, fn(i) { r[i] - alpha * ad[i] })
    let rs_new = dot(r_new, r_new)
    if rs_new < 1.0e-20 {
      return x_new
    }
    let beta = rs_new / rs
    let dir_new = Array::makei(k, fn(i) { r_new[i] + beta * dir[i] })
    x = x_new
    r = r_new
    dir = dir_new
    rs = rs_new
  }
  x
}

///|
pub fn natural_gradient_residual(
  p : Array[Double],
  grad : Array[Double],
  lambda : Double,
  x : Array[Double],
) -> Double {
  let lhs = fisher_damped_apply(p, x, lambda)
  vec_inf_diff(lhs, grad)
}

///|
/// Exact undamped natural gradient: the replicator step
/// p_i(t+1) ∝ p_i(t) exp(η (f_i - ⟨p,f⟩)).
pub fn replicator_step(
  p : Array[Double],
  fitness : Array[Double],
  eta : Double,
) -> Array[Double] {
  let mean = dot(p, fitness)
  normalize(
    Array::makei(p.length(), fn(i) {
      p[i] * @math.exp(clamp(eta * (fitness[i] - mean), -30.0, 30.0))
    }),
  )
}

///|
/// Reversed category order: a smooth reparameterisation of Δ^{K-1}.
pub fn reverse_categories(x : Array[Double]) -> Array[Double] {
  let k = x.length()
  Array::makei(k, fn(i) { x[k - 1 - i] })
}

///|
test "T13(a): F is PSD, singular, and equals the KL Hessian" {
  let p = normalize([0.1, 0.2, 0.3, 0.4])
  let f = fisher_matrix(p)
  let k = p.length()
  // symmetric
  for i in 0..<k {
    for j in 0..<k {
      assert_true(approx(f[i * k + j], f[j * k + i], tol=1.0e-14))
    }
  }
  assert_true(gram_min_eigenvalue(f, k) > -1.0e-9)
  assert_true(determinant(f, k).abs() < 1.0e-10) // rank K-1
  let v = [0.3, -0.4, 0.2, -0.1]
  assert_true(vec_inf_diff(fisher_apply(p, v), matvec(f, v, k)) < 1.0e-12)
  // KL(p+εu||p) + KL(p-εu||p) = ε² g_p(u,u) at second order
  let u = [0.002, -0.001, 0.0015, -0.0025]
  let pp = Array::makei(k, fn(i) { p[i] + u[i] })
  let pm = Array::makei(k, fn(i) { p[i] - u[i] })
  let hess = kl_divergence(pp, p) + kl_divergence(pm, p)
  assert_true(approx(hess, fisher_rao_quadratic(p, u), tol=1.0e-5))
}

///|
test "T13(b): replicator flow is the undamped natural gradient" {
  let p = normalize([0.15, 0.25, 0.35, 0.25])
  let fitness = [1.0, -0.5, 0.25, 0.75]
  let eta = 0.01
  let pr = replicator_step(p, fitness, eta)
  // linearisation: p + η F f = p + η p⊙(f - ⟨p,f⟩)
  let ng = fisher_apply(p, fitness)
  let lin = normalize(Array::makei(p.length(), fn(i) { p[i] + eta * ng[i] }))
  assert_true(vec_inf_diff(pr, lin) < 5.0e-3)
  assert_true(approx(mass(pr), 1.0, tol=1.0e-12))
}

///|
test "T13(c): damped CG solves (F+λI)x=g and is reparameterisation-equivariant" {
  let p = normalize([0.1, 0.35, 0.2, 0.35])
  let grad = [0.4, -0.2, 0.9, -0.5]
  let lambda = 0.05
  let x = natural_gradient_cg(p, grad, lambda, 40)
  assert_true(natural_gradient_residual(p, grad, lambda, x) < 1.0e-10)
  let pr = reverse_categories(p)
  let gr = reverse_categories(grad)
  let xr = natural_gradient_cg(pr, gr, lambda, 40)
  assert_true(vec_inf_diff(reverse_categories(xr), x) < 1.0e-10)
  // entropy increases under the replicator step toward a uniform fitness
  let grad_entropy = Array::makei(p.length(), fn(i) {
    -@math.ln(p[i] + 1.0e-15)
  })
  let pu = replicator_step(p, grad_entropy, 0.5)
  assert_true(entropy(pu) >= entropy(p) - 1.0e-12)
}
