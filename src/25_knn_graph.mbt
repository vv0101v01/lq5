// ============================================================================
// 25_knn_graph.mbt -- Coupling graphs for the agent population (frame1-Hexad).
//
// ABSTRACTION.  Round 1 hard-wired the agent graph to the N-cycle: every
// consensus channel (Kuramoto, Fisher beliefs, symbolic-neural unification)
// only talked to i-1 and i+1, so the algebraic connectivity is
//   λ_2(ring) = 2 - 2 cos(2π/N)  ~  4π²/N²,
// i.e. 0.0024 at N = 128.  A `CouplingGraph` abstracts the topology:
//
//   * `CouplingGraph::ring_pope` reproduces Round 1 exactly (bit-identical
//     traces: neighbour weights 1/2 = exact binary scaling, PoPE Kuramoto
//     weights copied verbatim);
//   * `CouplingGraph::knn` builds the graph from REAL embedding vectors
//     (LFM2.5-Embedding-350M, module 26): symmetric k-NN with heat-kernel
//     weights w_ij = exp((cos_ij - 1)/τ), plus one bridge from every agent to
//     its most similar agent of every other model family, plus a Kruskal
//     repair pass that guarantees connectivity.
//
// Theorem T47 (graph replacement is conservative).  For any connected
// CouplingGraph, L = D - W is symmetric PSD with L·1 = 0 and λ_2 > 0; the
// row-stochastic neighbour mean P = D^{-1} W preserves the simplex; and for
// the ring P_ij = 1/2 so `neighbour_mean` equals Round 1's `(l + r) * 0.5`
// to the last bit.
// Algorithm A38: cyclic Jacobi eigen-solver for the full Laplacian spectrum.
// Algorithm A39: k-NN + cross-family bridge + Kruskal repair graph builder.
// ============================================================================

///|
/// Agent topology shared by every consensus channel.
pub(all) struct CouplingGraph {
  n : Int32
  kind : String
  nbr : Array[Array[Int32]] // neighbour lists (ascending, ring: [left, right])
  wt : Array[Array[Double]] // matching edge weights
  kuramoto : Array[Double] // n*n multiplier: θ̇_i += (K/n)·kuramoto[i,j]·sin(θ_j-θ_i)
  laplacian : Array[Double] // dense L = D - W
  lambda2 : Double // algebraic connectivity
  lambda_max : Double // spectral radius of L
}

///|
/// Round-1 topology: the N-cycle for beliefs / concepts, PoPE all-to-all
/// translation-equivariant weights for Kuramoto.
pub fn CouplingGraph::ring_pope(n : Int32, base : Double) -> CouplingGraph {
  let weights = pope_coupling_weights(PolarPoPE::new(8, base), n)
  let kuramoto = Array::makei(n * n, fn(t) {
    let i = t / n
    let j = t % n
    let d = (j - i + n) % n
    let off = if d <= n - d { d } else { n - d }
    weights[off]
  })
  CouplingGraph::{
    n,
    kind: "ring",
    nbr: Array::makei(n, fn(i) { [(i - 1 + n) % n, (i + 1) % n] }),
    wt: Array::makei(n, fn(_) { [1.0, 1.0] }),
    kuramoto,
    laplacian: ring_laplacian(n),
    lambda2: algebraic_connectivity(n),
    lambda_max: ring_lambda_max(n),
  }
}

///|
/// CFL bound dt <= 2/(κ λ_max(L)) for the explicit heat integrator.
pub fn CouplingGraph::heat_cfl(self : CouplingGraph, kappa : Double) -> Double {
  2.0 / (kappa * self.lambda_max)
}

///|
/// Weighted degree of agent `i`.
pub fn CouplingGraph::degree(self : CouplingGraph, i : Int32) -> Double {
  let mut d = 0.0
  for w in self.wt[i] {
    d += w
  }
  d
}

///|
/// Number of undirected edges.
pub fn CouplingGraph::edge_count(self : CouplingGraph) -> Int32 {
  let mut e = 0
  for i in 0..<self.n {
    e += self.nbr[i].length()
  }
  e / 2
}

///|
/// Row-stochastic neighbour mean `Σ_j (w_ij / d_i) v_j` (P = D^{-1} W).
pub fn CouplingGraph::neighbour_mean(
  self : CouplingGraph,
  i : Int32,
  vecs : Array[Array[Double]],
) -> Array[Double] {
  let d = self.degree(i)
  let k = vecs[i].length()
  let acc = Array::make(k, 0.0)
  for e in 0..<self.nbr[i].length() {
    let j = self.nbr[i][e]
    let c = self.wt[i][e] / d
    for t in 0..<k {
      acc[t] = acc[t] + c * vecs[j][t]
    }
  }
  acc
}

///|
/// Kuramoto vector field on an arbitrary coupling matrix.
pub fn kuramoto_derivative_graph(
  theta : Array[Double],
  omega : Array[Double],
  g : CouplingGraph,
  coupling : Double,
) -> Array[Double] {
  let n = theta.length()
  Array::makei(n, fn(i) {
    let mut s = omega[i]
    for j in 0..<n {
      s += coupling /
        n.to_double() *
        g.kuramoto[i * n + j] *
        @math.sin(theta[j] - theta[i])
    }
    s
  })
}

///|
/// Heun (RK2) step on the coupling graph.
pub fn kuramoto_step_graph(
  theta : Array[Double],
  omega : Array[Double],
  g : CouplingGraph,
  coupling : Double,
  dt : Double,
) -> Array[Double] {
  let k1 = kuramoto_derivative_graph(theta, omega, g, coupling)
  let mid = Array::makei(theta.length(), fn(i) { theta[i] + dt * k1[i] })
  let k2 = kuramoto_derivative_graph(mid, omega, g, coupling)
  Array::makei(theta.length(), fn(i) { theta[i] + 0.5 * dt * (k1[i] + k2[i]) })
}

///|
/// Symbolic-neural consensus on the graph (generalises
/// `consensus_unification`, whose neighbourhood is the ring).
pub fn consensus_unification_graph(
  g : CouplingGraph,
  agents : Array[Hybrid],
  rounds : Int32,
  eta~ : Double = 0.5,
) -> Array[Hybrid] {
  let n = agents.length()
  let mut cur = agents.copy()
  for _ in 0..<rounds {
    let vecs = Array::makei(n, fn(i) { cur[i].vec })
    let next = Array::makei(n, fn(i) {
      let lead = cur[g.nbr[i][0]]
      let nb = Hybrid::{
        term: lead.term,
        vec: g.neighbour_mean(i, vecs),
        subst: lead.subst,
      }
      unification_operator(cur[i], nb, eta)
    })
    cur = next
  }
  cur
}

// ---------------------------------------------------------------------------
// A38: cyclic Jacobi eigenvalues of a dense symmetric matrix.
// ---------------------------------------------------------------------------

///|
/// All eigenvalues of the symmetric `n x n` matrix `a`, ascending.
pub fn sym_eigenvalues(
  a : Array[Double],
  n : Int32,
  sweeps~ : Int32 = 60,
) -> Array[Double] {
  let m = a.copy()
  for _ in 0..<sweeps {
    let mut off = 0.0
    for p in 0..<n {
      for q in (p + 1)..<n {
        off += m[p * n + q] * m[p * n + q]
      }
    }
    if off < 1.0e-26 {
      break
    }
    for p in 0..<n {
      for q in (p + 1)..<n {
        let apq = m[p * n + q]
        if apq.abs() > 1.0e-300 {
          let theta = (m[q * n + q] - m[p * n + p]) / (2.0 * apq)
          let sgn = if theta >= 0.0 { 1.0 } else { -1.0 }
          let t = sgn / (theta.abs() + (theta * theta + 1.0).sqrt())
          let c = 1.0 / (t * t + 1.0).sqrt()
          let s = t * c
          for k in 0..<n {
            let akp = m[k * n + p]
            let akq = m[k * n + q]
            m[k * n + p] = c * akp - s * akq
            m[k * n + q] = s * akp + c * akq
          }
          for k in 0..<n {
            let apk = m[p * n + k]
            let aqk = m[q * n + k]
            m[p * n + k] = c * apk - s * aqk
            m[q * n + k] = s * apk + c * aqk
          }
        }
      }
    }
  }
  sort_total(Array::makei(n, fn(i) { m[i * n + i] }))
}

// ---------------------------------------------------------------------------
// A39: k-NN graph from embedding vectors.
// ---------------------------------------------------------------------------

///|
/// Decode rows of signed bytes stored as two lowercase hex digits each
/// (module 26 stores the LFM vectors this way).
pub fn decode_hex_embeddings(rows : Array[String]) -> Array[Array[Double]] {
  rows.map(fn(row) {
    let out : Array[Double] = []
    let mut hi = -1
    for ch in row {
      let c = ch.to_int()
      let v = if c >= 97 { c - 87 } else { c - 48 }
      if hi < 0 {
        hi = v
      } else {
        let b = hi * 16 + v
        let sb = if b > 127 { b - 256 } else { b }
        out.push(sb.to_double())
        hi = -1
      }
    }
    out
  })
}

///|
/// Cosine-similarity matrix of the rows of `emb` (n*n, row-major).
pub fn cosine_matrix(emb : Array[Array[Double]]) -> Array[Double] {
  let n = emb.length()
  let unit = emb.map(fn(v) { vec_scale(v, 1.0 / dot(v, v).sqrt()) })
  Array::makei(n * n, fn(t) { dot(unit[t / n], unit[t % n]) })
}

///|
/// Connected-component label of every vertex of the dense weight matrix.
pub fn graph_components(w : Array[Double], n : Int32) -> Array[Int32] {
  let label = Array::make(n, -1)
  let mut c = 0
  for s in 0..<n {
    if label[s] < 0 {
      let stack = [s]
      label[s] = c
      while stack.length() > 0 {
        let u = stack.unsafe_pop()
        for v in 0..<n {
          if w[u * n + v] > 0.0 && label[v] < 0 {
            label[v] = c
            stack.push(v)
          }
        }
      }
      c += 1
    }
  }
  label
}

///|
/// Build the LFM k-NN coupling graph (A39).
/// * `labels[i]` = model family of agent i (cross-family bridges);
/// * `k` nearest neighbours by cosine, heat kernel `exp((cos - 1)/tau)`.
pub fn CouplingGraph::knn(
  emb : Array[Array[Double]],
  labels : Array[Int32],
  k~ : Int32 = 8,
  tau~ : Double = 0.25,
) -> CouplingGraph {
  let n = emb.length()
  let sim = cosine_matrix(emb)
  let w = Array::make(n * n, 0.0)
  fn link(i : Int32, j : Int32) {
    let v = @math.exp((sim[i * n + j] - 1.0) / tau)
    w[i * n + j] = v
    w[j * n + i] = v
  }
  // (1) symmetric k-NN (ties broken by the smaller index: deterministic)
  for i in 0..<n {
    let taken = Array::make(n, false)
    taken[i] = true
    for _ in 0..<k {
      let mut best = -1
      for j in 0..<n {
        if !taken[j] && (best < 0 || sim[i * n + j] > sim[i * n + best]) {
          best = j
        }
      }
      taken[best] = true
      link(i, best)
    }
  }
  // (2) one bridge to the most similar agent of every other model family
  let mut families = 0
  for l in labels {
    if l + 1 > families {
      families = l + 1
    }
  }
  for i in 0..<n {
    for f in 0..<families {
      if f != labels[i] {
        let mut best = -1
        for j in 0..<n {
          if labels[j] == f && (best < 0 || sim[i * n + j] > sim[i * n + best]) {
            best = j
          }
        }
        if best >= 0 {
          link(i, best)
        }
      }
    }
  }
  // (3) Kruskal repair: join components by their most similar pair
  let mut comp = graph_components(w, n)
  let mut done = false
  while !done {
    let mut bi = -1
    let mut bj = -1
    for i in 0..<n {
      for j in 0..<n {
        if comp[i] != comp[j] && (bi < 0 || sim[i * n + j] > sim[bi * n + bj]) {
          bi = i
          bj = j
        }
      }
    }
    if bi < 0 {
      done = true
    } else {
      link(bi, bj)
      comp = graph_components(w, n)
    }
  }
  // adjacency lists, Laplacian, Kuramoto multiplier, spectrum
  let nbr : Array[Array[Int32]] = []
  let wt : Array[Array[Double]] = []
  let deg = Array::make(n, 0.0)
  for i in 0..<n {
    let ni : Array[Int32] = []
    let wi : Array[Double] = []
    for j in 0..<n {
      if w[i * n + j] > 0.0 {
        ni.push(j)
        wi.push(w[i * n + j])
        deg[i] = deg[i] + w[i * n + j]
      }
    }
    nbr.push(ni)
    wt.push(wi)
  }
  let laplacian = Array::makei(n * n, fn(t) {
    if t / n == t % n {
      deg[t / n]
    } else {
      -w[t]
    }
  })
  let kuramoto = Array::makei(n * n, fn(t) { n.to_double() * w[t] / deg[t / n] })
  let ev = sym_eigenvalues(laplacian, n)
  CouplingGraph::{
    n,
    kind: "lfm-knn(k=\{k}, tau=\{tau})",
    nbr,
    wt,
    kuramoto,
    laplacian,
    lambda2: ev[1],
    lambda_max: ev[n - 1],
  }
}

///|
test "T47: ring_pope reproduces the Round-1 ring exactly" {
  let n = 16
  let g = CouplingGraph::ring_pope(n, 0.6)
  assert_true(g.edge_count() == n)
  let rng = Rng::new(5)
  let vs = Array::makei(n, fn(_) {
    Array::makei(4, fn(_) { rng.range(0.0, 1.0) })
  })
  for i in 0..<n {
    let a = g.neighbour_mean(i, vs)
    let b = vec_scale(vec_add(vs[(i - 1 + n) % n], vs[(i + 1) % n]), 0.5)
    for t in 0..<4 {
      assert_true(a[t] == b[t]) // bit-identical, not approximate
    }
  }
  assert_true(g.heat_cfl(0.7) == heat_cfl(0.7, n))
  let theta = Array::makei(n, fn(_) { rng.range(0.0, 6.0) })
  let omega = Array::make(n, 0.01)
  let w = pope_coupling_weights(PolarPoPE::new(8, 0.6), n)
  let x = kuramoto_step_graph(theta, omega, g, 3.2, 0.05)
  let y = kuramoto_step(theta, omega, w, 3.2, 0.05)
  for i in 0..<n {
    assert_true(x[i] == y[i])
  }
}

///|
test "A38: Jacobi spectrum matches the analytic ring spectrum" {
  let n = 16
  let ev = sym_eigenvalues(ring_laplacian(n), n)
  let exact = sort_total(
    Array::makei(n, fn(k) {
      2.0 - 2.0 * @math.cos(2.0 * @math.PI * k.to_double() / n.to_double())
    }),
  )
  for i in 0..<n {
    assert_true(approx(ev[i], exact[i], tol=1.0e-10))
  }
  assert_true(approx(ev[1], algebraic_connectivity(n), tol=1.0e-10))
  assert_true(approx(ev[n - 1], ring_lambda_max(n), tol=1.0e-10))
}

///|
test "A39: k-NN graph is connected, symmetric, PSD with L·1 = 0" {
  // two well-separated clusters: plain k-NN would disconnect them
  let rng = Rng::new(11)
  let emb = Array::makei(12, fn(i) {
    let c = if i < 6 { 1.0 } else { -1.0 }
    [c + 0.05 * rng.range(-1.0, 1.0), 0.05 * rng.range(-1.0, 1.0), 0.3]
  })
  let labels = Array::makei(12, fn(i) { if i < 6 { 0 } else { 1 } })
  let g = CouplingGraph::knn(emb, labels, k=2, tau=0.25)
  let n = g.n
  let comp = graph_components(
    Array::makei(n * n, fn(t) {
      if t / n == t % n {
        0.0
      } else {
        -g.laplacian[t]
      }
    }),
    n,
  )
  for c in comp {
    assert_true(c == 0)
  }
  for i in 0..<n {
    let mut row = 0.0
    for j in 0..<n {
      row += g.laplacian[i * n + j]
      assert_true(g.laplacian[i * n + j] == g.laplacian[j * n + i])
    }
    assert_true(row.abs() < 1.0e-12)
    let p = g.neighbour_mean(i, Array::makei(n, fn(_) { [1.0] }))
    assert_true(approx(p[0], 1.0, tol=1.0e-12)) // row-stochastic
  }
  assert_true(g.lambda2 > 0.0)
  let ev = sym_eigenvalues(g.laplacian, n)
  assert_true(ev[0].abs() < 1.0e-10)
  assert_true(
    approx(
      largest_eigenvalue_power(g.laplacian, n, 2000),
      g.lambda_max,
      tol=1.0e-6,
    ),
  )
  // without the cross-family bridges + repair, plain 2-NN is disconnected
  let sim = cosine_matrix(emb)
  assert_true(sim[0 * 12 + 11] < 0.0)
}
