// ============================================================================
// 27_hexad_emergence.mbt -- frame1-Hexad: 128 agents x 640 steps on the
// LFM2.5-Embedding k-NN graph versus the Round-1 ring.
//
// The 128 agents are the role cards of the six-way synthesis
// (40 K2-Horizon-0.9B, 32 Limite-1B, 24 DeepSeek-R1-Distill-1.5B,
// 16 Qwen3.5-0.8B, 16 LFM2.5-Embedding-350M; scripts/hexad/agents.json).
// Their REAL LFM vectors (module 26) define the coupling graph (module 25).
//
// Theorem T48 (replacement certificate).  On the same initial state and the
// same seed, the LFM k-NN topology
//   (a) is connected with λ_2(L_W) >= 0.05  (ring: 2-2cos(2π/128) = 0.0024),
//   (b) keeps R(T) > 0.999 and belief dispersion < 1e-6,
//   (c) ends with a strictly smaller symbolic-neural concept dispersion
//       than the ring.
// ============================================================================

///|
/// Hexad configuration: Round-1 dynamics, 128 agents, 640 steps.
pub fn hexad_config() -> EmergenceConfig {
  { ..default_config(), n: 128, steps: 640, }
}

///|
/// The LFM k-NN coupling graph of the 128 Hexad agents.
pub fn hexad_lfm_graph(k~ : Int32 = 8, tau~ : Double = 0.25) -> CouplingGraph {
  CouplingGraph::knn(
    decode_hex_embeddings(hexad_lfm_hex),
    hexad_agent_model,
    k~,
    tau~,
  )
}

///|
/// Ring vs LFM k-NN on identical initial conditions.
pub(all) struct HexadComparison {
  ring : EmergenceTrace
  knn : EmergenceTrace
  graph : CouplingGraph
  cross_edges : Int32 // edges joining two different model families
}

///|
pub fn run_hexad_comparison(steps~ : Int32 = 640) -> HexadComparison {
  let cfg = { ..hexad_config(), steps, }
  let graph = hexad_lfm_graph()
  let mut cross = 0
  for i in 0..<graph.n {
    for j in graph.nbr[i] {
      if hexad_agent_model[i] != hexad_agent_model[j] {
        cross += 1
      }
    }
  }
  HexadComparison::{
    ring: run_emergence(cfg),
    knn: run_emergence_on(cfg, graph),
    graph,
    cross_edges: cross / 2,
  }
}

///|
pub fn HexadComparison::summary(self : HexadComparison) -> String {
  let r = self.ring
  let k = self.knn
  let t = r.steps
  "agents: \{r.agents}, steps: \{t}, families: K2 40 | Limite 32 | R1 24 | Qwen3.5 16 | LFM 16\n" +
  "graph: \{self.graph.kind}, edges = \{self.graph.edge_count()}, cross-family edges = \{self.cross_edges}\n" +
  "lambda_2: ring = \{r.lambda2} -> knn = \{k.lambda2}\n" +
  "lambda_max: ring = 4 -> knn = \{self.graph.lambda_max}, heat CFL dt: ring = \{r.cfl} -> knn = \{k.cfl}\n" +
  "R(T): ring = \{r.order[t]} | knn = \{k.order[t]}\n" +
  "belief dispersion(T): ring = \{r.dispersion[t]} | knn = \{k.dispersion[t]}\n" +
  "concept dispersion(0) = \{r.concept[0]} -> (T): ring = \{r.concept[t]} | knn = \{k.concept[t]}\n" +
  "max Super-Jacobi defect: knn = \{max_of(k.jacobi_defect)}"
}

///|
/// Side-by-side CSV: step, R, belief and concept dispersion for both graphs.
pub fn HexadComparison::csv(self : HexadComparison) -> String {
  let r = self.ring
  let k = self.knn
  let mut s = "step,order_R_ring,order_R_knn,belief_disp_ring,belief_disp_knn,concept_disp_ring,concept_disp_knn,kl_ring,kl_knn\n"
  for t in 0..<=r.steps {
    s = s +
      "\{t},\{r.order[t]},\{k.order[t]},\{r.dispersion[t]},\{k.dispersion[t]},\{r.concept[t]},\{k.concept[t]},\{r.kl[t]},\{k.kl[t]}\n"
  }
  s
}

///|
pub fn max_of(xs : Array[Double]) -> Double {
  let mut m = 0.0
  for x in xs {
    if x > m {
      m = x
    }
  }
  m
}

///|
test "T48: LFM k-NN graph replaces the ring for 128 agents x 640 steps" {
  let h = run_hexad_comparison()
  let t = h.ring.steps
  assert_true(t == 640)
  assert_true(h.ring.agents == 128)
  // (a) connectivity
  assert_true(h.knn.lambda2 >= 0.05)
  assert_true(h.knn.lambda2 > 20.0 * h.ring.lambda2)
  assert_true(h.cross_edges > 0)
  // (b) the emergent order survives the replacement
  assert_true(h.knn.order[t] > 0.999)
  assert_true(h.knn.dispersion[t] < 1.0e-6)
  assert_true(max_of(h.knn.jacobi_defect) < 1.0e-12)
  // (c) faster symbolic-neural consensus than the ring
  assert_true(h.knn.concept[t] < h.ring.concept[t])
  let rows = h.csv().split("\n").collect()
  assert_true(rows[0].split(",").collect().length() == 9)
}

///|
test "T48: decoded LFM vectors are 128 x 1024 with five model families" {
  let emb = decode_hex_embeddings(hexad_lfm_hex)
  assert_true(emb.length() == 128)
  for v in emb {
    assert_true(v.length() == 1024)
  }
  let counts = Array::make(5, 0)
  for m in hexad_agent_model {
    counts[m] = counts[m] + 1
  }
  assert_true(counts == [40, 32, 24, 16, 16])
}
