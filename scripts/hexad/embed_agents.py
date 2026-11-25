#!/usr/bin/env python3
"""Embed the 128 frame1-Hexad agent role cards with the REAL LFM2.5-Embedding-350M.

Requires a llama.cpp server started with the official GGUF:
  llama-server -m LFM2.5-Embedding-350M-Q8_0.gguf --embeddings --port 8090
Outputs:
  scripts/hexad/agents.json              role cards + metadata (model, sha256, dims)
  src/26_lfm_agent_embeddings.mbt        128 x 1024 vectors, int8 per-row scaled, hex
"""
from __future__ import annotations
import argparse, hashlib, json, urllib.request
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]

K2 = ["multi-step task planning and decomposition", "Python programming and debugging",
      "algorithm design and complexity analysis", "data structures and graph algorithms",
      "SQL queries and database schema design", "JavaScript and TypeScript web development",
      "systems programming in C and Rust", "unit testing and test-driven development",
      "API tool calling with JSON arguments", "shell scripting and Linux command line",
      "classical mechanics and kinematics problems", "electromagnetism and circuits",
      "thermodynamics and statistical physics", "general chemistry and stoichiometry",
      "molecular biology and genetics", "numerical methods and floating point error",
      "linear algebra for machine learning", "probability and statistics for data analysis",
      "long document reading and summarisation", "scientific literature question answering",
      "code review and refactoring", "concurrency, locks and race conditions",
      "compiler construction and parsing", "functional programming and type systems",
      "MoonBit language syntax and tooling", "WebAssembly compilation targets",
      "software architecture and design patterns", "version control with git",
      "performance profiling and optimisation", "security vulnerabilities in code",
      "spreadsheet formulas and data cleaning", "unit conversion and dimensional analysis",
      "astronomy and orbital mechanics", "earth science and climate data",
      "economics and quantitative finance", "formal logic and puzzle solving",
      "agentic web search planning", "regular expressions and text processing",
      "numerical optimisation and gradient descent", "operations research and scheduling"]
LIMITE = ["olympiad number theory: divisibility and modular arithmetic",
          "Diophantine equations and integer solutions", "prime numbers and multiplicative functions",
          "Euclidean geometry: triangle centers and circles", "coordinate geometry and conic sections",
          "trigonometric identities and equations", "combinatorics: counting and bijections",
          "combinatorics: pigeonhole principle and extremal arguments",
          "graph theory competition problems", "probability puzzles and expected value",
          "polynomials: roots, Vieta formulas and factorisation", "functional equations",
          "inequalities: AM-GM, Cauchy-Schwarz and Jensen", "sequences, recurrences and generating functions",
          "complex numbers and roots of unity", "calculus: limits, derivatives and integrals",
          "series convergence and summation", "linear algebra: determinants and eigenvalues",
          "abstract algebra: groups and symmetry", "invariants and monovariants in games",
          "3D solid geometry and volumes", "vectors and analytic geometry",
          "logarithms and exponential equations", "floor functions and digit problems",
          "AIME style short-answer problems", "HMMT style team round problems",
          "proof writing by induction", "proof by contradiction and infinite descent",
          "modular inverses and Chinese remainder theorem", "lattice points and Pick's theorem",
          "geometric transformations: rotation, reflection, homothety",
          "Zeta functions and Euler products"]
R1 = ["check every arithmetic step for calculation errors", "verify the logical validity of each inference",
      "search for counterexamples to a claimed statement", "detect unjustified assumptions in a proof",
      "re-derive the final answer by an independent method", "check edge cases and boundary conditions",
      "verify units and orders of magnitude", "check code output by mentally executing it",
      "compare two conflicting solutions and decide which is correct",
      "identify circular reasoning", "check that all problem constraints are satisfied",
      "estimate the answer and test plausibility", "verify algebraic simplifications",
      "check probability calculations for double counting", "audit a chain of thought for skipped steps",
      "verify geometric constructions and diagrams", "detect hallucinated facts or citations",
      "check off-by-one errors in loops and indices", "critique the efficiency of an algorithm",
      "verify a formal definition is applied correctly", "arbitrate between a polarised group of answers",
      "summarise remaining doubts and open gaps", "rate the confidence of a final answer",
      "self-reflect and correct a previous mistake"]
QWEN = ["translate Cantonese user requests into precise English", "translate Traditional Chinese into English",
        "translate the final answer back into Cantonese", "write the final answer in Traditional Chinese",
        "parse an image or screenshot into text", "read a scanned document with OCR",
        "extract the user's intent and constraints", "split a request into subtasks",
        "detect the language and dialect of the input", "ask the user a clarifying question",
        "format the answer as markdown with tables", "call tools and format tool results",
        "summarise a long conversation history", "write friendly conversational replies",
        "explain a result to a non-expert", "render a short certificate summary for the user"]
LFM = ["retrieve relevant documents from long-term memory", "store a verified result into memory",
       "cluster candidate answers by semantic similarity", "route a subtask to the best expert model",
       "deduplicate near-identical candidate solutions", "measure similarity between two agents' outputs",
       "rank documents for a search query", "detect topic drift in a conversation",
       "build a nearest-neighbour graph of agents", "cross-lingual retrieval across languages",
       "find previously solved similar problems", "estimate novelty of a new candidate",
       "index code snippets for retrieval", "index mathematical theorems for retrieval",
       "match a question to a benchmark category", "compress context into a dense vector"]
ROLES = [("K2-Horizon-0.9B", "general reasoning, planning, coding and STEM expert", K2),
         ("Limite-1B-Violetto", "competition mathematics solver", LIMITE),
         ("DeepSeek-R1-Distill-Qwen-1.5B", "long chain-of-thought critic and verifier", R1),
         ("Qwen3.5-0.8B", "multilingual multimodal interface agent", QWEN),
         ("LFM2.5-Embedding-350M", "embedding, memory and routing agent", LFM)]

def embed(url: str, text: str) -> np.ndarray:
    req = urllib.request.Request(url, data=json.dumps({"input": text}).encode(),
                                 headers={"Content-Type": "application/json"})
    v = np.array(json.load(urllib.request.urlopen(req))["data"][0]["embedding"], dtype=np.float64)
    return v / np.linalg.norm(v)

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="http://127.0.0.1:8090/v1/embeddings")
    ap.add_argument("--gguf", default="/tmp/models/lfm-emb-q8.gguf")
    a = ap.parse_args()
    cards = []
    for model, role, specs in ROLES:
        for s in specs:
            cards.append({"id": len(cards), "model": model,
                          "text": f"document: A {role} agent ({model}) specialised in {s}."})
    assert len(cards) == 128, len(cards)
    E = np.stack([embed(a.url, c["text"]) for c in cards])          # 128 x 1024, unit rows
    # int8 per-row symmetric quantisation (cosine is scale-invariant)
    scale = np.abs(E).max(axis=1, keepdims=True) / 127.0
    Q = np.clip(np.rint(E / scale), -127, 127).astype(np.int8)
    Eq = Q.astype(np.float64); Eq /= np.linalg.norm(Eq, axis=1, keepdims=True)
    S, Sq = E @ E.T, Eq @ Eq.T
    k = 8
    def knn(M):
        M = M.copy(); np.fill_diagonal(M, -np.inf); return [set(np.argsort(-r)[:k]) for r in M]
    agree = np.mean([len(x & y) / k for x, y in zip(knn(S), knn(Sq))])
    sha = hashlib.sha256(Path(a.gguf).read_bytes()).hexdigest() if Path(a.gguf).exists() else "n/a"
    meta = {"embedder": "LiquidAI/LFM2.5-Embedding-350M-GGUF (Q8_0)", "gguf_sha256": sha,
            "server": "llama.cpp b12676 llama-server --embeddings", "dim": int(E.shape[1]),
            "agents": len(cards), "int8_max_cosine_error": float(np.abs(S - Sq).max()),
            "int8_knn8_agreement": float(agree)}
    (ROOT / "scripts/hexad/agents.json").write_text(
        json.dumps({"meta": meta, "agents": cards}, ensure_ascii=False, indent=1), encoding="utf-8")
    hexrows = ["".join(f"{(int(x) & 0xFF):02x}" for x in row) for row in Q]
    out = ["// ============================================================================",
           "// 26_lfm_agent_embeddings.mbt -- GENERATED by scripts/hexad/embed_agents.py.",
           "// Real LFM2.5-Embedding-350M (Q8_0 GGUF, llama.cpp b12676) vectors of the",
           "// 128 frame1-Hexad agent role cards (scripts/hexad/agents.json).",
           f"// dim = 1024, int8 per-row quantised; max cosine error vs float = {meta['int8_max_cosine_error']:.2e},",
           f"// 8-NN agreement with float = {agree:.4f}.  GGUF sha256 = {sha}",
           "// Do not edit by hand.",
           "// ============================================================================",
           "", "///|", "/// Model owning each agent: 0 K2, 1 Limite, 2 R1, 3 Qwen3.5, 4 LFM.",
           "pub let hexad_agent_model : Array[Int32] = ["]
    mid = {m: i for i, (m, _, _) in enumerate(ROLES)}
    ids = [str(mid[c["model"]]) for c in cards]
    for i in range(0, 128, 32):
        out.append("  " + ", ".join(ids[i:i + 32]) + ",")
    out += ["]", "", "///|", "/// 128 rows of 1024 signed bytes, two hex digits each.",
            "pub let hexad_lfm_hex : Array[String] = ["]
    out += [f'  "{h}",' for h in hexrows]
    out += ["]", ""]
    (ROOT / "src/26_lfm_agent_embeddings.mbt").write_text("\n".join(out), encoding="utf-8")
    print(json.dumps(meta, indent=1))
    # role-block similarity summary
    names = [m for m, _, _ in ROLES]; lab = np.array([mid[c["model"]] for c in cards])
    print("mean cosine between role blocks:")
    for i in range(5):
        print(f"  {names[i][:14]:14s}", " ".join(f"{S[np.ix_(lab==i, lab==j)].mean():.3f}" for j in range(5)))

if __name__ == "__main__":
    main()
