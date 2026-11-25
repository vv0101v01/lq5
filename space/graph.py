"""frame1-Hexad coupling graphs in numpy (port of frame1 src/25_knn_graph.mbt, A39).

ring  : the Round-1 N-cycle, lambda2 = 2 - 2cos(2pi/N)
knn   : symmetric k-NN on cosine, heat-kernel w = exp((cos-1)/tau),
        + one bridge per agent to its most similar agent of every other family,
        + Kruskal repair until connected.
"""
from __future__ import annotations

import numpy as np


def cosine_matrix(emb: np.ndarray) -> np.ndarray:
    x = np.asarray(emb, dtype=np.float64)
    x = x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-300)
    return x @ x.T


def components(w: np.ndarray) -> np.ndarray:
    n = w.shape[0]
    label = -np.ones(n, dtype=int)
    c = 0
    for s in range(n):
        if label[s] >= 0:
            continue
        stack = [s]
        label[s] = c
        while stack:
            u = stack.pop()
            for v in np.nonzero((w[u] > 0) & (label < 0))[0]:
                label[v] = c
                stack.append(int(v))
        c += 1
    return label


def knn_graph(emb, labels, k: int = 8, tau: float = 0.25, bridges: bool = True, repair: bool = True):
    sim = cosine_matrix(emb)
    n = sim.shape[0]
    labels = np.asarray(labels, dtype=int)
    w = np.zeros((n, n))

    def link(i, j):
        v = np.exp((sim[i, j] - 1.0) / tau)
        w[i, j] = v
        w[j, i] = v

    k = max(1, min(k, n - 1))
    for i in range(n):
        s = sim[i].copy()
        s[i] = -np.inf
        # stable argsort on -s == "ties broken by the smaller index"
        for j in np.argsort(-s, kind="stable")[:k]:
            link(i, int(j))
    if bridges:
        fams = int(labels.max()) + 1 if n else 0
        for i in range(n):
            for f in range(fams):
                if f == labels[i]:
                    continue
                idx = np.nonzero(labels == f)[0]
                if idx.size:
                    j = idx[np.argmax(sim[i, idx])]  # argmax returns first max = smaller index
                    link(i, int(j))
    if repair:
        comp = components(w)
        while comp.max() > 0:
            mask = comp[:, None] != comp[None, :]
            s = np.where(mask, sim, -np.inf)
            bi, bj = np.unravel_index(int(np.argmax(s)), s.shape)
            link(int(bi), int(bj))
            comp = components(w)
    return w


def ring_graph(n: int) -> np.ndarray:
    w = np.zeros((n, n))
    for i in range(n):
        w[i, (i + 1) % n] = w[(i + 1) % n, i] = 1.0
    return w


def laplacian(w: np.ndarray) -> np.ndarray:
    return np.diag(w.sum(1)) - w


def spectrum(w: np.ndarray) -> np.ndarray:
    return np.linalg.eigvalsh(laplacian(w))


def graph_stats(w: np.ndarray, labels=None, kappa: float = 0.7) -> dict:
    ev = spectrum(w)
    n = w.shape[0]
    e = int(np.count_nonzero(np.triu(w, 1)))
    out = {
        "agents": n,
        "edges": e,
        "components": int(components(w).max() + 1),
        "lambda2": float(ev[1]) if n > 1 else 0.0,
        "lambda_max": float(ev[-1]),
        "heat_CFL_dt": float(2.0 / (kappa * ev[-1])) if ev[-1] > 0 else float("inf"),
    }
    if labels is not None:
        labels = np.asarray(labels)
        iu = np.triu_indices(n, 1)
        cross = (w[iu] > 0) & (labels[iu[0]] != labels[iu[1]])
        out["cross_family_edges"] = int(cross.sum())
    return out


def diffuse(w: np.ndarray, x0: np.ndarray, steps: int, eta: float = 0.5) -> np.ndarray:
    """Random-walk heat flow x <- (1-eta) x + eta * D^-1 W x (stable for eta <= 1 on any graph).
    Returns the per-step dispersion (mean distance of agents to the population mean)."""
    p = w / np.maximum(w.sum(1, keepdims=True), 1e-300)
    x = np.array(x0, dtype=np.float64)
    out = np.empty(steps + 1)
    for t in range(steps + 1):
        out[t] = float(np.linalg.norm(x - x.mean(0), axis=1).mean())
        if t < steps:
            x = (1.0 - eta) * x + eta * (p @ x)
    return out


def kuramoto(w: np.ndarray, steps: int, coupling: float = 3.2, dt: float = 0.05, seed: int = 0) -> np.ndarray:
    """Kuramoto on the graph, frame1 normalisation: theta_i' = omega_i + (K/n) sum_j n*w_ij/d_i sin(theta_j-theta_i)."""
    n = w.shape[0]
    rng = np.random.default_rng(seed)
    th = rng.uniform(0, 2 * np.pi, n)
    om = rng.normal(0, 0.1, n)
    m = w / np.maximum(w.sum(1, keepdims=True), 1e-300)

    def f(t):
        return om + coupling * (m * np.sin(t[None, :] - t[:, None])).sum(1)

    out = np.empty(steps + 1)
    for s in range(steps + 1):
        out[s] = float(np.abs(np.exp(1j * th).mean()))
        if s < steps:
            k1 = f(th)
            k2 = f(th + dt * k1)
            th = th + 0.5 * dt * (k1 + k2)
    return out


def pca2(emb: np.ndarray) -> np.ndarray:
    x = np.asarray(emb, dtype=np.float64)
    x = x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-300)
    x = x - x.mean(0)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    return u[:, :2] * s[:2]
