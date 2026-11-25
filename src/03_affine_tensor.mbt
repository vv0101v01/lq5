// ============================================================================
// 03_affine_tensor.mbt -- Affine tensor space A(T) and affine connections.
//
// ABSTRACTION.  A tensor of type (r,s) is a multilinear map, i.e. an element
//   T ∈ V^{⊗r} ⊗ (V*)^{⊗s}.
// On a computer it is (shape, variance, flat data): `variance[i] = +1`
// for upper (contravariant) and -1 for lower (covariant) indices.  An affine
// map Φ(x) = Ax + b acts by pushforward on upper and pullback on lower
// indices.  A constant A has vanishing derivative, hence Γ' = 0 and R' = 0;
// flatness is an affine invariant (Theorem T6).
// ============================================================================

///|
pub(all) struct Tensor {
  shape : Array[Int32]
  variance : Array[Int32] // +1 upper, -1 lower
  data : Array[Double] // row-major over the full multi-index
}

///|
pub fn Tensor::new(
  shape : Array[Int32],
  variance : Array[Int32],
  data : Array[Double],
) -> Tensor {
  Tensor::{ shape, variance, data, }
}

///|
pub fn Tensor::from_matrix(
  a : Array[Double],
  n : Int32,
  var0 : Int32,
  var1 : Int32,
) -> Tensor {
  Tensor::{ shape: [n, n], variance: [var0, var1], data: a, }
}

///|
pub fn Tensor::size(self : Tensor) -> Int32 {
  self.data.length()
}

///|
pub fn Tensor::index(self : Tensor, idx : Array[Int32]) -> Int32 {
  let mut k = 0
  let mut stride = 1
  let mut i = self.shape.length() - 1
  while i >= 0 {
    k += idx[i] * stride
    stride *= self.shape[i]
    i -= 1
  }
  k
}

///|
pub fn Tensor::at(self : Tensor, idx : Array[Int32]) -> Double {
  self.data[self.index(idx)]
}

///|
/// Contraction of a rank-2 tensor with one upper and one lower index.
/// Returns a rank-0 tensor (scalar).  This is the trace of the associated
/// endomorphism, an affine invariant (Theorem T6).
pub fn Tensor::contract_rank2(self : Tensor) -> Double {
  if self.shape.length() != 2 {
    return 0.0
  }
  if self.shape[0] != self.shape[1] {
    return 0.0
  }
  if self.variance[0] + self.variance[1] != 0 {
    return 0.0
  }
  let n = self.shape[0]
  let mut s = 0.0
  for k in 0..<n {
    s += self.data[k * n + k]
  }
  s
}

///|
/// Tensor (outer) product: (S ⊗ T)_{i...,j...} = S_{i...} T_{j...}.
pub fn Tensor::product(self : Tensor, other : Tensor) -> Tensor {
  let shape = self.shape.copy()
  shape.append(other.shape[:])
  let variance = self.variance.copy()
  variance.append(other.variance[:])
  let data = Array::makei(self.size() * other.size(), fn(k) {
    self.data[k / other.size()] * other.data[k % other.size()]
  })
  Tensor::{ shape, variance, data, }
}

///|
/// Affine map Φ(x) = A x + b on body coordinates.
pub fn affine_map(
  x : Array[Double],
  a : Array[Double],
  b : Array[Double],
  n : Int32,
) -> Array[Double] {
  let ax = matvec(a, x, n)
  Array::makei(n, fn(i) { ax[i] + b[i] })
}

///|
/// Pullback of a covariant 2-tensor: g' = Aᵀ g A.
/// For upper (contravariant) tensors it would be A^{-1} g A^{-T}.
pub fn pullback_metric(
  g : Array[Double],
  a : Array[Double],
  n : Int32,
) -> Array[Double] {
  let at = transpose(a, n, n)
  matrix_product(at, n, n, matrix_product(g, n, n, a, n), n)
}

///|
pub fn rotation(phi : Double) -> Array[Double] {
  let c = @math.cos(phi)
  let s = @math.sin(phi)
  [c, -s, s, c]
}

///|
/// Γ'^k_{ij} = (A^{-1})^k_l A^p_i A^q_j Γ^l_{pq} + (A^{-1})^k_l ∂_i A^l_j.
/// For a CONSTANT A the derivative term vanishes.  Stored as a flat array
/// gamma[(i*n + j)*n + k].
pub fn transformed_connection(
  gamma : Array[Double],
  a : Array[Double],
  n : Int32,
) -> Array[Double] {
  let ai = inverse(a, n)
  Array::makei(n * n * n, fn(t) {
    let k = t % n
    let ij = t / n
    let i = ij / n
    let j = ij % n
    let mut s = 0.0
    for l in 0..<n {
      for p in 0..<n {
        for q in 0..<n {
          s += ai[k * n + l] *
            a[p * n + i] *
            a[q * n + j] *
            gamma[(p * n + q) * n + l]
        }
      }
    }
    s
  })
}

///|
/// Action of the connection on a contravariant vector:
/// (∇_v X)^k = Γ^k_{ij} v^i X^j  (constant field ⇒ ∂X = 0).
pub fn connection_action(
  gamma : Array[Double],
  v : Array[Double],
  x : Array[Double],
  n : Int32,
) -> Array[Double] {
  Array::makei(n, fn(k) {
    let mut s = 0.0
    for i in 0..<n {
      for j in 0..<n {
        s += gamma[(i * n + j) * n + k] * v[i] * x[j]
      }
    }
    s
  })
}

///|
/// Torsion T^k_{ij} = Γ^k_{ij} - Γ^k_{ji}; 0 for a symmetric connection.
pub fn torsion_norm(gamma : Array[Double], n : Int32) -> Double {
  let mut s = 0.0
  for i in 0..<n {
    for j in 0..<n {
      for k in 0..<n {
        let d = gamma[(i * n + j) * n + k] - gamma[(j * n + i) * n + k]
        if d.abs() > s {
          s = d.abs()
        }
      }
    }
  }
  s
}

///|
/// R^l_{ijk} = Γ^l_{ip}Γ^p_{jk} - Γ^l_{jp}Γ^p_{ik} (constant Γ).
pub fn curvature_norm(gamma : Array[Double], n : Int32) -> Double {
  let mut s = 0.0
  for i in 0..<n {
    for j in 0..<n {
      for k in 0..<n {
        for l in 0..<n {
          let mut r = 0.0
          for p in 0..<n {
            r += gamma[(i * n + p) * n + l] * gamma[(j * n + k) * n + p] -
              gamma[(j * n + p) * n + l] * gamma[(i * n + k) * n + p]
          }
          if r.abs() > s {
            s = r.abs()
          }
        }
      }
    }
  }
  s
}

///|
/// Parallel transport by explicit Euler along a constant direction v.
pub fn parallel_transport(
  x : Array[Double],
  gamma : Array[Double],
  v : Array[Double],
  dt : Double,
  steps : Int32,
) -> Array[Double] {
  let cur = x.copy()
  for _ in 0..<steps {
    let d = connection_action(gamma, v, cur, x.length())
    for k in 0..<cur.length() {
      cur[k] -= dt * d[k]
    }
  }
  cur
}

///|
test "T6: contraction is the invariant trace; tensor product Leibniz" {
  let id2 = Tensor::from_matrix(mat_identity(2), 2, 1, -1)
  assert_true(approx(id2.contract_rank2(), 2.0))
  let g = Tensor::from_matrix([2.0, 0.0, 0.0, 5.0], 2, -1, -1)
  let v = Tensor::new([2], [1], [1.0, 2.0])
  let prod = g.product(v)
  assert_eq(prod.size(), 8)
  assert_true(approx(prod.at([1, 1, 1]), 10.0))
  assert_eq(prod.index([1, 1, 0]), 6)
}

///|
test "T6: affine isometries and flatness are affine invariants" {
  let n = 2
  let id = mat_identity(n)
  let r = rotation(0.7)
  // rotation is an isometry of the Euclidean metric
  let g2 = pullback_metric(id, r, n)
  for k in 0..<4 {
    assert_true(approx(g2[k], id[k], tol=1.0e-12))
  }
  // affine map is genuinely affine
  let p = [1.0, 2.0]
  let q = affine_map(p, r, [3.0, -1.0], n)
  assert_true(approx(q[0], r[0] * 1.0 + r[1] * 2.0 + 3.0))
  // flat connection stays flat under any constant affine change
  let zero = Array::make(n * n * n, 0.0)
  let a = [2.0, 0.5, 0.0, 1.5]
  let shifted = transformed_connection(zero, a, n)
  assert_eq(curvature_norm(shifted, n), 0.0)
  assert_eq(torsion_norm(shifted, n), 0.0)
  let x = [1.0, 0.0]
  let v = [0.0, 1.0]
  let out = parallel_transport(x, zero, v, 0.01, 10)
  assert_true(approx(out[0], 1.0))
}

///|
test "T6: symmetric connection is torsion-free" {
  let n = 2
  let gamma = Array::makei(n * n * n, fn(t) {
    let i = t / (n * n)
    let j = t / n % n
    let k = t % n
    if i == j {
      0.5 * (k + 1).to_double()
    } else {
      0.25
    }
  })
  assert_true(approx(torsion_norm(gamma, n), 0.0))
  // Leibniz rule for the connection action on a product direction
  let v = [1.0, 0.5]
  let x = [2.0, -1.0]
  let dx = connection_action(gamma, v, x, n)
  let scaled = connection_action(gamma, v, vec_scale(x, 2.0), n)
  for k in 0..<n {
    assert_true(approx(scaled[k], 2.0 * dx[k], tol=1.0e-12))
  }
}
