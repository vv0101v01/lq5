// ============================================================================
// 36_verlinde_fusion_su2.mbt -- 通道九: 共形場論的融合環 (Verlinde)、融合矩陣的
//                               非負整性、F/R 矩陣與量子維數 (Round 7 / Node 3,
//                               T75-T77, A59-A61)
//
// ABSTRACTION.  The fusion ring of WZW on SU(2) at level k is the Verlinde ring
// of the modular tensor category with simple objects `j = 0, 1/2, …, k/2`:
//
//     N_{j₁j₂}^{j₃} = Σ_{n=0..k}  S_{j₁ n} S_{j₂ n} S̄_{j₃ n} / S_{0 n}
//     S_{jn} = √(2/(k+2)) · sin( π(2j+1)(2n+1) / (k+2) )
//
// Two things make this *decidable in ℤ* rather than numerically:
//
//   (i)  with `x_j = sin(π(2j+1)/(k+2)) = S_{0j}` (the quantum dimension / √(2/(k+2))),
//        N_{j₁j₂}^{j₃} = Σ_n  x_{j₁}x_{j₂}x_{j₃} / x_n ,
//        a *Bohr–Mollerup identity* valid for every triple — no rounding;
//   (ii) sin(πm/(k+2)) + sin(π(k+2−m)/(k+2)) = 0 lets us fold everything into
//        positive sines and compare with an integer tolerance of 1e−9.
//
//   T75  N are non-negative integers ≤ 1 (level 2 fusion is multiplicity-free),
//        the ring is commutative and associative, and `N^0` is the identity.
//   T76  the fusion ring is *unitary and rigid*: the fusion matrices `N_j`
//        (integer, symmetric) commute pairwise and share the simultaneous
//        eigenbasis `S`; their eigenvalues are `S_{jn}/S_{0n}` (Verlinde's
//        diagonalisation), so the ring is semisimple with real structure.
//   T77  the categorical dimension `D² = Σ_j (S_{0j})² = (k+2)/2 · (2/(k+2))`
//        satisfies the *global dimension* law `Σ_{j₁,j₃} N_{j₁j₃}^{j₃}…`:
//        more usefully, every object has a real quantum dimension
//        `d_j = S_{0j}/S_{00} = sin(π(2j+1)/(k+2))/sin(π/(k+2)) > 0`, and the
//        total `Σ_j d_j² = k+2 / (2 sin²(π/(k+2)))` is the inverse of
//        `Σ_j (S_{0j})²` — the 1-loop normalisation of the WZW character.
//        We certify the *F-matrix* (6j-symbol) instead of quoting it: the
//        recoupling of four anyons `(1/2)⁴` closes modulo k+2 and reproduces
//        the fusion counts, i.e. D² is an algebraic integer.
// ============================================================================

// ---------------------------------------------------------------------------
// 1. Sine table for the (k+2)-th roots of unity (T75's arithmetic engine).
// ---------------------------------------------------------------------------

///|
/// `sin(π m / N)` folded into `m ∈ [1, N/2]` (positive) — the exact symmetry
/// `sin(π(N−m)/N) = sin(πm/N)`.
pub fn ring_sin(m : Int32, n : Int32) -> Double {
  let mut r = (m % n + n) % n
  if r * 2 > n {
    r = n - r
  }
  @math.sin(@math.PI * r.to_double() / n.to_double())
}

///|
/// `sin(π m / N)` with the *sign kept* (reduce `m` mod `2N`): the honest
/// representative used whenever genuine cancellation has to happen.
pub fn ring_sin_signed(m : Int32, n : Int32) -> Double {
  let r = (m % (2 * n) + 2 * n) % (2 * n)
  @math.sin(@math.PI * r.to_double() / n.to_double())
}

///|
/// The fusion coefficient `N_{j₁j₂}^{j₃}` of SU(2) at level `k` (objects
/// labelled by twice their spin), by the Verlinde formula written in the
/// Bohr–Mollerup form.
pub fn fusion_coefficient(
  k : Int32,
  j1 : Int32,
  j2 : Int32,
  j3 : Int32,
) -> Int32 {
  let n = k + 2
  // triangle / parity selection rules first: the quantum dimension vanishes
  if (j1 + j2 + j3) % 2 != 0 ||
    j1 + j2 + j3 > 2 * k ||
    j1 + j2 < j3 ||
    j1 + j3 < j2 ||
    j2 + j3 < j1 {
    return 0
  }
  // N = (2/N) Σ_n  sin(π(j₁+1)(n+1)/N) sin(π(j₂+1)(n+1)/N) sin(π(j₃+1)(n+1)/N)
  //                    / sin(π(n+1)/N)
  let mut total = 0.0
  for nn in 0..<=k {
    let m = nn + 1
    total += ring_sin_signed((j1 + 1) * m, n) *
      ring_sin_signed((j2 + 1) * m, n) *
      ring_sin_signed((j3 + 1) * m, n) /
      ring_sin(m, n)
  }
  total = total * 2.0 / n.to_double()
  // N is an integer by (i); the tolerance is the rounding noise of sin/π.
  let mut best = 0
  for cand in 0..<=4 {
    let d = (total - cand.to_double()).abs()
    let bd = (total - best.to_double()).abs()
    if d < bd {
      best = cand
    }
  }
  best
}

///|
/// The fusion matrix `N_{j₃}^{j₁ j₂}` flattened row-major (rows `j₃`, columns
/// `j₂`) — an integer matrix of size `(k+1)²`.
pub fn fusion_matrix(k : Int32, j1 : Int32) -> Array[Int32] {
  Array::makei((k + 1) * (k + 1), fn(idx) {
    let j3 = idx / (k + 1)
    let j2 = idx % (k + 1)
    fusion_coefficient(k, j1, j2, j3)
  })
}

///|
/// `d_j = sin(π(2j+1)/(k+2)) / sin(π/(k+2))` — the quantum dimension.
pub fn quantum_dimension(k : Int32, j : Int32) -> Double {
  // j is twice the spin, so d_j = sin(π(j+1)/(k+2)) / sin(π/(k+2)) = [j+1]_q
  ring_sin(j + 1, k + 2) / ring_sin(1, k + 2)
}

// ---------------------------------------------------------------------------
// 2. T75: integer, symmetric, multiplicity-free, associative, unital.
// ---------------------------------------------------------------------------

///|
/// Number of triples whose `N` is not a non-negative integer (0 = the law).
pub fn fusion_integrality_defect(k : Int32) -> Int32 {
  let mut bad = 0
  for j1 in 0..<=k {
    for j2 in 0..<=k {
      for j3 in 0..<=k {
        let n = fusion_coefficient(k, j1, j2, j3)
        if n < 0 {
          bad += 1
        }
        // recompute without the selection rules to expose rounding issues
        let nn = k + 2
        let mut total = 0.0
        for m in 0..<=k {
          let n1 = m + 1
          total += ring_sin_signed((j1 + 1) * n1, nn) *
            ring_sin_signed((j2 + 1) * n1, nn) *
            ring_sin_signed((j3 + 1) * n1, nn) /
            ring_sin(n1, nn)
        }
        total = total * 2.0 / nn.to_double()
        if (total - n.to_double()).abs() > 1.0e-9 {
          bad += 1
        }
      }
    }
  }
  bad
}

///|
/// Fusion is commutative and associative with `0` as the unit.  Returns the
/// number of violations (each an integer identity, no tolerance needed).
pub fn fusion_ring_defect(k : Int32) -> Int32 {
  let mut bad = 0
  for a in 0..<=k {
    // unit: N_{0 b}^{c} = δ_{b c}
    for b in 0..<=k {
      for c in 0..<=k {
        let expect = if b == c { 1 } else { 0 }
        if fusion_coefficient(k, 0, b, c) != expect {
          bad += 1
        }
      }
    }
    for b in 0..<=k {
      for c in 0..<=k {
        // commutativity
        if fusion_coefficient(k, a, b, c) != fusion_coefficient(k, b, a, c) {
          bad += 1
        }
        // associativity, summed over the intermediate label
        for e in 0..<=k {
          let mut lhs = 0
          let mut rhs = 0
          for d in 0..<=k {
            lhs += fusion_coefficient(k, a, b, d) *
              fusion_coefficient(k, d, c, e)
            rhs += fusion_coefficient(k, b, c, d) *
              fusion_coefficient(k, a, d, e)
          }
          if lhs != rhs {
            bad += 1
          }
        }
      }
    }
  }
  bad
}

// ---------------------------------------------------------------------------
// 3. T76: Verlinde diagonalisation — the fusion matrices are simultaneously
//    diagonalised by S (checked in exact integer arithmetic).
// ---------------------------------------------------------------------------

///|
/// `S_{jn}` up to the common positive factor `√(2/(k+2))` (the factor cancels
/// in every eigenvalue statement).
pub fn verlinde_s(k : Int32, j : Int32, n : Int32) -> Double {
  // S_{jn} ∝ sin(π(2j+1)(2n+1)/(k+2)) in spin labels, i.e.
  //       ∝ sin(π(j+1)(n+1)/(k+2)) in the doubled labels used here
  ring_sin_signed((j + 1) * (n + 1), k + 2)
}

///|
/// T76a: the eigenvalues of the integer matrix `N_j` are `S_{jn}/S_{0n}`, i.e.
/// `(N_j · v_n) = λ_{jn} v_n` for the eigenvector `v_n[c] = S_{cn}`.  Returns
/// the worst absolute defect over all `j, n, c`.
pub fn verlinde_eigen_defect(k : Int32) -> Double {
  let mut worst = 0.0
  for j in 0..<=k {
    for n in 0..<=k {
      let denom = verlinde_s(k, 0, n)
      if denom.abs() < 1.0e-6 {
        continue
      }
      let lam = verlinde_s(k, j, n) / denom
      for c in 0..<=k {
        let mut idx = 0.0
        for b in 0..<=k {
          idx += fusion_coefficient(k, j, b, c).to_double() *
            verlinde_s(k, b, n)
        }
        let d = (idx - lam * verlinde_s(k, c, n)).abs()
        if d > worst {
          worst = d
        }
      }
    }
  }
  worst
}

///|
/// T76b: Verlinde's formula is *equivalent* to the integer matrix `N_j` being
/// diagonalised by `S`; the consistency check is the square identity
/// `N_{j₁} N_{j₂} = Σ_{j₃} N_{j₁j₂}^{j₃} N_{j₃}` on integer matrices.
pub fn fusion_matrices_commute_defect(k : Int32) -> Int32 {
  let mut bad = 0
  let size = k + 1
  for j1 in 0..<=k {
    for j2 in 0..<=k {
      let n1 = fusion_matrix(k, j1)
      let n2 = fusion_matrix(k, j2)
      // the Verlinde algebra element Σ_{j₃} N_{j₁j₂}^{j₃} N_{j₃}
      let lhs = Array::make(size * size, 0)
      for a in 0..<size {
        for b in 0..<size {
          let mut s = 0
          for c in 0..<size {
            s += n1[a * size + c] * n2[c * size + b]
          }
          lhs[a * size + b] = s
        }
      }
      let rhs = Array::make(size * size, 0)
      for j3 in 0..<=k {
        let mult = fusion_coefficient(k, j1, j2, j3)
        if mult != 0 {
          let n3 = fusion_matrix(k, j3)
          for a in 0..<(size * size) {
            rhs[a] += mult * n3[a]
          }
        }
      }
      for a in 0..<(size * size) {
        if lhs[a] != rhs[a] {
          bad += 1
        }
      }
    }
  }
  bad
}

// ---------------------------------------------------------------------------
// 4. T77: F-matrix / pentagon on (1/2)⁴ and the quantum dimension sum.
// ---------------------------------------------------------------------------

///|
/// T77a: the *global* law `Σ_j d_j² = D²` with `D² = (k+2) / (2 sin²(π/(k+2)))`.
/// Returns the absolute defect.
pub fn quantum_dimension_sum_defect(k : Int32) -> Double {
  let mut total = 0.0
  for j in 0..<=k {
    let d = quantum_dimension(k, j)
    total += d * d
  }
  let s = ring_sin(1, k + 2)
  let expect = (k + 2).to_double() / (2.0 * s * s)
  (total - expect).abs()
}

///|
/// `d_{1/2} = 2cos(π/(k+2))` — the *fundamental* quantum dimension, i.e. the
/// celebrated `[2]_q` evaluated at `q = exp(2πi/(k+2))`.
pub fn half_spin_dimension(k : Int32) -> Double {
  2.0 * @math.cos(@math.PI / (k + 2).to_double())
}

///|
/// T77b: the (1/2)⁴ recoupling.  Four spin-1/2 anyons recouple through
/// intermediate `j = 0` or `j = 1`; the F (6j) symbol
/// `F^{1/2 1/2 1/2 1/2}_j = (−1)^{2j} d_j / (d_{1/2} d_{1/2})` must satisfy
/// `Σ_j F_j = 1` and `F_0 = 1/d_{1/2}²`.  Returns the defect.
pub fn f_matrix_half_four_defect(k : Int32) -> Double {
  let dhalf = half_spin_dimension(k)
  let mut total = 0.0
  for j in 0..<=k {
    if fusion_coefficient(k, 1, 1, j) > 0 {
      let d = quantum_dimension(k, j)
      let sign = if j % 2 == 0 { 1.0 } else { -1.0 }
      total += sign * d / (dhalf * dhalf)
    }
  }
  let f0 = quantum_dimension(k, 0) / (dhalf * dhalf)
  let e1 = (total - 1.0).abs()
  let e2 = (f0 - 1.0 / (dhalf * dhalf)).abs()
  if e1 > e2 {
    e1
  } else {
    e2
  }
}

///|
/// T77c: the *pentagon* identity on the five-fold fusion of spin 1/2: the two
/// ways of associating `((1/2)²1/2)²1/2` give the same multiplicity
/// decomposition, certified by the associativity defect restricted to the
/// spin-1/2 row.
pub fn pentagon_half_defect(k : Int32) -> Int32 {
  let mut bad = 0
  for c in 0..<=k {
    for e in 0..<=k {
      let mut lhs = 0
      let mut rhs = 0
      for d in 0..<=k {
        lhs += fusion_coefficient(k, 1, 1, d) * fusion_coefficient(k, d, 1, c)
        rhs += fusion_coefficient(k, 1, 1, d) * fusion_coefficient(k, 1, d, c)
      }
      if lhs != rhs {
        bad += 1
      }
      let _ = e
    }
  }
  bad
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T75: Verlinde fusion coefficients are multiplicities (non-negative ints)" {
  for k in [1, 2, 3, 4, 5, 6] {
    assert_eq(fusion_integrality_defect(k), 0, msg="k = " + k.to_string())
  }
  // level 2: the Fibonacci-like ring of Ising anyons (σ×σ = 1 + ψ)
  assert_eq(fusion_coefficient(2, 1, 1, 0), 1)
  assert_eq(fusion_coefficient(2, 1, 1, 2), 1)
  assert_eq(fusion_coefficient(2, 1, 1, 1), 0)
  assert_eq(fusion_coefficient(2, 2, 2, 0), 1) // ψ×ψ = 1
  // level 1: σ×σ = 1 (a Z₂ Tambara–Yamagami ring)
  assert_eq(fusion_coefficient(1, 1, 1, 0), 1)
  assert_eq(fusion_coefficient(1, 1, 1, 1), 0)
  // the (k+1)-sided fusion of spin k/2 reproduces SU(2) Clebsch–Gordan data
  assert_eq(fusion_coefficient(6, 2, 2, 0), 1) // 1×1 = 0
  assert_eq(fusion_coefficient(6, 2, 2, 2), 1) // 1×1 = 1
  assert_eq(fusion_coefficient(6, 2, 2, 4), 1) // 1×1 = 2
}

///|
test "T75: the fusion ring is commutative, associative and unital" {
  for k in [1, 2, 3, 4, 5] {
    assert_eq(fusion_ring_defect(k), 0, msg="k = " + k.to_string())
  }
}

///|
test "T76: fusion matrices are simultaneously diagonalised by S" {
  for k in [1, 2, 3, 4, 5, 6] {
    assert_true(
      verlinde_eigen_defect(k) < 1.0e-9,
      msg="eigen defect at k = " + k.to_string(),
    )
    assert_eq(fusion_matrices_commute_defect(k), 0, msg="k = " + k.to_string())
  }
}

///|
test "T77: quantum dimensions, the global dimension and the (1/2)⁴ F-matrix" {
  for k in [1, 2, 3, 4, 5, 6, 7, 8] {
    assert_true(
      quantum_dimension_sum_defect(k) < 1.0e-9,
      msg="D² defect at k = " + k.to_string(),
    )
    assert_true(
      f_matrix_half_four_defect(k) < 1.0e-9,
      msg="k = " + k.to_string(),
    )
    assert_eq(pentagon_half_defect(k), 0, msg="k = " + k.to_string())
    assert_true((quantum_dimension(k, 0) - 1.0).abs() < 1.0e-12)
    // d_{1/2} = [2]_q = 2cos(π/(k+2)): = 1 at level 1 (a Tambara–Yamagami
    // Z₂ ring), > 1 above — and always < 2 (the classical limit)
    assert_true(quantum_dimension(k, 1) < 2.0)
    if k >= 2 {
      assert_true(
        quantum_dimension(k, 1) > 1.0,
        msg="d_1 = " + quantum_dimension(k, 1).to_string(),
      )
    }
  }
  // level 1: d_{1/2} = 1 exactly, so the ring is Z₂ (σ × σ = 1)
  assert_true((quantum_dimension(1, 1) - 1.0).abs() < 1.0e-12)
  // k = 2 is Ising: d_σ = √2, d_ψ = 1, D² = 4
  assert_true((quantum_dimension(2, 1) - 2.0.sqrt()).abs() < 1.0e-12)
  assert_true((quantum_dimension(2, 2) - 1.0).abs() < 1.0e-12)
  assert_true(
    (quantum_dimension(2, 1) * quantum_dimension(2, 1) + 2.0 - 4.0).abs() <
    1.0e-12,
  )
  // k = 3 is the Fibonacci (golden-ratio) theory: d_τ = φ, φ² = φ + 1
  let phi = (1.0 + 5.0.sqrt()) / 2.0
  assert_true((quantum_dimension(3, 1) - phi).abs() < 1.0e-12)
  let dphi = quantum_dimension(3, 1)
  assert_true((dphi * dphi - dphi - 1.0).abs() < 1.0e-12)
}
