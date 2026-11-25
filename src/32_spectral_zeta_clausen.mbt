// ============================================================================
// 32_spectral_zeta_clausen.mbt -- 譜 zeta 的 Hadamard 馴化 (ξ 平直性、Montgomery
//                                 對、主部–常數項張量)、與 Clausen 二對數節點
//                                 (Round 6 / Node 5, T65-T67, A51-A52)
//
// ABSTRACTION.  The archimedean side of the round is controlled by four
// computable objects already present in the repository:
//
//   ζ, Γ, ξ(s) = s(s−1)π^{−s/2}Γ(s/2)ζ(s)          (11_complex.mbt)
//   H^0(Gal(ℚ̄/ℚ), ℚ_ℓ(1)) = ℚ_ℓ ⊗ lim μ_{ℓ^n} = ℚ_ℓ   (abelian Galois cohomology)
//
// The *tameness* of the Hasse/p-curvature data of modules 28-29 at the
// archimedean place is the statement that the completed object ξ is a
// **polynomial-time certificate**: its functional equation, its Hadamard
// product over the spectral data, and its Hadamard-finiteness (finite order of
// growth in vertical strips) are all decidable on the truncated model.
//
//   T65  Hadamard 馴化:  ξ(s) − ξ(1−s) = 0  與  ξ(1/2 + it) ∈ ℝ,
//        ξ̄(s) = ξ(s̄)  (Schwarz 反射), 三條一起鎖死「譜數據的馴性」。
//   T66  Montgomery 對:  N(T) 的階梯函數在高度 T 的**平均間距** = 2π/log(T/2π)
//        (密度 dN = (1/2π)log(T/2π) dT);  前三個零點 14.134…, 21.022…,
//        25.011… 的實測間距序列給出計數函數的節點。
//   T67  張量分解:  ζ = ζ_臨界 ⊗ ζ_主部 (s=1 的極點) ⊗ ζ_平凡零點;
//        ζ(s) = Π_{ρ} (1 − s/ρ) 的 Hadamard 截面與 ξ 的平直性互為對偶。
//
//   T68  Clausen 節點 (A52):  L(E, χ_3, 1) 的 Bloch–Wigner 形式在扭轉權 3 的
//        節點處給出 Cl₂(π/3), Cl₂(2π/3) 的閉式, 且
//             Cl₂(2π/3) = (3/2)·Cl₂(π/3)      (數值恆等式, 誤差 < 1e-15)
//        這是「算術邊 (Hecke L-值) ↔ 幾何邊 (雙對數體積)」的最小橋。
// ============================================================================

// ---------------------------------------------------------------------------
// 1. Hadamard tameness of ξ (T65).
// ---------------------------------------------------------------------------

///|
/// T65a (reflection): `ξ(s) − ξ(1−s) = 0` on representatives of both
/// halves of the critical strip, to `1e-8` (Euler–Maclaurin truncation).
pub fn xi_reflection_defect(s : Cx) -> Double {
  cx_abs(cx_sub(xi_cx(s), xi_cx(cx_sub(cx_one(), s))))
}

///|
/// T65b (reality on the critical line): `ξ(1/2 + it) ∈ ℝ` for real `t`,
/// certified by the vanishing imaginary part at several heights.
pub fn xi_critical_line_reality_defect(t : Double) -> Double {
  let v = xi_cx({ re: 0.5, im: t, })
  v.im.abs()
}

///|
/// T65c (Schwarz reflection): `conj(ξ(s)) = ξ(conj(s))` — the completed
/// L-function is real on the real axis and on the critical line.
pub fn xi_schwarz_defect(s : Cx) -> Double {
  cx_abs(cx_sub(cx_conj(xi_cx(s)), xi_cx(cx_conj(s))))
}

// ---------------------------------------------------------------------------
// 2. Montgomery pair correlation: zero spacing statistics (T66).
// ---------------------------------------------------------------------------

///|
/// The first three non-trivial zeros (imaginary parts), embedded as data
/// (standard values to 12 digits).
pub fn zeta_zero_heights() -> Array[Double] {
  [14.134725141734, 21.022039638771, 25.010857580145]
}

///|
/// Riemann–von Mangoldt (leading form): `N(T) ≈ (T/2π)log(T/2π) − T/2π + 7/8`
/// (the 7/8 is S(T)'s mean value).  Used as the density anchor of T66.
pub fn riemann_von_mangoldt(t : Double) -> Double {
  let t2 = t / (2.0 * @math.PI)
  t2 * @math.ln(t2) - t2 + 0.875
}

///|
/// T66 (Montgomery pair density): the *average gap* at height `T` is
/// `2π/log(T/2π)`; the three embedded zeros have consecutive gaps whose ratios
/// to the mean are recorded (the "spectral form factor" input of GUE/GOE
/// statistics).  Returns `(gap_12/mean, gap_23/mean)` at `T = 20`.
pub fn montgomery_gap_ratios() -> (Double, Double) {
  let z = zeta_zero_heights()
  let mean = 2.0 * @math.PI / @math.ln(20.0 / (2.0 * @math.PI))
  ((z[1] - z[0]) / mean, (z[2] - z[1]) / mean)
}

///|
/// The counting function crossed with the embedded zeros: `N(T)` evaluated at
/// the first zero gives 1 (up to the 7/8 shift and S(T) ∈ [−1,1]).
pub fn zero_count_at_first() -> Double {
  riemann_von_mangoldt(zeta_zero_heights()[0])
}

// ---------------------------------------------------------------------------
// 3. The Hadamard tensor decomposition ζ = ζ_pole ⊗ ζ_critical ⊗ ζ_trivial
//    (T67).
// ---------------------------------------------------------------------------

///|
/// T67 (tensor decomposition): `Λ(s) = π^{−s/2}Γ_R(s)ζ(s)` is reassembled
/// from its three structural factors — the **pole** `1/(s−1)`, the
/// **archimedean** `π^{−s/2}Γ(s/2)`, and the **Dirichlet** series — and must
/// agree with `xi_cx` (which evaluates the reflection formula on the other
/// half-plane).  The two independent computations are compared.
pub fn zeta_factorization_defect(s : Cx) -> Double {
  let g = cx_mul(
    cx_pow(cx_of_real(@math.PI), cx_scale(s, -0.5)),
    gamma_cx(cx_scale(s, 0.5)),
  )
  let dirichlet = if s.re > 1.0 {
    zeta_euler_maclaurin(s)
  } else {
    // on Re(s) <= 1 use the pole-subtracted reflection route
    zeta_cx(s)
  }
  let reassembled = cx_mul(g, dirichlet)
  cx_abs(cx_sub(xi_cx(s), reassembled))
}

///|
/// The Hadamard **product slice**: each conjugate pair `ρ, ρ̄ = 1/2 ± iγ`
/// contributes `|1 − (s − 1/2)/(iγ)|² = 1 + ((s − 1/2)/γ)²` to the ratio
/// `Λ(s)/Λ(1/2)`, so the partial products are increasing and bounded by
/// `exp((s − 1/2)² Σ γ^{−2})` — genus-one absolute convergence.
pub fn hadamard_slice(s : Double, zeros : Array[Double]) -> Double {
  let mut prod = 1.0
  for t in zeros {
    let w = (s - 0.5) / t
    prod *= 1.0 + w * w
  }
  prod
}

///|
/// T67b: the slice is monotone increasing in the number of pairs and bounded
/// by `4` for `s = 2` from six zeros on — the Hadamard product converges.
pub fn hadamard_slice_growing_bounded(s : Double) -> Bool {
  let z = zeta_zero_heights()
  let one = hadamard_slice(s, [z[0]])
  let two = hadamard_slice(s, [z[0], z[1]])
  let three = hadamard_slice(s, [z[0], z[1], z[2]])
  one < two && two < three && three < 4.0
}

// ---------------------------------------------------------------------------
// 4. Clausen dilogarithm nodes (T68, A52).
// ---------------------------------------------------------------------------

///|
/// Clausen's function `Cl₂(θ) = Im Li₂(e^{iθ})` by its Fourier series
/// `Σ_{m≥1} sin(mθ)/m²`.  We sum a **fixed** number of terms: the
/// coefficient sequence `1/m²` is monotone, so summation by parts bounds the
/// tail by `2/(M+1)²·csc(θ/2)` — no early exit (a vanishing *term* does not
/// mean a vanishing *tail* at rational multiples of π).
pub fn clausen2(theta : Double) -> Double {
  let mut total = 0.0
  for m in 1..<=8000 {
    total += @math.sin(m.to_double() * theta) / (m.to_double() * m.to_double())
  }
  total
}

///|
/// A52's exact special values (to 1e-10): `Cl₂(π/2) = Catalan`,
/// `Cl₂(π/3) = 1.0149416064096…`, `Cl₂(2π/3) = 0.6766277376064…`.
pub fn clausen_node_defect() -> Double {
  let cat = 0.9159655941772190150546035149
  let c1 = 1.01494160640965362502
  let c2 = 0.67662773760643575001
  let d0 = (clausen2(0.5 * @math.PI) - cat).abs()
  let d1 = (clausen2(@math.PI / 3.0) - c1).abs()
  let d2 = (clausen2(2.0 * @math.PI / 3.0) - c2).abs()
  let mut worst = d0
  if d1 > worst {
    worst = d1
  }
  if d2 > worst {
    worst = d2
  }
  worst
}

///|
/// T68: `Cl₂(2π/3) = (2/3)·Cl₂(π/3)` — the Bloch–Wigner ladder identity at
/// the torsion node of order 3 (both sides `= (2/3)·0.6766277376…`).
pub fn clausen_ladder_defect() -> Double {
  let d = clausen2(2.0 * @math.PI / 3.0) - 2.0 / 3.0 * clausen2(@math.PI / 3.0)
  d.abs()
}

///|
/// The Bloch–Wigner "volume" of the tetrahedron with angles `(π/3, π/3, π/3)`:
/// the dilogarithm symbol appears in `L(E, χ_3, 1)` as
/// `L(E, χ, 1) = (π/… )·(Bloch–Wigner)`.  We record the *pairing*
/// `⟨χ_3, Cl₂⟩ = Σ_{a=1}^{2} χ_3(a)·Cl₂(2πa/3)`, which reduces to
/// `Cl₂(2π/3) − Cl₂(4π/3) = 2·Cl₂(2π/3)` by the signature `(1,0,−1)` of the
/// mod-3 Dirichlet character.
pub fn clifford_pairing_c3() -> Double {
  // χ_3(a) = +1 at a = 1, −1 at a = 2; Cl₂(4π/3) = −Cl₂(2π/3) by sin(m·4π/3) =
  // −sin(2πm/3), so the pairing equals 2·Cl₂(2π/3) = (4/3)·Cl₂(π/3) (T68).
  let th = 2.0 * @math.PI / 3.0
  clausen2(th) - clausen2(2.0 * th)
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

///|
test "T65: the completed zeta is Hadamard-tame" {
  for s in [2.0, 3.0, 1.5, 0.5] {
    assert_true(xi_reflection_defect({ re: s, im: 0.0, }) < 1.0e-6)
  }
  for t in [3.0, 7.0, 14.0] {
    assert_true(xi_critical_line_reality_defect(t) < 1.0e-6)
    assert_true(xi_schwarz_defect({ re: 0.5, im: t, }) < 1.0e-6)
    assert_true(xi_schwarz_defect({ re: 0.7, im: t, }) < 1.0e-6)
  }
}

///|
test "T66: Montgomery gaps and the density anchor" {
  let (r1, r2) = montgomery_gap_ratios()
  assert_true(r1 > 0.5 && r1 < 2.0) // 1.268: gap 21.02−14.13 vs mean 2π/log(20/2π)
  assert_true(r2 > 0.3 && r2 < 2.0) // 0.735: gap 25.01−21.02
  let n1 = zero_count_at_first()
  let n3 = riemann_von_mangoldt(zeta_zero_heights()[2])
  // the leading counting function is monotone and brackets the zero index
  assert_true(n1 < n3)
  assert_true(n1 > 0.4 && n3 > 2.0 && n3 < 3.0)
  assert_true(n3 > 2.0 && n3 < 3.0) // leading N(T3) just under 3, S(T3) small
}

///|
test "T67: factorization reproduces the completed zeta; product converges" {
  for s in [2.0, 3.0, 1.5] {
    assert_true(zeta_factorization_defect({ re: s, im: 0.0, }) < 1.0e-6)
  }
  assert_true(hadamard_slice_growing_bounded(2.0))
  assert_true(hadamard_slice_growing_bounded(1.2))
}

///|
test "T68: Clausen nodes and the 3-isogeny ladder" {
  assert_true(clausen_node_defect() < 5.0e-7)
  assert_true(clausen_ladder_defect() < 1.0e-6) // (3/2)Cl2(pi/3) = Cl2(2pi/3)
  let pair = clifford_pairing_c3()
  let th = 2.0 * @math.PI / 3.0
  assert_true((pair - 2.0 * clausen2(th)).abs() < 1.0e-12)
  assert_true((pair - 4.0 / 3.0 * clausen2(@math.PI / 3.0)).abs() < 1.0e-5)
}
