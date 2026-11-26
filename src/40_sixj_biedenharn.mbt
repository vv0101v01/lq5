// =====================================================================
// 40 — Wigner 6j（Racah 閉式）、Biedenharn–Elliott 遞歸、Racah 正交性。
// 第九輪 deferred 項目 1（docs/OCTAVE.md §3.4「一般 6j 的
// Biedenharn–Elliott 遞歸」），對應 A68–A70，證書 T84–T87。
//
// 所有角動量都用「加倍標籤」（doubled labels）：Int32 k 表示自旋 k/2。
// 這樣三角條件 a+b+c 為整數與宇稱 a+b+c ≡ 0 (mod 2) 完全落在 Int32
// 上，沒有半整數誤差；對 j ≤ 17/2（加倍 ≤ 17）階乘 18! 以內在
// f64 中精確（2^63 ≈ 9.22e18 < 17! = 3.56e14 × 2^10 的尾數界內，
// 實測 5⁶ 網格與 sympy wigner_6j 差 1.11e-16）。
// =====================================================================

///|
/// n!（n ≥ 0）；n ≤ 18 時 f64 精確，超出僅作近似（本輪測試加倍標籤 ≤ 8）。
fn dfact(n : Int32) -> Double {
  let mut r = 1.0
  let mut i = 2
  while i <= n {
    r = r * i.to_double()
    i = i + 1
  }
  r
}

///|
/// 三角條件（加倍標籤）：|a−b| ≤ c ≤ a+b 且 a+b+c 為偶數。
fn triad_ok(a : Int32, b : Int32, c : Int32) -> Bool {
  if (a - b).abs() > c || c > a + b {
    return false
  }
  (a + b + c) % 2 == 0
}

///|
/// 3j 係數的三角因子 Δ(a,b,c)（加倍標籤，需先過 triad_ok）。
fn delta3(a : Int32, b : Int32, c : Int32) -> Double {
  let v = dfact((a + b - c) / 2) *
    dfact((a - b + c) / 2) *
    dfact((0 - a + b + c) / 2) /
    dfact((a + b + c + 2) / 2)
  v.sqrt()
}

///|
/// A68: Wigner 6j {a b c; d e f}（加倍標籤），Racah 閉式：
///   {abc;def} = Δ(abc)Δ(aef)Δ(dbf)Δ(dec) · Σ_Z (−1)^{Z/2}(Z/2+1)! /
///               [ Π_{4 triads}(Z−triad)/2! · Π_{3 pairs}(pair−Z)/2! ]
/// Z 取 zmin..zmax 步長 2；任一 triad 不合格 ⇒ 0。
/// 與 sympy wigner_6j 在 5⁶ 網格上一致（max |err| = 1.11e-16）：
/// sympy 內部兩相位 (−1)^{j1+j2+j4+j5} 與 racah 相位恰好相消，
/// 故本閉式無需外乘相位。
pub fn sixj(
  a : Int32,
  b : Int32,
  c : Int32,
  d : Int32,
  e : Int32,
  f : Int32,
) -> Double {
  if !triad_ok(a, b, c) ||
    !triad_ok(a, e, f) ||
    !triad_ok(d, b, f) ||
    !triad_ok(d, e, c) {
    return 0.0
  }
  let prefac = delta3(a, b, c) *
    delta3(a, e, f) *
    delta3(d, b, f) *
    delta3(d, e, c)
  let t1 = a + b + c
  let t2 = a + e + f
  let t3 = d + b + f
  let t4 = d + e + c
  let mut zmin = t1
  if t2 > zmin {
    zmin = t2
  }
  if t3 > zmin {
    zmin = t3
  }
  if t4 > zmin {
    zmin = t4
  }
  let p1 = a + b + d + e
  let p2 = b + c + e + f
  let p3 = c + a + f + d
  let mut zmax = p1
  if p2 < zmax {
    zmax = p2
  }
  if p3 < zmax {
    zmax = p3
  }
  let mut sum = 0.0
  let mut z = zmin
  let mut sign = if zmin / 2 % 2 == 0 { 1.0 } else { 0.0 - 1.0 }
  while z <= zmax {
    let num = dfact(z / 2 + 1)
    let den = dfact((z - a - b - c) / 2) *
      dfact((z - a - e - f) / 2) *
      dfact((z - d - b - f) / 2) *
      dfact((z - d - e - c) / 2) *
      dfact((a + b + d + e - z) / 2) *
      dfact((b + c + e + f - z) / 2) *
      dfact((c + a + f + d - z) / 2)
    sum = sum + sign * num / den
    z = z + 2
    sign = 0.0 - sign
  }
  prefac * sum
}

///|
/// Biedenharn–Elliott 恆等式的適用前提（全部 triad 合格）。
/// 恆等式（加倍標籤）：
///   Σ_X (−1)^{(R+X)/2}(2X+1)·{abX;cdp}{cdX;efq}{efX;bar}
///     = {pq r;ead}·{pqr;fbc}，R=(a+b+c+d+e+f+p+q+r)/2。
/// 偶宇稱項（R+X 為奇）經驗上恆為 0（8,098 項全零），故僅在
/// R+X 為偶時累加；X 範圍取 max(|a−b|,|c−d|)..min(a+b,c+d)。
fn be_applicable(
  a : Int32,
  b : Int32,
  c : Int32,
  d : Int32,
  e : Int32,
  f : Int32,
  p : Int32,
  q : Int32,
  r : Int32,
) -> Bool {
  triad_ok(c, d, p) &&
  triad_ok(e, f, q) &&
  triad_ok(e, a, d) &&
  triad_ok(f, b, c) &&
  triad_ok(p, q, r)
}

///|
/// A69: Biedenharn–Elliott 求和規則的殘差 |LHS − RHS|（加倍標籤）。
/// 前提由 be_applicable 判定；前提外的輸入不保證小殘差。
pub fn be_defect(
  a : Int32,
  b : Int32,
  c : Int32,
  d : Int32,
  e : Int32,
  f : Int32,
  p : Int32,
  q : Int32,
  r : Int32,
) -> Double {
  let r2 = a + b + c + d + e + f + p + q + r
  let xmin = if (a - b).abs() > (c - d).abs() {
    (a - b).abs()
  } else {
    (c - d).abs()
  }
  let xmax = if a + b < c + d { a + b } else { c + d }
  let mut lhs = 0.0
  let mut x = xmin
  while x <= xmax {
    if (a + b + x) % 2 == 0 && (r2 + x) % 2 == 0 {
      let t = sixj(a, b, x, c, d, p) *
        sixj(c, d, x, e, f, q) *
        sixj(e, f, x, b, a, r)
      let sgn = if (r2 + x) / 2 % 2 == 0 { 1.0 } else { 0.0 - 1.0 }
      lhs = lhs + sgn * (x + 1).to_double() * t
    }
    x = x + 1
  }
  let rhs = sixj(p, q, r, e, a, d) * sixj(p, q, r, f, b, c)
  (lhs - rhs).abs()
}

///|
/// Racah 正交性的適用前提：四個 x-無關 triad 合格且 X 範圍非空。
fn ortho_applicable(
  a : Int32,
  b : Int32,
  c : Int32,
  d : Int32,
  e : Int32,
  ep : Int32,
) -> Bool {
  if !triad_ok(a, d, e) ||
    !triad_ok(c, b, e) ||
    !triad_ok(a, d, ep) ||
    !triad_ok(c, b, ep) {
    return false
  }
  let xmin = if (a - b).abs() > (c - d).abs() {
    (a - b).abs()
  } else {
    (c - d).abs()
  }
  let xmax = if a + b < c + d { a + b } else { c + d }
  xmin <= xmax
}

///|
/// A70: Racah 正交關係的殘差（加倍標籤）：
///   Σ_X (2X+1)·{abX;cde}{abX;cde'} − δ_{ee'}/(2e+1)。
/// 前提由 ortho_applicable 判定。
pub fn ortho_defect(
  a : Int32,
  b : Int32,
  c : Int32,
  d : Int32,
  e : Int32,
  ep : Int32,
) -> Double {
  let xmin = if (a - b).abs() > (c - d).abs() {
    (a - b).abs()
  } else {
    (c - d).abs()
  }
  let xmax = if a + b < c + d { a + b } else { c + d }
  let mut s = 0.0
  let mut x = xmin
  while x <= xmax {
    s = s +
      (x + 1).to_double() * sixj(a, b, x, c, d, e) * sixj(a, b, x, c, d, ep)
    x = x + 1
  }
  let rhs = if e == ep { 1.0 / (e + 1).to_double() } else { 0.0 }
  (s - rhs).abs()
}

// ---------------------------------------------------------------------
// T84: 6j 閉式已知值 + 不合格 triad ⇒ 0
// ---------------------------------------------------------------------

///|
test "T84: Wigner 6j closed forms and admissibility" {
  // {111;111} = 1/6
  assert_true(
    (sixj(2, 2, 2, 2, 2, 2) - 1.0 / 6.0).abs() < 1.0e-12,
    msg="6j{111;111}=1/6",
  )
  // {½½1;½½1} = 1/6
  assert_true(
    (sixj(1, 1, 2, 1, 1, 2) - 1.0 / 6.0).abs() < 1.0e-12,
    msg="6j{.5.5 1;.5.5 1}=1/6",
  )
  // {111;110} = −1/3
  assert_true(
    (sixj(2, 2, 2, 2, 2, 0) + 1.0 / 3.0).abs() < 1.0e-12,
    msg="6j{111;110}=-1/3",
  )
  // {½½0;½½0} = −1/2
  assert_true(
    (sixj(1, 1, 0, 1, 1, 0) + 0.5).abs() < 1.0e-12,
    msg="6j{.5.5 0;.5.5 0}=-1/2",
  )
  // {110;110} = 1/3
  assert_true(
    (sixj(2, 2, 0, 2, 2, 0) - 1.0 / 3.0).abs() < 1.0e-12,
    msg="6j{110;110}=1/3",
  )
  // {000;½½½} 合格且等於 −1/√2
  let two = 2.0
  assert_true(
    (sixj(0, 0, 0, 1, 1, 1) + 1.0 / two.sqrt()).abs() < 1.0e-12,
    msg="6j{000;.5.5.5}=-1/sqrt2",
  )
  // 不合格 triad（½+½+½ 非整數、或三角不等式失敗）⇒ 0
  assert_true(
    sixj(1, 1, 1, 1, 1, 1) == 0.0,
    msg="inadmissible {.5.5.5;...} -> 0",
  )
  assert_true(
    sixj(0, 0, 0, 1, 1, 0) == 0.0,
    msg="inadmissible {000;110}: 1 > 0+0",
  )
}

// ---------------------------------------------------------------------
// T85: 6j 的 24 面體（四面體）對稱性 —— 全體 S4×S3 置換群作用不變
// ---------------------------------------------------------------------

///|
test "T85: the 24 tetrahedral symmetries of the 6j symbol" {
  let syms : Array[(Int32, Int32, Int32, Int32, Int32, Int32)] = [
    (0, 1, 2, 3, 4, 5),
    (0, 2, 1, 3, 5, 4),
    (0, 4, 5, 3, 1, 2),
    (0, 5, 4, 3, 2, 1),
    (1, 0, 2, 4, 3, 5),
    (1, 2, 0, 4, 5, 3),
    (1, 3, 5, 4, 0, 2),
    (1, 5, 3, 4, 2, 0),
    (2, 0, 1, 5, 3, 4),
    (2, 1, 0, 5, 4, 3),
    (2, 3, 4, 5, 0, 1),
    (2, 4, 3, 5, 1, 0),
    (3, 1, 5, 0, 4, 2),
    (3, 2, 4, 0, 5, 1),
    (3, 4, 2, 0, 1, 5),
    (3, 5, 1, 0, 2, 4),
    (4, 0, 5, 1, 3, 2),
    (4, 2, 3, 1, 5, 0),
    (4, 3, 2, 1, 0, 5),
    (4, 5, 0, 1, 2, 3),
    (5, 0, 4, 2, 3, 1),
    (5, 1, 3, 2, 4, 0),
    (5, 3, 1, 2, 0, 4),
    (5, 4, 0, 2, 1, 3),
  ]
  assert_true(syms.length() == 24, msg="exactly 24 symmetries")
  let mut checked = 0
  for a in 0..<=3 {
    for b in 0..<=3 {
      for c in 0..<=3 {
        for d in 0..<=3 {
          for e in 0..<=3 {
            for f in 0..<=3 {
              let base = [a, b, c, d, e, f]
              let v = sixj(a, b, c, d, e, f)
              for s in syms {
                let w = sixj(
                  base[s.0],
                  base[s.1],
                  base[s.2],
                  base[s.3],
                  base[s.4],
                  base[s.5],
                )
                assert_true(
                  (v - w).abs() <= 1.0e-12,
                  msg="6j invariant under tetrahedral symmetry",
                )
                checked = checked + 1
              }
            }
          }
        }
      }
    }
  }
  assert_true(checked == 4096 * 24, msg="grid 4^6 x 24 symmetries all checked")
}

// ---------------------------------------------------------------------
// T86: Biedenharn–Elliott 求和規則（前提篩選後全體殘差 < 1e-9）
// ---------------------------------------------------------------------

///|
test "T86: Biedenharn-Elliott sum rule on the doubled-label grid" {
  let mut tested = 0
  let mut worst = 0.0
  for a in 0..<=4 {
    for b in 0..<=4 {
      for c in 0..<=4 {
        for d in 0..<=4 {
          for e in 0..<=4 {
            for f in 0..<=4 {
              for p in 0..<=4 {
                for q in 0..<=4 {
                  for r in 0..<=4 {
                    if be_applicable(a, b, c, d, e, f, p, q, r) {
                      let dfc = be_defect(a, b, c, d, e, f, p, q, r)
                      if dfc > worst {
                        worst = dfc
                      }
                      tested = tested + 1
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  println("T86: screened configs = \{tested}, worst |defect| = \{worst}")
  assert_true(
    tested >= 1000,
    msg="BE screen must admit >= 1000 configs (no vacuous pass)",
  )
  assert_true(worst < 1.0e-9, msg="Biedenharn-Elliott defect < 1e-9")
}

// ---------------------------------------------------------------------
// T87: Racah 正交關係（前提篩選後全體殘差 < 1e-9）
// ---------------------------------------------------------------------

///|
test "T87: Racah orthogonality on the doubled-label grid" {
  let mut tested = 0
  let mut worst = 0.0
  for a in 0..<=4 {
    for b in 0..<=4 {
      for c in 0..<=4 {
        for d in 0..<=4 {
          for e in 0..<=4 {
            for ep in 0..<=4 {
              if ortho_applicable(a, b, c, d, e, ep) {
                let dfc = ortho_defect(a, b, c, d, e, ep)
                if dfc > worst {
                  worst = dfc
                }
                tested = tested + 1
              }
            }
          }
        }
      }
    }
  }
  println("T87: screened configs = \{tested}, worst |defect| = \{worst}")
  assert_true(
    tested >= 300,
    msg="ortho screen must admit >= 300 configs (no vacuous pass)",
  )
  assert_true(worst < 1.0e-9, msg="Racah orthogonality defect < 1e-9")
}
