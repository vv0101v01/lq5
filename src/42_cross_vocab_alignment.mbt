// =====================================================================
// 42 — 跨詞表（cross-vocab）對齊融合 A73，第九輪 deferred 項目 3
// （HEXAD 計劃 §5.4/§11.1，算法 A40），證書 T92–T94。
//
// 與 src/39 的同詞表 fuse_* 對偶：α 混合共享候選、(1−α) 混合對側
// 兜底。token 以 String 表示（ASCII 測試詞表；定律與表示無關——
// 真實 tokenizer 的非 UTF-8 位元組 token 可原樣換成 Bytes，屬記錄
// 在案的邊界）。前綴相容兜底：|候選| 與 token 任一方為對方前綴
// （equality ∨ starts_with），取最大 logprob；無相容候選時退化到
// min(對側) − margin（floor，不再硬否決 —— 對比 T81 fuse_product）。
// 輸出序 = left 序在前、right 獨有在後（同詞表情形與 fuse_product
// 逐位對應，故 T92 可做精確歸約斷言）。
// =====================================================================

///|
/// token 是否前綴相容（相等或互為前綴）。
fn align_prefix(a : String, b : String) -> Bool {
  a == b || a.has_prefix(b) || b.has_prefix(a)
}

///|
/// token 在表中的下標（無則 −1）。
fn align_index(list : Array[(String, Double)], tok : String) -> Int32 {
  let mut idx = 0 - 1
  for i in 0..<list.length() {
    if list[i].0 == tok {
      idx = i
    }
  }
  idx
}

///|
/// 對側兜底分：前綴相容候選的最大 logprob；一個都沒有則
/// min(對側 logprob) − margin。對側為空時返回 −∞（安全值）。
fn align_fallback(
  other : Array[(String, Double)],
  tok : String,
  margin : Double,
) -> Double {
  if other.length() == 0 {
    return fuse_neg_inf()
  }
  let mut best = fuse_neg_inf()
  let mut found = false
  let mut minv = other[0].1
  for cand in other {
    if cand.1 < minv {
      minv = cand.1
    }
    if align_prefix(cand.0, tok) {
      if !found || cand.1 > best {
        best = cand.1
        found = true
      }
    }
  }
  if found {
    best
  } else {
    minv - margin
  }
}

///|
/// A73（raw）: 對齊打分 ——
///   共享 token：  α·own + (1−α)·對側同 token
///   left 獨有：  α·own + (1−α)·fallback(對側)
///   right 獨有： (1−α)·own + α·fallback(對側)
/// α 永遠配 left 自有項（swap 反變的前提）。輸出序 = left 序 + right 獨有。
pub fn align_scores(
  left : Array[(String, Double)],
  right : Array[(String, Double)],
  alpha : Double,
  margin : Double,
) -> Array[(String, Double)] {
  let om = 1.0 - alpha
  let out : Array[(String, Double)] = []
  for t in left {
    let j = align_index(right, t.0)
    let s = if j >= 0 {
      alpha * t.1 + om * right[j].1
    } else {
      alpha * t.1 + om * align_fallback(right, t.0, margin)
    }
    out.push((t.0, s))
  }
  for t in right {
    if align_index(left, t.0) < 0 {
      let s = om * t.1 + alpha * align_fallback(left, t.0, margin)
      out.push((t.0, s))
    }
  }
  out
}

///|
/// A73（normalized）: align_scores + fuse_log_softmax（無硬否決）。
pub fn align_fuse(
  left : Array[(String, Double)],
  right : Array[(String, Double)],
  alpha : Double,
  margin : Double,
) -> Array[(String, Double)] {
  let raw = align_scores(left, right, alpha, margin)
  let lps = Array::makei(raw.length(), fn(i) { raw[i].1 })
  let fused = fuse_log_softmax(lps)
  Array::makei(raw.length(), fn(i) { (raw[i].0, fused[i]) })
}

// ---------------------------------------------------------------------
// T92: 同詞表歸約 —— align(L,L,α) 與 fuse_product([row,row],[α,1−α])
// 逐位一致（同浮點運算次序 ⇒ 精確）
// ---------------------------------------------------------------------

///|
test "T92: alignment reduces to the shared-vocabulary product" {
  let left : Array[(String, Double)] = [
    ("a", 0.0 - 0.1),
    ("b", 0.0 - 0.5),
    ("c", 0.0 - 1.2),
    ("d", 0.0 - 2.0),
  ]
  let row = [0.0 - 0.1, 0.0 - 0.5, 0.0 - 1.2, 0.0 - 2.0]
  let alpha = 0.3
  let om = 1.0 - alpha
  let got = align_fuse(left, left, alpha, 1.0)
  let want = fuse_product([row, row], [alpha, om])
  assert_true(got.length() == want.length(), msg="same output length")
  let mut i = 0
  while i < got.length() {
    assert_true(got[i].0 == left[i].0, msg="output order = left order")
    assert_true(
      (got[i].1 - want[i]).abs() <= 1.0e-15,
      msg="exact reduction to fuse_product",
    )
    i = i + 1
  }
}

// ---------------------------------------------------------------------
// T93: swap 反變 —— (L,R,α) ↔ (R,L,1−α) 逐 token raw 精確（α=0.5）
// 或 1e-12 內（α=0.3，因 1−(1−α) ≠ α 的浮點 ulp），normalized 1e-12 內
// ---------------------------------------------------------------------

///|
test "T93: swap covariance of the cross-vocabulary alignment" {
  let left : Array[(String, Double)] = [
    ("th", 0.0 - 0.5),
    ("the", 0.0 - 0.6),
    ("xy", 0.0 - 1.0),
  ]
  let right : Array[(String, Double)] = [
    ("th", 0.0 - 0.4),
    ("t", 0.0 - 0.2),
    ("xy", 0.0 - 0.9),
    ("q", 0.0 - 2.0),
  ]
  let union = ["th", "the", "xy", "t", "q"]
  // α = 0.5：own-weight 恆配自有項、加法可交換 ⇒ raw 精確相等
  let r1 = align_scores(left, right, 0.5, 1.0)
  let r2 = align_scores(right, left, 0.5, 1.0)
  assert_true(
    r1.length() == union.length() && r2.length() == union.length(),
    msg="union support both ways",
  )
  for tok in union {
    let i1 = align_index(r1, tok)
    let i2 = align_index(r2, tok)
    assert_true(i1 >= 0 && i2 >= 0, msg="token present in both orders")
    assert_true(r1[i1].1 == r2[i2].1, msg="raw swap-exact at alpha=0.5")
  }
  // α = 0.3：1−(1−0.3) ≠ 0.3 有一個 ulp，raw 差 < 1e-12
  let r3 = align_scores(left, right, 0.3, 1.0)
  let r4 = align_scores(right, left, 0.7, 1.0)
  for tok in union {
    let i3 = align_index(r3, tok)
    let i4 = align_index(r4, tok)
    assert_true(
      (r3[i3].1 - r4[i4].1).abs() <= 1.0e-12,
      msg="raw swap at alpha=0.3 < 1e-12",
    )
  }
  // normalized：logsumexp 輸出序不同 ⇒ 只保證 1e-12
  let f1 = align_fuse(left, right, 0.5, 1.0)
  let f2 = align_fuse(right, left, 0.5, 1.0)
  for tok in union {
    let i1 = align_index(f1, tok)
    let i2 = align_index(f2, tok)
    assert_true(
      (f1[i1].1 - f2[i2].1).abs() <= 1.0e-12,
      msg="normalized swap < 1e-12",
    )
  }
}

// ---------------------------------------------------------------------
// T94: 前綴兜底精確值、floor 精確值、無硬否決（對比 T81）
// ---------------------------------------------------------------------

///|
test "T94: prefix fallback, floor value, and no hard veto" {
  // (a) 完全不相容詞表：雙側都保留、分數有限（fuse_product 會全 −∞）
  let la : Array[(String, Double)] = [("aa", 0.0 - 0.1)]
  let ra : Array[(String, Double)] = [("bb", 0.0 - 0.2)]
  let rawa = align_scores(la, ra, 0.5, 1.0)
  assert_true(rawa.length() == 2, msg="both sides kept: no veto")
  let saa = rawa[align_index(rawa, "aa")].1
  let sbb = rawa[align_index(rawa, "bb")].1
  // floor = min(對側) − margin：aa → 0.5·(−0.1) + 0.5·(−0.2−1.0) = −0.65
  assert_true(
    (saa - (0.0 - 0.65)).abs() <= 1.0e-12,
    msg="aa floor score = -0.65",
  )
  // bb → 0.5·(−0.2) + 0.5·(−0.1−1.0) = −0.65
  assert_true(
    (sbb - (0.0 - 0.65)).abs() <= 1.0e-12,
    msg="bb floor score = -0.65",
  )
  let fused = align_fuse(la, ra, 0.5, 1.0)
  assert_true(fused.length() == 2, msg="fused keeps both")
  for p in fused {
    assert_true(
      !fuse_is_zero(p.1),
      msg="all fused scores finite (T81 would veto)",
    )
    assert_true(p.1 > 0.0 - 10.0, msg="finite, not -inf")
  }
  // (b) 前綴兜底精確值："the" 獨有 → 對側前綴 "th"/"t" 的最大 = −0.2
  let lb : Array[(String, Double)] = [("the", 0.0 - 0.1)]
  let rb : Array[(String, Double)] = [
    ("th", 0.0 - 0.5),
    ("t", 0.0 - 0.2),
    ("x", 0.0 - 0.01),
  ]
  let rawb = align_scores(lb, rb, 0.5, 1.0)
  assert_true(rawb.length() == 4, msg="union size 4")
  let sthe = rawb[align_index(rawb, "the")].1
  // 0.5·(−0.1) + 0.5·max(−0.5, −0.2) = −0.05 − 0.1 = −0.15
  assert_true(
    (sthe - (0.0 - 0.15)).abs() <= 1.0e-12,
    msg="prefix fallback = -0.15",
  )
  let sth = rawb[align_index(rawb, "th")].1
  // "th" 獨有（right）→ 自有 −0.5、兜底 into {"the"} = −0.1：−0.25 − 0.05 = −0.3
  assert_true(
    (sth - (0.0 - 0.3)).abs() <= 1.0e-12,
    msg="right-only prefix = -0.3",
  )
  // (c) floor = min(對側) − margin 精確值
  let lc : Array[(String, Double)] = [("zz", 0.0 - 0.3)]
  let rc : Array[(String, Double)] = [("th", 0.0 - 0.5), ("t", 0.0 - 0.2)]
  let rawc = align_scores(lc, rc, 0.5, 0.4)
  let szz = rawc[align_index(rawc, "zz")].1
  // floor = min(−0.5, −0.2) − 0.4 = −0.9：0.5·(−0.3) + 0.5·(−0.9) = −0.6
  assert_true((szz - (0.0 - 0.6)).abs() <= 1.0e-12, msg="floor = -0.6")
}
