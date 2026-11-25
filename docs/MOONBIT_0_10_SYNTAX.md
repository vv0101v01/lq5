# MoonBit 0.10 最新語法核查筆記

> 核查日期：2026-11-17（Asia/Hong_Kong）。以官方 `latest` 穩定通道實測，不是憑記憶寫語法。

## 0. 版本證據

```
moon 0.1.20261104 (940aa18 2026-11-04) /tmp/mb/moon_home/bin/moon
moonc v0.10.20+65efd5b93 (2026-11-05) /tmp/mb/moon_home/bin/moonc
moonrun 0.1.20261104 (940aa18 2026-11-04) /tmp/mb/moon_home/bin/moonrun
Feature flags enabled: rr_moon_mod,rr_moon_pkg
```

官方 [Updates 索引](https://www.moonbitlang.com/updates/) 與 [2026-11-05 Release](https://www.moonbitlang.com/updates/2026/11/05/index) 標示 `moonc version: v0.10.19`；CLI `latest` 通道提供同日補丁 `v0.10.20+65efd5b93`。兩者視為同一穩定世代。v0.10.19 的語言／工具鏈變更包含：

1. `derive(...)` 從型別定義後方移到**型別上方的屬性** `#derive(Debug, Eq, Hash)`，`moon fmt` 可自動遷移。
2. `Int` / `UInt` 正名為 `Int32` / `UInt32`（別名仍可用）。
3. `compare` 遇 NaN 結果未定義；改用 `Float::total_compare` / `Double::total_compare`。
4. `extenum` 拆分為 `pub extenum E` + `extenum E += { A }`。
5. 新增 byte regex `reb"..."`、頂層 `declare let`、`pkgtype(kind: "foreign_library")` 原生靜態庫。
6. **`moon.pkg.json` / `moon.mod.json` 不再被接受**；`moon fmt` 是遷移入口。
7. `moon test` 預設平行執行測試執行檔，可用 `-j N` / `--no-parallelize` 調整；`moon check` 會自動產生 `pkg.generated.mbti`。
8. 使用 core 套件（`@math`、`@debug`、`@test` 等）必須在 `moon.pkg` 明確 `import`，否則視為錯誤（不是警告）。

來源：[2](https://www.moonbitlang.com/updates/2026/11/05/index)、[1](https://www.moonbitlang.com/updates/)、[stable 文件：套件設定](https://docs.moonbitlang.com/en/stable/toolchain/moon/package.html)、[stable 文件：專案與套件](https://docs.moonbitlang.com/en/stable/language/packages.html)、[stable 文件：模組設定](https://docs.moonbitlang.com/en/stable/toolchain/moon/module.html)。

## 1. 本專案實際使用、且經 `moon check --deny-warn` 驗證的語法

### 1.1 設定檔 DSL（`moon.mod` / `moon.pkg`）

```text
# moon.mod
name = "arena/superemergence"
version = "0.1.0"
source = "src"
preferred_target = "wasm"
```

```text
# src/moon.pkg
import {
  "moonbitlang/core/math",
  "moonbitlang/core/debug",
}
```

```text
# src/cmd/main/moon.pkg
import {
  "arena/superemergence" @geo,
}

pkgtype(kind: "executable")
```

注意：舊的 `options("is-main": true)` 已被 `pkgtype(kind: "executable")` 取代；`source = "src"` 之後，可執行套件也必須放在 `src/` 底下，否則 `moon run cmd/main` 會找不到套件（`moon run src/cmd/main` 才對）。

### 1.2 新 `#derive` 屬性 + `pub extend`（避免隱式方法提升警告）

```moonbit
#derive(Debug, Eq)
pub(all) struct GradedDim {
  m : Int32
  n : Int32
}

pub extend GradedDim with Eq::{equal, not_equal}

pub extend GradedDim with @moonbitlang/core/debug.Debug::{to_repr}
```

若省略 `pub extend`，0.10 會發出 `implicit_impl_as_method` 警告（`--deny-warn` 下會變成錯誤）。`trait` + `pub impl ... with ...` 也一樣：

```moonbit
pub trait PosKernel {
  eval(Self, Int32) -> Double
}

pub impl PosKernel for PolarPoPE with eval(self, delta) { ... }

pub extend PolarPoPE with PosKernel::{eval}
```

### 1.3 數值與位元

- `Int32`（32-bit）、`Double`（IEEE754）、`i.to_double()`、`x.sqrt()`、`x.abs()`、`@math.cos/exp/ln/acosh/pow/atan2`、`@math.PI`。
- 位元運算子：`&`、`|`、`^`、`<<`、`>>`。`popcount` 直接用 `v & 1` + `v >> 1` 實作。
- 浮點字面值**必須有小數點**：`1.0e-12` 合法，`1e-12` 會被解析成 `1` 後接識別字而錯。
- 比較浮點一律用自訂 `approx(a, b, tol~=...)`；`assert_eq` 會要求 `Eq + @debug.Debug`。

### 1.4 有預設值的參數必須是 labelled argument

```moonbit
pub fn approx(a : Double, b : Double, tol~ : Double = 1.0e-9) -> Bool { ... }
approx(1.0, 1.0)              // OK
approx(1.0, 1.0, tol=1e-12)  // OK
```

寫成 `tol : Double = 1.0e-9`（positional + default）會被編譯器拒絕：`Only labelled arguments can have default value.`

### 1.5 迴圈與範圍

```moonbit
for i in 0..<n { ... }          // 半開
for j in 0..=steps { ... }      // 含端點
while i >= 0 { ...; i -= 1 }
```

`for i in k + 1..<n` 有優先權歧義，必須寫 `for i in (k + 1)..<n`。`let mut` 只在變數本身被重新指派時需要；`Array` 是參考型別，`push`／`arr[i] = v`／`arr.sort()` 不需要 `mut`（多寫會觸發 `unused_mut` 警告）。

### 1.6 測試、ADT、Option、labelled args、泛型

```moonbit
test "T8: Super-Jacobi Identity for Mat(2|1)" {
  assert_true(super_jacobiator_m(x, y, z).norm() < 1.0e-11)
  assert_eq(popcount(0b101101), 4)
  assert_true(approx(a, b, tol=1.0e-12), msg="...")  // msg? 是 labelled optional
}

pub(all) enum Term {
  Var(Int32)
  Const(Double)
  Sym(String)
  Add(Term, Term)
  Mul(Term, Term)
}

pub fn[T : PosKernel] kernel_gram(k : T, len : Int32) -> Array[Double] { ... }

fn probe(xs : Array[(Int32, Double)], key : Int32) -> Double? { ... Some(v) ... None }
```

`test "..."` 是 white-box 測試（也可放 `_wbtest.mbt`）；`_test.mbt` 是 black-box，需 `@pkg.` 明確限定。來源：[Writing Tests](https://docs.moonbitlang.com/en/latest/language/tests.html)、[Fundamentals](https://docs.moonbitlang.com/en/stable/language/fundamentals.html)。

### 1.7 指令

| 指令 | 用途 |
| --- | --- |
| `moon new <path>` | 產生新模組（0.10 直接給 `moon.mod` / `moon.pkg`） |
| `moon fmt` / `moon fmt --check` | 格式化／檢查（同時遷移舊 JSON 設定） |
| `moon check --deny-warn` | 型別檢查，警告視為錯誤 |
| `moon test --deny-warn -j 4` | 測試（0.10.19 起預設平行） |
| `moon build --target wasm` | 產出 `.wasm` |
| `moon run src/cmd/main` | 執行 |
| `moon version --all` | 完整版本元組 |

## 2. 實測踩過的坑（已修正，留作備忘）

| 現象 | 原因 | 修正 |
| --- | --- | --- |
| `unexpected token id (lowercase start)` | `1e-12` 沒有小數點 | `1.0e-12` |
| `Cannot find package to build ... cmd/main` | `source = "src"` 時套件搜尋根目錄是 `src/` | 移到 `src/cmd/main` |
| `Only labelled arguments can have default value` | positional 參數給預設值 | 改成 `name~ : T = ...` |
| `Operator + and ..< are ambiguous` | 範圍運算子優先權 | 加括號 `(k + 1)..<n` |
| `unused_package` / `unused_mut` / `unused_value` | `--deny-warn` 嚴格 | 移除未用 import、`mut`、參數改名 `_self` |
| power iteration 收斂到 0 | 初始向量 `1/√n` 正好落在圖拉普拉斯的零特徵空間 | 改用 `sin(0.7(i+1))` 初始向量 |

## 3. 工具鏈下載與完整性

```bash
curl -fsSL https://cli.moonbitlang.com/install/unix.sh | bash
curl -O -L https://cli.moonbitlang.com/binaries/latest/moonbit-linux-x86_64.tar.gz.sha256
curl -O -L https://cli.moonbitlang.com/binaries/latest/moonbit-linux-x86_64.tar.gz
sha256sum -c moonbit-linux-x86_64.tar.gz.sha256
```

來源：[MoonBit Download / Verifying binaries](https://www.moonbitlang.com/download)。本專案實測 `80d987523308fa490eba9f1baf98df996a5a7891e118357553b7c2608c0597c5` → `OK`。

---

## 11. 第二輪重裝：`moonc v0.10.22+124d61d05`（2026-11-19）

第二輪開發時 `latest` 通道已前進到：

```
moon    0.1.20261119 (738607b 2026-11-19)   /tmp/mb/moon_home/bin/moon
moonc   v0.10.22+124d61d05 (2026-11-19)     /tmp/mb/moon_home/bin/moonc
moonrun 0.1.20261119 (738607b 2026-11-19)   /tmp/mb/moon_home/bin/moonrun
```

安裝方式不變（`scripts/install-moonbit.sh` → 官方 `install/unix.sh` + `binaries/latest/*.sha256` 驗證）。
第一輪的 `v0.10.20+65efd5b93` 與第二輪的 `v0.10.22+124d61d05` **都能編譯本專案**
（兩輪之間只用了同一組 `0.10.x` 語法）。

### 11.1 這輪新踩到的坑（全部已被 `--deny-warn` 逼出來）

| 症狀（原文） | 原因 | 正確寫法 |
|---|---|---|
| `Value sqrt not found in package math` | `@math.sqrt` / `@math.abs` 不是函數 | `x.sqrt()`、`x.abs()`（方法形式） |
| `Warning (deprecated): Use PI instead` | `@math.pi()` 已棄用 | `@math.PI` |
| `Parse error, unexpected token *`（跨行算式） | 0.10 不允許在行尾斷開二元運算 | 收成一行，或先算中間變數 |
| `Parse error, unexpected token id`（`1e-12`） | 浮點字面量需有小數點 | `1.0e-12` |
| `Error Warning (unused_mut)` | 只改內容（`a[i] = ...`、`push`）的陣列不需要 `mut` | 去掉 `let mut` |
| `assert_eq!` / `print` → `unbound` | 巨集／函數不存在 | `assert_true(a == b)`、`println` |
| `Type _/0 has no method op_get`（匿名函數參數） | 閉包參數型別無法推導 | `fn(q : Array[Double]) { ... }` |
| `Array<Double>` 解析失敗 | 泛型用方括號 | `Array[Double]` |
| `unexpected token if`（作運算元） | `if` 需放在獨立 `let` | `let d = if c { a } else { b }` |
| `The value identifier ~+ is unbound`（跨行 `+` 串字串） | 只能用行尾 `+` 續行 | `"a" +` 換行 `"b"`，或 `\{expr}` 插值 |

### 11.2 本輪好用、值得記下來的 0.10 寫法

```moonbit
// 字串插值 + 行尾 `+` 續行（EmergenceTrace::summary 的風格）
let s = "agents: \{cfg.agents}, steps: \{cfg.steps}\n" +
  "R(0) = \{first} -> R(T) = \{last}\n"

// 標籤參數附預設值
pub fn project_to_leaf(p0 : Array[Double], p1 : Array[Double], iters~ : Int32 = 3) -> Array[Double]

// 閉包傳遞（Poisson 括號需要「函數作為引數」）
pub fn poisson_bracket(f : (Array[Double]) -> Double, g : (Array[Double]) -> Double, p : Array[Double], h~ : Double = 1.0e-5) -> Double

// 測試內取得 Int32 → Int
assert_true(tr.rows.length() == cfg.steps.to_int() + 1)

// CSV 解析（測試內）
tr.rows[0].split(",").collect().length()
```

---

## 12. 第三輪前沿擴充核查（2026-11-24）：`moonc v0.10.22+124d61d05` 原生硬體與最新標準庫特性

核查日期：2026-11-24（Asia/Hong_Kong），SHA-256（`moonbit-linux-x86_64.tar.gz`）：`aeb78ded1d5e25185ed788328cacef5e843a9ab923145523b7652e8d082f74fb`。

本輪從 `~/.moon/lib/core/CHANGELOG.md` 與官方 Updates 挖掘並實測落地的 **MoonBit 0.10.19–0.10.22 最新語法與標準庫特性**：

1. **硬體位元內建指令（`%i32_popcnt`, `%u64.popcnt`, `%u64.ctz`）**：
   `Int32::popcnt()`、`UInt64::popcnt()`、`UInt64::ctz()` 直接編譯為 Wasm/Native 底層位元指令，讓 `koszul_sign64` 只需迭代非零位元 `m = m & (m - 1UL)`（`17_sparse_super_m_n.mbt`）。
2. **`Double::total_compare`（IEEE-754 全序比較）**：
   消除 0.10.19 起 `compare` 遇 `NaN` 未定義的風險，實作 `sort_total` 並用於 `Hybrid::canonical` 軌道投影（`00_core.mbt`, `09_unify.mbt`）。
3. **`@sorted_map.SortedMap` 有界區間導航 `range_bounds` + `@cmp.Bound` + `for k, v in map`（`Iter2` 模式解構）**：
   0.10.22 新增 `@cmp.Bound::{Included, Excluded, Unbounded}` 與 `SortedMap::range_bounds`，配合 0.10.14+ 的 `for mask, v in self.terms.range_bounds(@cmp.Bound::Excluded(0UL), @cmp.Bound::Unbounded)` 直接跳過 `0UL`（body 項）掃描稀疏靈魂（`17_sparse_super_m_n.mbt`）。
4. **結構展開更新語法 `{ ..base_struct, field: new_val }`**：
   於複數 `SL_2(C)` 矩陣構造 `Mat2Cx`（`19_opers_nonabelian_hodge.mbt`）無痛更新行列式正規化欄位。
5. **非同步標準庫遷移備忘**：`moonbitlang/core/async` 已併入標準庫（需 `allow-experimental-packages = true`）；管道語法已棄用 `x |> fn(y) { ... }`，統一改用 `x |> y => { ... }`。
