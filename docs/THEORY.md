# 由計算機抽象導出的定理與演算法

> 目標：把「數學名詞」變成「可執行的數學結構」，再從抽象本身（而不是從外部硬套）導出定律、定理與演算法，最後用於類智能體湧現（agent-like emergence）。
>
> 所有定理都附有 MoonBit 0.10 的可執行憑證：`moon test --deny-warn` 共 82 個測試全部通過（wasm + js，0 warnings）。檔案錨點格式為 `src/xx.mbt :: 函式名 / test 名稱`。第四輪（T35–T46 / A29–A37）詳見 [`docs/UNIVERSE.md`](UNIVERSE.md)。

---

## 0. 方法論：五層抽象階梯

```
數學物件  ──(離散化/編碼)──▶  資料型別(載體)
    │                            │
    │                            ▼
    │                     運算 (operations)
    │                            │
    ▼                            ▼
 定律(laws)  ◀──(可判定化)──  數值斷言 (test)
    │
    ▼
 定理 T*  ──▶  演算法 A*  ──▶  湧現角色
```

關鍵在於**選對編碼**：一旦把 Grassmann 單項式編成位元遮罩，`θ_i² = 0` 就從「一條公理」變成「一個位元測試 `(a & b) != 0`」，而 Koszul 符號就變成 `popcount` 的奇偶性。公理變成可計算的事實，定理就變成可驗證的恆等式，演算法則是這些恆等式的構造性證明。

| # | 數學名詞 | 計算機抽象 | MoonBit 載體 | 湧現角色 |
| --- | --- | --- | --- | --- |
| D1 | 龐加萊–黎曼超流形 $\mathcal{M}^{m\mid n}$ | 偶部 = 取樣點 + 度量；奇部 = Λ(n) 纖維 | `SuperManifold`, `SuperNumber` | agent 狀態空間（相位 + 信念 + 概念） |
| D2 | 仿射張量空間 | `(shape, variance, flat data)` | `Tensor` | 座標變換不變量、慣量 |
| D3 | PoPE 位置編碼 | 位置核 = 自伴算子 J 的譜測度 | `trait PosKernel`, `PolarPoPE`, `LegendrePoPE` | 只依賴相對位置的耦合 |
| D4 | 非交換幾何／超矩陣方程 | 分塊矩陣 + 分次標籤 + 有限譜三元組 | `SuperMatrix`, `SpectralTriple` | 距離 = 通訊強度倒數 |
| D5 | 偏微分演化流 | 圖拉普拉斯離散化 Δ + RK2 | `ring_laplacian`, `kuramoto_flow` | 對稱破缺動力學 |
| D6 | 超雅可比恆等式 | 分次結合代數上的超交換子 | `super_jacobiator`, `super_jacobiator_m` | 代數封閉性憑證 |
| D7 | 二階超質量張量軌道空間 | 排序特徵值 = Weyl chamber 投影 | `SuperMassTensor`, `orbit_distance` | 巨觀序參量（含向列序） |
| D8 | 符號–神經合一算子 | 閉包算子 = 正規式 + 排序 | `Hybrid::canonical`, `unification_operator` | 概念共識 |
| D9 | 自然梯度／Fisher 度量 | cometric `F = diag(p) − ppᵀ`，無矩陣 CG | `fisher_apply`, `natural_gradient_cg` | 信念共識 |

---

## 1. 形式定義

**D1（超圖／超函數層）**. 取 $\mathcal{M}^{2\mid 2}$：$\mathcal{M}$ 的偶部是 Poincaré 上半平面 $\mathbb{H}^2=\{(x,y):y>0\}$，$g = y^{-2}I_2$，曲率 $-1$；奇纖維是 $\Lambda(\theta_1,\theta_2)$，其係數向量長度 $2^n$。函數層為 $C^\infty(\mathbb{H}^2)\otimes\Lambda(n)$。超維度 $\operatorname{sdim}=m-n$，超度量 $G=\operatorname{diag}(g_{\text{even}},J_{\text{odd}})$，$J$ 為標準辛形式。

**D2（Grassmann 積）**. 對遮罩 $A,B\subseteq\{0,\dots,n-1\}$：

$$\theta_A\theta_B = \begin{cases}0,& A\cap B\neq\varnothing\\ (-1)^{\#\{(i\in A,j\in B):j<i\}}\,\theta_{A\cup B},& A\cap B=\varnothing\end{cases}$$

$\theta_A:=\theta_{i_1}\cdots\theta_{i_k}$（索引遞增）。Berezin 積分 $\int d\theta_n\cdots d\theta_1\, f := [\,\theta_1\cdots\theta_n\,]f$。

**D3（PoPE 統一）**. 位置編碼是譜三元組 $(\mathcal{A},\mathcal{H},J)$：$\mathcal{A}$ 由自伴算子 $J$ 生成，核 $K=\langle\Phi(\text{pos}),\Phi(\text{pos}')\rangle$ 是 $J$ 的譜測度。

- **Legendre PoPE**：$J$ = 正規化 Legendre 三項遞推的 Jacobi 三對角矩陣，$J_{k,k+1}=k/\sqrt{(2k-1)(2k+1)}$；$\operatorname{spec}(J_N)=\{P_N\text{ 的根}\}$；$K(a,b)=\frac{2}{2a+1}\delta_{ab}$。
- **Polar PoPE**：$J=\operatorname{diag}(\theta_1..\theta_d)$，$\theta_c=\text{base}^{c/d}$；$K(\Delta)=\frac1d\sum_c\cos(\Delta\theta_c+\delta_c)$，分數 $a(t,s)=\sum_c\mu_q\mu_k\cos((s-t)\theta_c+\delta_c)$。

**D4（分次 Sylvester）**. $AX+\sigma XB=C$，$\sigma=(-1)^{|A||X|}$。向量化後 $(I\otimes A+\sigma B^{\mathsf T}\otimes I)\operatorname{vec}X=\operatorname{vec}C$。

**D5（離散 PDE）**. $\dot u=-\kappa Lu$（$L=2I-S-S^{\mathsf T}$）；Kuramoto $\dot\theta_i=\omega_i+\frac{K}{N}\sum_j w(|i-j|)\sin(\theta_j-\theta_i)$。

**D7（超質量張量）**. $G_{\text{even}}=\langle xx^{\mathsf T}\rangle$，$x=(\cos\theta,\sin\theta)$；$G_{\text{odd}}=\langle x_i\wedge x_{i+1}\rangle J$。軌道 = $\{R^{\mathsf T}G_{\text{even}}R:R\in SO(2)\}$，商空間 $\mathbb R^2/\mathfrak S_2$。

**D9（Fisher）**. $g_p(u,v)=\sum_i u_iv_i/p_i$；cometric $F=\operatorname{diag}(p)-pp^{\mathsf T}$；阻尼自然梯度 $x=(F+\lambda I)^{-1}\nabla L$。

---

## 2. 定理（皆附證明綱要與可執行憑證）

### T1 — 纖維維度與冪零性
$\dim\Lambda(n)=2^n$ 且 $\theta_i^2=0$；任意含重複奇生成元的乘積為零。

*證明綱要*：單項式與 $\{0,1\}^n$ 的遮罩一一對應，共 $2^n$ 個。$\theta_i\theta_i$ 對應 $A\cap B=\{i\}\neq\varnothing$，由 D2 得 $0$。在程式裡這等價於 `koszul_sign(a,b) == 0.0` 當 `(a & b) != 0`。
*憑證*：`01_super_algebra.mbt :: koszul_sign / test "T1: Λ(n) is 2^n-dimensional and θ_i^2 = 0"`。

### T2 — 分次交換性與 Koszul 符號
$\theta_A\theta_B=(-1)^{|A||B|}\theta_B\theta_A$（同質元素）。

*證明綱要*：交換 $A,B$ 時逆序對數變為 $|A||B|-\#\{(i\in A,j\in B):j<i\}$（因為 $A\cap B=\varnothing$ 時所有有序對恰分兩類），故符號差 $(-1)^{|A||B|}$。
*憑證*：`test "T2: graded commutativity and the Koszul sign"`（驗證 $\theta_0\theta_2=-\theta_2\theta_0$）。

### T3 — 超 Stokes 定理與超雅可比恆等式
(a) $\int Df=0$ 對任意奇導數 $D=\partial/\partial\theta_i$；(b) $(-1)^{|a||c|}[a,[b,c]]+(-1)^{|b||a|}[b,[c,a]]+(-1)^{|c||b|}[c,[a,b]]=0$；(c) $\operatorname{ad}_a$ 是分次導子：$[a,[b,c]]=[[a,b],c]+(-1)^{|a||b|}[b,[a,c]]$。

*證明綱要*：(a) $\partial_{\theta_i}$ 把遮罩 $A\ni i$ 送到 $A\setminus\{i\}$，永遠得不到全遮罩 $2^n-1$，故 Berezin 頂係數為 0。(b)(c) 在任一**結合**分次代數中展開八項，利用 $|[b,c]|=|b|+|c|$ 逐對抵消。⟹ Chevalley–Eilenberg 微分 $d^2=0$。
*憑證*：`test "T3: Berezin integral is a total derivative"`, `test "T3: super-Jacobi identity and ad as a graded derivation"`（8 種奇偶組合全測），以及矩陣版 `05_supermatrix.mbt :: test "T8: Super-Jacobi Identity for Mat(2|1)"`。

### T4 — 超指數有限終止
對偶元素 $x=b+s$（$s$ 為冪零 soul），$e^x=e^b\sum_{k=0}^{n}s^k/k!$（有限和）。

*證明綱要*：$b$ 與 $s$ 交換；$s$ 中每個單項式至少含一個 $\theta$，故 $s^{n+1}=0$。
*憑證*：`test "T4: exp terminates on nilpotent soul"`（$e^{\theta_0\theta_1}=1+\theta_0\theta_1$）。

### T5 — 超體積與 Berezinian 變數變換
(a) 對區塊對角 $G=\operatorname{diag}(g,J)$，$\operatorname{Ber}(G)=\det g/\det J$；當 $\det J=1$ 時 $\sqrt{\operatorname{Ber}G}=\sqrt{\det g}=y^{-2}$，超體積退化為雙曲體積。
(b) 在 $x\mapsto\lambda x,\ \theta\mapsto\mu\theta$ 下 $\operatorname{Ber}G\mapsto\lambda^{2m}\mu^{-2n}\operatorname{Ber}G$；體積不變 $\iff \mu=\lambda^{m/n}$。
(c) 奇變數線性替換 $\theta\mapsto B\theta$ 給出 $\int d\theta\, f(B\theta)=\det B\int d\theta\, f(\theta)$。

*證明綱要*：(a) Schur 互補 $\operatorname{Ber}M=\det(A-BD^{-1}C)/\det D$，$B=C=0$ 即得。(b) 協變二張量在座標縮放下 $\det$ 取 $\lambda^{2m}$ 與 $\mu^{2n}$。(c) 取 $f=\theta_1\cdots\theta_n$，$\prod_i(\sum_j B_{ij}\theta_j)$ 的頂係數正是 $\det B$（外代數的定義）。
*憑證*：`02_supermanifold.mbt :: berezin_linear_substitution / ber_scaling_law / test "T5: Berezinian = det ratio..."`、`05_supermatrix.mbt :: SuperMatrix::berezinian / test "T8: supertrace cyclicity..."`。

### T6 — 仿射不變量與平坦性
(a) 縮併 $\operatorname{tr}$ 在仿射座標變換下不變；(b) 常數 $A$ 的仿射映射 $x\mapsto Ax+b$ 使 $\Gamma'^{\,k}_{ij}=(A^{-1})^k_lA^p_iA^q_j\Gamma^l_{pq}$，故 $\Gamma=0\Rightarrow\Gamma'=0\Rightarrow R=0$；(c) 對稱 $\Gamma$ 無撓；$R^l_{ijk}=\Gamma^l_{ip}\Gamma^p_{jk}-\Gamma^l_{jp}\Gamma^p_{ik}$。

*證明綱要*：$\partial A=0$ 使聯絡變換式只剩齊次項；曲率是 $\Gamma$ 的二次式。旋轉 $R\in SO(2)$ 滿足 $R^{\mathsf T}IR=I$，是等距同構。
*憑證*：`03_affine_tensor.mbt :: transformed_connection / curvature_norm / test "T6: affine isometries and flatness are affine invariants"`。

### T7 — 兩個 PoPE 的譜統一與解耦合
(a) Gauss–Legendre 求積下 $\langle P_a,P_b\rangle=\frac{2}{2a+1}\delta_{ab}$；Jacobi 矩陣 $J_N$ 的特徵值恰為 $P_N$ 的根。
(b) Polar PoPE 核為 Toeplitz（平移等變）、正定；分數可分解 $a(t,s)=\langle \mu_q\odot\mu_k,\ \cos((s-t)\theta+\delta)\rangle$，不含內容相位。RoPE 的二維分數為 $q^{\mathsf T}R(\Delta\theta)k=\mu_q\mu_k\cos(\Delta\theta+\phi_k-\phi_q)$，含 what–where 交叉項。
(c) 兩者都是譜三元組：相對位置由線性算子（三項遞推／相位旋轉）復原，不需絕對座標。

*證明綱要*：(a) Favard 定理：正交多項式必有 $xP_n=a_nP_{n+1}+b_nP_n+a_{n-1}P_{n-1}$；截斷 Jacobi 矩陣的特徵多項式即 $P_N$。(b) $\cos$ 是 U(1) 特徵，混合仍正定；$K(t,s)=K(t-s)$。(c) 兩者的核都是某個自伴 $J$ 的譜測度。
*憑證*：`04_pope.mbt :: gauss_legendre / legendre_jacobi / char_poly_at / pope_score / rope_2d_score / test "T7(L): Favard recurrence..." / test "T7(P): translation equivariance, PSD Gram, what/where decoupling"`。

### T8 — 超矩陣代數
(a) $\operatorname{str}(AB)=(-1)^{|A||B|}\operatorname{str}(BA)$；(b) $\operatorname{Mat}(m|n)$ 配上 $[A,B]=AB-(-1)^{|A||B|}BA$ 是 Lie 超代數（超雅可比）；(c) 兩個奇矩陣的括號是**反交換子** $\{Y,Z\}=YZ+ZY$，結果為偶。

*證明綱要*：超跡是分次循環跡；結合性 ⇒ T3；(c) $(-1)^{1\cdot1}=-1$。
*憑證*：`05_supermatrix.mbt :: test "T8: Super-Jacobi Identity for Mat(2|1)"`；模擬中每個時間步的 `closure_certificate` 也驗證此式（最大缺陷 1.6e-15）。

### T9 — Connes 距離 = 耦合倒數（幾何湧現）
對兩點譜三元組 $A=\mathbb C^2,\ D=\begin{pmatrix}0&t\\ t&0\end{pmatrix}$，$\|[D,f]\|=|f_0-f_1||t|$，故
$$d_{\text{Connes}}(0,1)=\sup\{|f_0-f_1|:\|[D,f]\|\le1\}=1/|t|.$$

*證明綱要*：Lipschitz 半範數直接可算；線性規劃的一維上確界在約束飽和時取得。
*推論*：agent 間通訊耦合 $K=|t|$ 越強，非交換幾何距離越小 ⇒ **度量由交互湧現，而非預設**。
*憑證*：`05_supermatrix.mbt :: SpectralTriple::connes_distance_grid / test "T9: Connes distance = inverse coupling"`。

### T10 — 分次 Sylvester 方程可解性
$AX+\sigma XB=C$ 可解 $\iff \lambda_i(A)+\sigma\lambda_j(B)\neq0\ \forall i,j$；等價於向量化係數矩陣 $K=I\otimes A+\sigma B^{\mathsf T}\otimes I$ 非奇異。

*證明綱要*：$\operatorname{spec}(K)=\{\lambda_i(A)+\sigma\lambda_j(B)\}$（Kronecker 和）；高斯消去即可構造解。
*憑證*：`05_supermatrix.mbt :: graded_sylvester_solve / graded_sylvester_residual / test "T10: graded Sylvester equation..."`，含 $\sigma=+1$ 與 $\sigma=-1$ 兩例。

### T11 — 演化流的守恆、耗散與 CFL
(a) $\sum_i u_i$ 守恆；(b) $\frac{d}{dt}E=-\kappa u^{\mathsf T}Lu\le0$，$E=\frac12\|u\|^2$；(c) 顯式 Euler 穩定 $\iff dt\le 2/(\kappa\lambda_{\max}(L))$，其中環圖 $\lambda_k=2-2\cos(2\pi k/N)$；(d) 共識速率由 $\lambda_2$（algebraic connectivity）控制。

*證明綱要*：$L\mathbf 1=0$ ⇒ 質量守恆；$L\succeq0$ ⇒ 耗散；更新矩陣 $I-\kappa\,dt\,L$ 的譜半徑 $\le1$ 給出 CFL。
*憑證*：`06_flow.mbt :: heat_flow / heat_cfl / algebraic_connectivity / test "T11: heat flow conserves mass..." / test "T11: algebraic connectivity controls convergence rate"`。

### T12 — 軌道不變量與向列序
(a) $SO(2)$ 共軛軌道由 $(\lambda_1,\lambda_2)$（等價地 $(\operatorname{tr},\det)$）完全分類，商空間 $\cong\mathbb R^2/\mathfrak S_2$，商度量 = 排序後特徵值的最佳匹配距離；
(b) 各向異性 $\alpha=\frac{\lambda_1-\lambda_2}{\lambda_1+\lambda_2}=|\langle e^{2i\theta}\rangle|$，即**第二諧波**（向列／反相簇）序參量；一階 Kuramoto 序參量 $R=|\langle e^{i\theta}\rangle|$ 看不到它；
(c) $\operatorname{Ber}G=\det G_{\text{even}}/\det G_{\text{odd}}=\det G_{\text{even}}/w^2$，$w=\langle\sin(\theta_{i+1}-\theta_i)\rangle$ 為手徵流。

*證明綱要*：令 $A=\langle\cos2\theta\rangle,B=\langle\sin2\theta\rangle$，則 $G_{\text{even}}=\frac12\begin{pmatrix}1+A&B\\B&1-A\end{pmatrix}$，$\det=\frac{1-(A^2+B^2)}4$，判別式 $=\frac{\sqrt{A^2+B^2}}2$，故 $\alpha=\sqrt{A^2+B^2}$。反相雙簇（0 與 π）給 $R=0,\alpha=1$。
*憑證*：`07_orbit.mbt :: eigen2 / orbit_distance / anisotropy / nematic_order / test "T12: anisotropy = |<exp(2iθ)>| detects nematic order invisible to R"`。

### T13 — Fisher 度量、KL Hessian 與自然梯度
(a) $F=\operatorname{diag}(p)-pp^{\mathsf T}\succeq0$，$\operatorname{rank}F=K-1$，且 $\mathrm{KL}(p+\varepsilon u\|p)+\mathrm{KL}(p-\varepsilon u\|p)=\varepsilon^2 g_p(u,u)+O(\varepsilon^4)$；
(b) 平均適應度 $L(p)=\langle p,f\rangle$ 的未阻尼自然梯度是 $Ff=p\odot(f-\langle p,f\rangle)$，即**replicator 方程**；
(c) 自然梯度在重參數化（此處為類別置換）下等變：$\mathcal N_{\pi p}(\pi g)=\pi\mathcal N_p(g)$。

*證明綱要*：(a) 對 $\sum u_i=0$ 有 $u^{\mathsf T}Fu=\sum p_iu_i^2>0$，而 $F\mathbf1=0$；KL 二階展開。(b) $g_p$ 的逆在餘切空間正是 $F$。(c) $F_{\pi p}=\Pi F_p\Pi^{\mathsf T}$，阻尼 CG 與置換可交換。
*憑證*：`08_fisher.mbt :: fisher_matrix / natural_gradient_cg / replicator_step / test "T13(a/b/c)"`。

### T14 — 符號–神經合一算子是閉包算子
令 $\operatorname{can}(h)=(\text{正規式}(t),\ \operatorname{sort}(v),\ \operatorname{dedup}(\sigma))$，$U(a,b)=\operatorname{can}(\text{join}(t_a,t_b),\ (1-\eta)v_a+\eta v_b+\gamma(\operatorname{emb}(t_a)-\operatorname{emb}(t_b)),\ \sigma_a\uplus\sigma_b)$。

(a) 重寫系統終止且合流（多項式解釋：$\deg$ 遞減），排序冪等 ⇒ $\operatorname{can}^2=\operatorname{can}$，$U$ 是閉包算子，$\operatorname{Fix}(U)$ 是完備格（「湧現概念」）。
(b) $U(x,x)=x$ 對 $x\in\operatorname{Fix}(U)$：偶部（向量混合）對稱，奇部（符號殘差）反對稱 ⇒ 自反時消失。
(c) 在連通圖上對 $U$ 做混沌迭代（chaotic iteration）收斂到共同不動點，收斂率由 $\lambda_2(L)$ 主導（與 T11 同一譜間隙）。

*證明綱要*：(a) Tarski 不動點定理 + 冪等投影；(c) 迭代算子是雙隨機矩陣的凸組合，Perron–Frobenius 給出 $\le(1-c\lambda_2)^t$。
*憑證*：`09_unify.mbt :: rewrite_fix / unify / unification_operator / consensus_unification / test "T14(i,ii)/(iii)"`；實測概念離散度 0.1819 → 0.0026。

### T15 — 湧現主定理（整合）
在連通環圖、正的平移等變 PoPE 權重 $w$、近似相同內在頻率、Fisher 自然梯度共識與 $U$ 混沌迭代下，耦合流
$$(\dot\theta,\dot p,\dot h)=(\text{Kuramoto}_w,\ F(p)\nabla L,\ U\text{-consensus})$$
滿足：

1. $R(t)\to1$（$U(1)$ 相位對稱自發破缺）；
2. 信念離散度 $\to0$、$\mathrm{KL}\to0$、概念離散度 $\to0$；
3. 對所有 $t$，超雅可比缺陷 $=0$（代數封閉性不因湧現而破壞）；
4. 巨觀序完全由二階超質量張量的**軌道不變量** $(\lambda_1,\lambda_2,\alpha,\operatorname{Ber})$ 捕捉；可觀測的只有商空間 $\mathbb R^2/\mathfrak S_2$，這就是「湧現」的自由度；
5. 有效度量半徑 $=1/K$（T9），擴散時間尺度由 CFL 與 $\lambda_2$ 決定（T11）。

*證明綱要*：Kuramoto 在 $w>0,\ \omega_i\approx0$ 時是勢函數 $-\frac K{2N}\sum_{ij}w_{ij}\cos(\theta_j-\theta_i)$ 的梯度流，全同步為其唯一穩定極小（$R=1$）；信念與概念通道各由 T11/T13 與 T14 收斂；三個通道只透過標量 $R$、均值耦合，故可分離收斂。T3/T8 是代數恆等式，與狀態無關。
*憑證*：`10_emergence.mbt :: run_emergence / test "T15: coupled flow produces emergent order and closed algebra"`；`emergence_trace.csv`。

---

## 3. 演算法（皆由上述定理構造性導出）

| ID | 演算法 | 依據 | 構想 | 複雜度 | 程式錨點 |
| --- | --- | --- | --- | --- | --- |
| A1 | 位元遮罩 Grassmann 乘法 | T1/T2 | `popcount` 算 Koszul 符號，`&` 測重複 | $O(4^n)$ | `SuperNumber::mul` |
| A2 | 左奇導數 + Berezin 積分 | T3/T4 | 取頂遮罩係數；導數清掉一個 bit | $O(2^n)$ | `odd_derivative`, `berezin` |
| A3 | Gauss–Legendre 節點／權重（Newton + Favard） | T7 | $x_{k+1}=x_k-P_N(x_k)/P_N'(x_k)$，$w_k=2/((1-x_k^2)P_N'^2)$ | $O(N^2)$ | `gauss_legendre` |
| A4 | Jacobi 譜驗證 | T7 | $\det(J_N-xI)$ 在節點上為零 | $O(N^4)$ | `legendre_jacobi`, `char_poly_at` |
| A5 | PoPE 相對位置核 + PSD 檢驗 | T7 | Toeplitz Gram；位移冪迭代求 $\lambda_{\min}$ | $O(dN^2)+O(N^2 I)$ | `kernel_gram`, `gram_min_eigenvalue` |
| A6 | Berezinian（Schur 互補）| T5/T8 | $\det(A-BD^{-1}C)/\det D$，含非方塊乘積 | $O((m+n)^3)$ | `SuperMatrix::berezinian` |
| A7 | 分次 Sylvester 求解 | T10 | 向量化 $(I\otimes A+\sigma B^{\mathsf T}\otimes I)x=c$ + 高斯消去 | $O((rc)^3)$ | `graded_sylvester_solve` |
| A8 | CFL 自適應熱流／Heun Kuramoto | T11 | $dt=0.9\cdot2/(\kappa\lambda_{\max})$；RK2 | $O(N^2T)$ | `heat_flow`, `kuramoto_step` |
| A9 | 阻尼 CG 自然梯度 | T13 | $(F+\lambda I)x=g$，僅用 matrix-free $Fv=p\odot(v-\langle p,v\rangle)$ | $O(K\cdot I)$ | `natural_gradient_cg` |
| A10 | 合一閉包 + 圖上混沌迭代 | T14 | 重寫至正規式、合一（含 occurs check）、排序投影、鄰域平均 | $O(R\cdot|\text{Term}| + Nd)$ | `unification_operator`, `consensus_unification` |
| A11 | 軌道商度量 | T12 | 排序特徵值後最佳匹配 | $O(1)$（2×2） | `orbit_distance` |
| A12 | 湧現觀測管線 | T15 | 每步輸出 8 個巨觀量至 CSV | $O(T(N^2+NKI))$ | `run_emergence`, `EmergenceTrace::csv` |

### 代表性偽代碼

**A3 + A4（由三項遞推導出的求積與譜驗證）**

```
c_0 = [1]; c_1 = [0,1]
for k = 1..N-1:  c_{k+1} = ((2k+1)/ (k+1))·x·c_k - (k/(k+1))·c_{k-1}
for i = 0..N-1:
    x = cos(π(i+3/4)/(N+1/2))
    repeat 64: x -= P_N(x)/P_N'(x)
    w_i = 2 / ((1-x²) P_N'(x)²)
J_N[k,k+1] = (k+1)/sqrt((2k+1)(2k+3))
assert det(J_N - x_i I) ≈ 0
```

**A7（分次 Sylvester）**

```
σ = (-1)^{|A||X|};  n = rows(A)·cols(B)
for i,j:  K[i*cc+j, h*cc+j] += A[i,h]
          K[i*cc+j, i*cc+l] += σ·B[l,j]
x = gauss_solve(K, vec(C))
```

**A9（無矩陣自然梯度）**

```
F(v) = p ⊙ (v - ⟨p,v⟩) + λv
CG: r=g, d=r, x=0
    α = ⟨r,r⟩/⟨d,F(d)⟩; x+=αd; r-=αF(d); d=r+βd
```

---

## 4. 實驗結果（32 agents / 240 步 / seed 20261117）

| 觀測量 | $t=0$ | $t=T$ | 說明 |
| --- | --- | --- | --- |
| $R=|\langle e^{i\theta}\rangle|$ | 0.18998 | **0.999956** | 極性同步（T15-1） |
| $\alpha=|\langle e^{2i\theta}\rangle|$ | 0.21297 | **0.999825** | 向列序（T12） |
| 信念離散度 | 0.28917 | **3.11e-9** | Fisher 共識（T13） |
| 概念離散度 | 0.18187 | **0.00261** | U 閉包共識（T14） |
| mean KL to population | 0.05957 | ≈ 0 | 資訊幾何收斂 |
| Super-Jacobi 缺陷（全程最大） | — | **1.55e-15** | 代數封閉（T3/T8） |
| Connes 距離 $=1/K$ | — | 0.3125 | 幾何湧現（T9） |
| $\lambda_2(L)$ / CFL $dt$ | — | 0.03843 / 0.71429 | T11 |
| 軌道特徵值 $[\lambda_1,\lambda_2]$ | — | [0.999912, 0.0000876] | T12 |
| $\operatorname{Ber}G$（終態） | — | 8756.82 | 手徵流趨零 ⇒ 超體積發散 |

完整 241 列軌跡：`../emergence_trace.csv`；視覺化：`report.html`。

---

## 5. 驗證清單

```bash
moon version --all                 # v0.10.20+65efd5b93
moon fmt --check                   # OK
moon check --deny-warn             # OK
moon test --deny-warn              # Total tests: 35, passed: 35, failed: 0
moon build --target wasm           # main.wasm
moon run src/cmd/main              # summary + CSV
```

定理 ↔ 測試對照：T1/T2/T3/T4 → `01`; T5 → `02`; T6 → `03`; T7 → `04`; T8/T9/T10 → `05`; T11 → `06`; T12 → `07`; T13 → `08`; T14 → `09`; T15 → `10`。

---

## 6. 假設、邊界與誠實聲明

1. **離散化**：連續流形以取樣點 + 圖拉普拉斯近似；定理中關於「連續 PDE」的陳述在此為離散類比（T11 的 CFL 是精確的離散結果）。
2. **編碼限制**：Grassmann 纖維用 `Int32` 位元遮罩，建議 $n\le16$；`2^n` 記憶體。
3. **具體化**：Poincaré–Riemann 偶部取 $\mathbb H^2$（$m=2$，曲率 $-1$）；辛奇部取 $n=2$。推廣到 $m|n$ 只需更換度量與區塊維度，`GradedDim` 已備。
4. **Berezin 慣例**：$\int d\theta_n\cdots d\theta_1\,\theta_1\cdots\theta_n=1$；奇線性替換權重為 $\det B$（與偶部的 $\det A^{-1}$ 相反），合起來即 $\operatorname{Ber}^{-1}$。
5. **數值**：所有等式以 `approx` 相對容差判定；0.10.19 起 `compare` 遇 NaN 未定義，故不把 NaN 當鍵。位移冪迭代的初始向量不能落在 $L$ 的零空間（見語法筆記§2）。
6. **「神經」的含義**：`embed` 是決定性特徵映射（代數同態），不是訓練過的網路；此處關注的是符號與向量兩個通道的**代數耦合結構**，可替換成任意可微編碼器而不改定理。
7. **PoPE 同名**：兩個不同構造都屬實，見 `README.md` §6；本專案以「譜三元組」統一，非宣稱兩篇論文等價。

---

## 7. 參考來源

- MoonBit 官方 Updates 索引：<https://www.moonbitlang.com/updates/>
- MoonBit v0.10.19 Release（2026-11-05）：<https://www.moonbitlang.com/updates/2026/11/05/index>
- MoonBit Download / Verifying binaries：<https://www.moonbitlang.com/download>
- MoonBit stable docs — Fundamentals / Packages / Module / Package config / Tests：<https://docs.moonbitlang.com/en/stable/>
- A. Aggarwal, *PoPE: Legendre Orthogonal Polynomials Based Position Encoding for LLMs*, arXiv:2405.04585：<https://arxiv.org/abs/2405.04585>
- A. Gopalakrishnan et al., *Decoupling the "What" and "Where" With Polar Coordinate Positional Embeddings*, arXiv:2509.10534 (ICML 2026)：<https://arxiv.org/abs/2509.10534>
- 數學背景：Berezin / Leites（超流形）、Connes（非交換幾何譜三元組）、Favard（正交多項式遞推）、Amari（資訊幾何）、Kuramoto（同步）、Strogatz（耦合振子）。
