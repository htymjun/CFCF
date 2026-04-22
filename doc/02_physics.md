# 物理モデルと数値手法

## 1. ソッド衝撃波管問題

### 問題設定

ソッド衝撃波管 (Sod shock tube) は圧縮性流れの数値解法の標準的な検証問題である．
長さ $L_x$ の 1 次元管の中央に隔膜があり，左右で圧力・密度が異なる初期状態から始める：

|  | 左側 ($x \leq L_x/2$) | 右側 ($x > L_x/2$) |
|---|---|---|
| 密度 $\rho$ | $\rho_l = 1.293\ \mathrm{kg/m^3}$ | $0.125\,\rho_l$ |
| 速度 $u$ | $0$ | $0$ |
| 圧力 $p$ | $\rho_l R T_{lr}$ | $0.1\,\rho_l R T_{lr}$ |

隔膜を取り除くと，衝撃波・接触不連続・膨張波が伝播する．解析解が既知であるため，数値解法の精度検証に広く用いられる．

### 無次元化

出力ファイル `keep.dat` の値は以下のように無次元化されている：

- 位置: $x / L_x$
- 密度: $\rho / \rho_l$
- 速度: $u / a$（$a = \sqrt{R T_{lr}}$ : 音速）
- 圧力: $p / (\rho_l R T_{lr})$

---

## 2. 支配方程式：1次元圧縮性 Navier-Stokes 方程式

連続式・運動量式・エネルギー式を保存形式で書くと：

$$\frac{\partial \mathbf{q}}{\partial t} + \frac{\partial (\mathbf{e} - \mathbf{e}_v)}{\partial x} = 0$$

ここで保存変数ベクトル $\mathbf{q}$，対流フラックス $\mathbf{e}$，粘性フラックス $\mathbf{e}_v$ はそれぞれ：

$$
\mathbf{q} = \begin{pmatrix} \rho \\ \rho u \\ E \end{pmatrix}, \quad
\mathbf{e} = \begin{pmatrix} \rho u \\ \rho u^2 + p \\ (E + p) u \end{pmatrix}, \quad
\mathbf{e}_v = \begin{pmatrix} 0 \\ \tau \\ \tau u - q_T \end{pmatrix}
$$

$E = \rho e + \frac{1}{2}\rho u^2$ は全エネルギー，$e = c_v T$ は内部エネルギー，$c_v = R/(\gamma-1)$ は定積比熱．

圧力は理想気体の状態方程式 $p = \rho R T$ および

$$p = (\gamma - 1)\!\left(E - \frac{1}{2}\rho u^2\right)$$

から計算できる（$\gamma = 1.4$）．

### 粘性応力と熱流束

粘性応力と熱流束（フーリエ則）：

$$\tau = \frac{4}{3}\mu \frac{\partial u}{\partial x}, \qquad q_T = -\frac{\mu C_p}{\mathrm{Pr}} \frac{\partial T}{\partial x}$$

$\mathrm{Pr} = 0.72$（プラントル数）．

### サザーランドの粘性則

温度依存の動粘性係数：

$$\mu(T) = \mu_0 \left(\frac{T_0 + S}{T + S}\right)\!\left(\frac{T}{T_0}\right)^{3/2}$$

パラメータ: $\mu_0 = 1.716 \times 10^{-5}\ \mathrm{Pa{\cdot}s}$，$T_0 = 273.2\ \mathrm{K}$，$S = 111\ \mathrm{K}$（サザーランド定数）．

---

## 3. 空間離散化：有限差分法

格子点 $j = 1, \ldots, N_x$ を等間隔 $\Delta x = L_x / (N_x - 1)$ で配置する．
フラックスは隣接格子点の**中間点**（界面）$j+\frac{1}{2}$ で定義する：

$$\frac{\partial \mathbf{q}_j}{\partial t} = -\frac{\mathbf{e}_{j+1/2} - \mathbf{e}_{j-1/2}}{\Delta x} + \frac{\mathbf{e}_{v,j+1/2} - \mathbf{e}_{v,j-1/2}}{\Delta x}$$

コード中の配列 `e(3, nx-1)`, `ev(3, nx-1)` の添字 $j$ が界面 $j+\frac{1}{2}$ に対応する．

---

## 4. KEEP スキーム（対流フラックス）

### 動機

通常の中心差分で対流項を離散化すると，計算が進むにつれて運動エネルギーが非物理的に増大し，数値的に発散することがある．
**KEEP (Kinetic Energy and Enthalpy Preserving) スキーム**は，離散レベルで運動エネルギーとエントロピーの保存性を厳密に満たすように対流フラックスを設計した手法であり，安定性が高い．

### フラックスの定義

界面 $j+\frac{1}{2}$ における対流フラックスを以下のように定義する：

**質量フラックス（補助量）**

$$c_{j+1/2} = \frac{1}{4}(\rho_j + \rho_{j+1})(u_j + u_{j+1})$$

**運動量フラックス**

$$e_{\rho u, j+1/2} = \frac{1}{2} c_{j+1/2}(u_j + u_{j+1}) + \frac{1}{2}(p_j + p_{j+1})$$

**エネルギーフラックス**

$$e_{E, j+1/2} = \underbrace{c_{j+1/2} \cdot \frac{R}{2(\gamma-1)}(T_j + T_{j+1})}_{\text{内部エネルギー輸送}} + \underbrace{\frac{1}{2} c_{j+1/2}\, u_j u_{j+1}}_{\text{運動エネルギー輸送}} + \underbrace{\frac{1}{2}(u_j p_{j+1} + u_{j+1} p_j)}_{\text{圧力仕事}}$$

コード中の変数名との対応：

| 数式 | コード変数 |
|---|---|
| $c_{j+1/2}$ | `c` |
| $\frac{1}{2}(p_j + p_{j+1})$ | `g`（圧力勾配項） |
| 内部エネルギー輸送 | `ie` |
| 運動エネルギー輸送 | `ke` |
| 圧力仕事 | `Lp` |

### 粘性フラックスの離散化

粘性応力と熱流束を 2 次精度の中心差分で近似する（界面での粘性係数は算術平均）：

$$\tau_{j+1/2} = \frac{4}{3} \cdot \frac{\mu_j + \mu_{j+1}}{2} \cdot \frac{u_{j+1} - u_j}{\Delta x}$$

$$q_{T,j+1/2} = -\frac{C_p(\mu_j + \mu_{j+1})}{2\,\mathrm{Pr}} \cdot \frac{T_{j+1} - T_j}{\Delta x}$$

エネルギー方程式の粘性項には速度と粘性応力の積（粘性散逸）と熱流束が入る：

$$e_{v,E,j+1/2} = \frac{1}{2}\tau_{j+1/2}(u_j + u_{j+1}) - q_{T,j+1/2}$$

コード変数：`tw`（$\tau$），`twu`（粘性散逸），`kT`（熱流束）．

---

## 5. 時間積分：1次精度オイラー法

時間方向は最もシンプルな前進オイラー法を採用している：

$$\mathbf{q}_j^{n+1} = \mathbf{q}_j^n - \frac{\Delta t}{\Delta x}\left[(\mathbf{e}_{j+1/2}^n - \mathbf{e}_{j-1/2}^n) - (\mathbf{e}_{v,j+1/2}^n - \mathbf{e}_{v,j-1/2}^n)\right]$$

**CFL 条件**により時間刻み幅を決定する：

$$\Delta t = \mathrm{CFL} \cdot \frac{\Delta x}{a}, \quad \mathrm{CFL} = 0.1$$

コードでは `dtdx = CFL * dx / (a * dx) = CFL / a` と定義されており，更新式は

```fortran
q(j) = q(j) - dtdx * (-e(j-1) + e(j) - (-ev(j-1) + ev(j)))
```

と書ける．

---

## 6. 計算パラメータのまとめ

| パラメータ | 記号 | 値 |
|---|---|---|
| 格子点数 | $N_x$ | 4096 |
| 時間ステップ数 | $N_t$ | 2500 |
| レイノルズ数 | $\mathrm{Re}$ | 25000 |
| CFL 数 | $\mathrm{CFL}$ | 0.1 |
| プラントル数 | $\mathrm{Pr}$ | 0.72 |
| 比熱比 | $\gamma$ | 1.4 |
| 気体定数 | $R$ | 287.03 J/(kg·K) |
| 参照温度 | $T_{lr}$ | 300 K |
| 参照密度 | $\rho_l$ | 1.293 kg/m³ |

領域長さ $L_x$ はレイノルズ数から逆算される：

$$L_x = \frac{\mathrm{Re}\cdot\mu_0}{\rho_l \cdot a}$$
