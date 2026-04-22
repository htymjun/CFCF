# コード詳細解説と3バリアントの比較

## 1. 全体構成

サンプルコードは 3 つのバリアントに分かれており，それぞれ同じ物理問題を異なる GPU プログラミングスタイルで解く：

| バリアント | ディレクトリ | GPU の使い方 | メモリレイアウト |
|---|---|---|---|
| AoS ディレクティブ | `keep_directive_AoS/` | `!$cuf kernel do` | AoS: `q(3, nx)` |
| SoA ディレクティブ | `keep_directive_SoA/` | `!$cuf kernel do` | SoA: `q(nx, 3)` |
| 明示的カーネル | `keep_kernel/` | `attributes(global)` | SoA: `q(nx, 3)` |

各バリアントは `parameters.f90`，`solver.f90`，`main.f90` の 3 ファイルから成る．

---

## 2. `parameters.f90`

物理定数・格子パラメータをすべてここで定義する．3 バリアントで共通の内容：

```fortran
module parameters
  implicit none
  integer,  parameter :: dp   = kind(1.0d0)  ! 倍精度の kind 値
  integer,  parameter :: nx   = 4096         ! 格子点数
  integer,  parameter :: nt   = 2500         ! 時間ステップ数
  real(dp), parameter :: gamma = 1.4d0
  real(dp), parameter :: Rgas  = 287.03d0    ! 気体定数 [J/(kg·K)]
  real(dp), parameter :: Pr    = 0.72d0      ! プラントル数
  real(dp), parameter :: mu0   = 1.716d-5   ! 参照粘性 [Pa·s]
  real(dp), parameter :: T0    = 273.2d0    ! 参照温度 [K]
  real(dp), parameter :: S     = 111d0      ! サザーランド定数 [K]
  real(dp), parameter :: Re    = 25000d0    ! レイノルズ数
  real(dp), parameter :: Tlr   = 300.d0     ! 初期温度 [K]
  real(dp), parameter :: CFL   = 0.1d0      ! CFL 数
  real(dp), parameter :: rhol  = 1.293d0    ! 左側初期密度 [kg/m^3]
  ! 以下は上記から計算される派生量
  real(dp), parameter :: Cp   = gamma * Rgas / (gamma - 1.d0)
  real(dp), parameter :: a    = sqrt(Rgas * 300d0)     ! 音速 [m/s]
  real(dp), parameter :: Lx   = Re * mu0 / (1.293d0 * a)  ! 領域長さ [m]
  real(dp), parameter :: dx   = Lx / dble(nx-1)        ! 格子幅 [m]
  real(dp), parameter :: dtdx = CFL * dx / (a * dx)    ! = CFL/a
end module parameters
```

**ポイント:** `dtdx = CFL / a` は時間更新 `q = q - dtdx * (flux_right - flux_left)` で使う $\Delta t / \Delta x$．
`CFL * dx / (a * dx)` の `dx` が約分されていることに注意．

---

## 3. メモリレイアウト：AoS vs SoA

多変数データをどの順序で配列に格納するかはGPUのメモリアクセス効率に影響する．

### AoS (Array of Structures)

変数ごとにインデックスが先に来る：`q(変数番号, 格子点番号)`

```
q(1, 1)  q(2, 1)  q(3, 1)   q(1, 2)  q(2, 2)  q(3, 2)  ...
 ←── 格子点 1 の 3 変数 ──→   ←── 格子点 2 の 3 変数 ──→
```

隣接格子点 `j` と `j+1` の同じ変数（例: `q(1,j)` と `q(1,j+1)`）はメモリ上で 3 要素離れている．

### SoA (Structure of Arrays)

格子点番号が先に来る：`q(格子点番号, 変数番号)`

```
q(1, 1)  q(2, 1)  q(3, 1)  ...  q(nx, 1)   q(1, 2)  ...
 ←──────────── 変数 1 (ρ) の全格子点 ────────→
```

隣接格子点 `j` と `j+1` の同じ変数（例: `q(j,1)` と `q(j+1,1)`）はメモリ上で隣接している．

### GPU での有利不利

GPU では，同一 Warp（32 スレッドのグループ）が連続するメモリ番地にアクセスすると **Coalesced アクセス** となり，メモリ帯域を最大限に利用できる．

本コードのカーネルでは，スレッド $j$ が格子点 $j$ を担当し，Warp 内で $j, j+1, j+2, \ldots$ を並列に処理する．
したがって SoA (`q(j,变量)`) の方が同じ変数へのアクセスが Coalesced となり，メモリ効率が良い．

---

## 4. `solver.f90` の解説（AoS バリアントを基準に）

### `init` サブルーチン（初期条件）

```fortran
subroutine init(nx, q)
  integer, intent(in)     :: nx
  real(dp), intent(inout) :: q(3,nx)   ! ← ホストメモリ上の配列（device 属性なし）
  real(dp) rho, u, p, en
  integer :: j
  do j = 1, nx
    if (j <= nx/2) then
      rho = rhol;  p = rhol * Rgas * Tlr  ! 左側: 高圧
    else
      rho = 0.125d0 * rhol;  p = 0.1d0 * rhol * Rgas * Tlr  ! 右側: 低圧
    end if
    u  = 0d0
    en = p / (gamma - 1.d0)   ! 内部エネルギー = p/(γ-1)（u=0 なので全エネルギー = 内部エネルギー）
    q(1:3, j) = (/ rho, rho * u, en /)
  enddo
end subroutine init
```

初期条件は CPU 上で設定し，その後 `main.f90` で GPU に転送する．

---

### `eev` サブルーチン（フラックス計算）

CUDA Fortran での GPU カーネルの書き方がバリアントによって異なる．

#### AoS バリアント（`keep_directive_AoS/solver.f90`）

```fortran
subroutine eev(nx, q, e, ev)
  integer, intent(in)             :: nx
  real(dp), intent(inout), device :: q(3,nx)       ! AoS レイアウト
  real(dp), intent(inout), device :: e(3,nx-1), ev(3,nx-1)
  real(dp), device :: qq(3,nx), T(nx), mu(nx)      ! GPU 上の作業配列
  real(dp) c, m, g, ie, ke, Lp, tw, twu, kT
  integer :: j

  ! ステップ 1: 保存変数 → 原始変数（密度, 速度, 圧力, 温度, 粘性）
  !$cuf kernel do(1)<<<32,4>>>
  do j = 1, nx
    qq(1,j) = q(1,j)                             ! ρ
    qq(2,j) = q(2,j) / qq(1,j)                  ! u = ρu/ρ
    qq(3,j) = (gamma-1)*(q(3,j) - 0.5*qq(1,j)*qq(2,j)**2)  ! p = (γ-1)(E - ½ρu²)
    T(j)   = qq(3,j) / (qq(1,j) * Rgas)         ! T = p/(ρR)
    mu(j)  = mu0 * ((T0+S)/(T(j)+S)) * (T(j)/T0)**1.5d0    ! サザーランド則
  enddo

  ! ステップ 2: KEEP スキームと粘性フラックス（界面 j+1/2 を計算）
  !$cuf kernel do(1)<<<32,4>>>
  do j = 1, nx-1
    c  = 0.25d0 * (qq(1,j)+qq(1,j+1)) * (qq(2,j)+qq(2,j+1))   ! 質量フラックス補助量
    m  = 0.5d0 * c * (qq(2,j)+qq(2,j+1))                        ! 運動量（対流部）
    g  = 0.5d0 * (qq(3,j)+qq(3,j+1))                            ! 圧力項
    ie = c * 0.5d0 * Rgas * (T(j)+T(j+1)) / (gamma-1.d0)        ! 内部エネルギー輸送
    ke = 0.5d0 * c * qq(2,j) * qq(2,j+1)                        ! 運動エネルギー輸送
    Lp = 0.5d0 * (qq(2,j)*qq(3,j+1) + qq(2,j+1)*qq(3,j))       ! 圧力仕事

    tw  = 4d0*(mu(j)+mu(j+1))*(-qq(2,j)+qq(2,j+1)) / (6d0*dx)  ! 粘性応力
    twu = 0.5d0 * tw * (qq(2,j)+qq(2,j+1))                      ! 粘性散逸
    kT  = Cp*(mu(j)+mu(j+1))*(-T(j)+T(j+1)) / (Pr*2d0*dx)      ! 熱流束

    e(1,j) = c;      e(2,j) = m + g;    e(3,j) = ie + ke + Lp   ! 対流フラックス
    ev(1,j) = 0.d0;  ev(2,j) = tw;      ev(3,j) = twu + kT      ! 粘性フラックス
  enddo
end subroutine eev
```

`qq` は `eev` 内のローカルな GPU 配列．`device` 属性をつけているのでグローバルメモリ上に確保される．

#### SoA バリアント（`keep_directive_SoA/solver.f90`）との違い

インデックスの順序だけが異なる：

```fortran
! AoS: q(変数番号, j)
qq(1,j) = q(1,j)
qq(2,j) = q(2,j) / qq(1,j)

! SoA: q(j, 変数番号)
rho(j) = q(j,1)
u(j)   = q(j,2) / rho(j)
```

SoA では原始変数も別々の 1D 配列 `rho(nx)`, `u(nx)`, `p(nx)`, `T(nx)`, `mu(nx)` に格納することで，隣接スレッドが同一変数配列の連続アドレスにアクセスできる．

#### 明示的カーネルバリアント（`keep_kernel/solver.f90`）

`eev` の2段目（フラックス計算）を明示的なカーネル関数 `kernel` に切り出している：

```fortran
! eev 内での呼び出し
!$cuf kernel do(1)<<<32,4>>>
do j = 1, nx
  rho(j) = q(j,1); ...     ! ステップ 1 はディレクティブのまま
enddo
call kernel<<<nx/128, 128>>>(nx, rho, u, p, T, mu, e, ev)  ! ステップ 2 を明示的カーネルで
```

カーネル本体：

```fortran
attributes(global) subroutine kernel(nx, rho, u, p, T, mu, e, ev)
  integer, intent(in), value    :: nx         ! スカラー引数は value 属性必須
  real(dp), intent(in),  device :: rho(nx), u(nx), p(nx), T(nx), mu(nx)
  real(dp), intent(out), device :: e(nx-1,3), ev(nx-1,3)
  real(dp) c, m, g, ie, ke, Lp, tw, twu, kT
  integer :: j, jt

  jt = threadIdx%x
  j  = (blockIdx%x - 1) * blockDim%x + jt    ! グローバルスレッドインデックス
  if (j > nx-1) return                         ! 境界チェック

  ! フラックス計算（SoA バリアントの内部ループと同一）
  c  = 0.25d0 * (rho(j)+rho(j+1)) * (u(j)+u(j+1))
  ...
  e(j,1) = c;  e(j,2) = m + g;  e(j,3) = ie + ke + Lp
  ev(j,1) = 0.d0;  ev(j,2) = tw;  ev(j,3) = twu + kT
end subroutine kernel
```

`<<<nx/128, 128>>>` で起動すると `nx/128 = 32` ブロック × 128 スレッド = 4096 スレッドが，4095 個の界面を 1 スレッド 1 界面で処理する．

---

### `time_step` サブルーチン（時間積分）

3 バリアントでほぼ同一の実装：

```fortran
subroutine time_step(nx, q)
  integer, intent(in)             :: nx
  real(dp), intent(inout), device :: q(3,nx)    ! GPU 上の保存変数
  real(dp), allocatable, device   :: e(:,:), ev(:,:)
  integer :: n, j

  allocate(e(3,nx-1), ev(3,nx-1))    ! GPU メモリを動的確保（cudaMalloc に相当）

  do n = 1, nt
    call eev(nx, q, e, ev)            ! フラックス計算

    !$cuf kernel do(1)<<<32,4>>>
    do j = 2, nx-1                    ! 両端を固定（境界条件代わり）
      q(1:3,j) = q(1:3,j) - dtdx * (-e(1:3,j-1) + e(1:3,j) &
                                    -(-ev(1:3,j-1) + ev(1:3,j)))
    enddo
  enddo

  deallocate(e, ev)                   ! cudaFree に相当
end subroutine time_step
```

`j = 2` から `nx-1` のみ更新しており，両端 `j=1`, `j=nx` はそのまま（この問題では境界での更新が不要）．

---

## 5. `main.f90` の解説

```fortran
program main
  use cudafor
  use parameters, only : dp, nx, Rgas, gamma, Lx, rhol, a, Tlr
  use solver
  implicit none
  real(dp) :: q(3,nx)                         ! CPU メモリ上の保存変数
  real(dp), allocatable, device :: q_gpu(:,:) ! GPU メモリ上の保存変数

  call init(nx, q)                             ! CPU 上で初期条件設定

  ! 計時しながら各フェーズを実行
  q_gpu = q                                    ! CPU → GPU メモリ転送
  call time_step(nx, q_gpu)                    ! GPU 上で時間積分
  q = q_gpu                                    ! GPU → CPU メモリ転送

  deallocate(q_gpu)

  ! 結果をファイルに書き出し（無次元化）
  open(10, file='keep.dat')
  do n = 1, nx
    rho = q(1,n)
    u   = q(2,n) / rho
    p   = (gamma-1) * (q(3,n) - 0.5*rho*u**2)
    x   = dble(n) * dx
    write(10,"(4es16.8e3)") x/Lx, rho/rhol, u/a, p/(rhol*Rgas*Tlr)
  enddo
  close(10)
end program main
```

`q_gpu` は `main.f90` 内で宣言されており，`time_step` に渡すことで GPU 上に留まったまま時間積分が実行される．
GPU → CPU 転送（`q = q_gpu`）は計算終了後に一度だけ行われる．

---

## 6. 3バリアントの比較まとめ

### コード量

`eev` サブルーチンの記述量はほぼ同じ．カーネルバリアントだけ `kernel` サブルーチンが追加される分が増える程度．

### メモリアクセスパターン

| バリアント | 同変数隣接アクセス | Coalesced？ |
|---|---|---|
| AoS | `qq(1,j)` と `qq(1,j+1)` が 3 要素おき | 非 Coalesced（ストライドアクセス） |
| SoA / Kernel | `rho(j)` と `rho(j+1)` が隣接 | Coalesced |

### グリッド・ブロック設定

| バリアント | 設定 | 総スレッド数 | 備考 |
|---|---|---|---|
| AoS ディレクティブ | `<<<32,4>>>` | 128 | 1 スレッドが複数格子点を担当 |
| SoA ディレクティブ | `<<<32,4>>>` | 128 | 同上 |
| 明示的カーネル（フラックス） | `<<<32,128>>>` | 4096 | 1 スレッドが 1 界面を担当 |

カーネルバリアントのフラックス計算は 1 スレッド 1 界面の設計で，より GPU の並列性を活用している．

### 書き方の難易度

| バリアント | 書きやすさ | 最適化の自由度 |
|---|---|---|
| ディレクティブ (AoS/SoA) | 既存 Fortran コードに近い | 低い（コンパイラ依存） |
| 明示的カーネル | やや複雑（スレッドインデックス等） | 高い（Shared Memory 等を活用可能） |

---

## 7. 発展的なトピック

本サンプルからさらに学ぶ上で参考になる方向性：

- **Shared Memory の活用**: 隣接格子点の値を Shared Memory にキャッシュすると，グローバルメモリアクセスを削減できる
- **高次精度スキーム**: 現在は 2 次精度の粘性項離散化と 1 次の時間積分．4 次 Runge-Kutta などへの拡張が考えられる
- **多次元化**: 2D/3D への拡張では `do(2)` / `do(3)` ディレクティブや 2D/3D グリッド (`dim3`) が必要になる
- **MPI + CUDA Fortran**: 複数 GPU への分散並列化では MPI との組み合わせが典型的
