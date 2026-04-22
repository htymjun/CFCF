# CUDA Fortran 入門

## 1. GPU アーキテクチャの概要

### CPU と GPU の違い

CPU は少数の高性能コアを持ち，複雑な制御フローや逐次処理に適している．
一方 GPU は数千～数万の小型コアを持ち，同じ演算を大量のデータに並列に適用する処理（SIMT: Single Instruction Multiple Threads）に特化している．

CFD のような格子点ごとに独立した計算が多いアプリケーションは，GPU の並列性を活かしやすい．

### CUDA のスレッド階層

GPU 上の計算は **スレッド (thread)** の集合として実行される．

```
Grid（グリッド）
└─ Block（ブロック） × 複数
   └─ Thread（スレッド） × 複数
```

- **Thread**: 最小の実行単位．各スレッドが格子点 1 点を担当するように設計するのが典型的．
- **Block**: Thread のグループ．同一 Block 内のスレッドは共有メモリ (shared memory) を使える．
- **Grid**: Block の集合．1 つのカーネル起動がひとつの Grid に対応する．

サンプルコードでは 1 次元の格子に対応して，1D の Grid / Block 構成を使っている．

### メモリの種類

| メモリ種別 | 容量 | 速度 | 備考 |
|---|---|---|---|
| グローバルメモリ | 数 GB | 遅い | GPU上のメインメモリ．`device`変数はここ |
| 共有メモリ | 数十 KB/Block | 速い | 同一Block内スレッドが共有 |
| レジスタ | 数十 KB/スレッド | 最速 | スレッドローカル変数 |
| ホストメモリ (CPU) | 数十 GB | 非常に遅い | CPU⇔GPU転送が必要 |

本サンプルではグローバルメモリのみを使用している．

---

## 2. CUDA Fortran の基礎

### `cudafor` モジュール

CUDA Fortran の機能は `cudafor` モジュールで提供される．

```fortran
use cudafor
```

これにより，`device` 属性，CUDA 組み込み変数 (`threadIdx`, `blockIdx`, `blockDim`)，カーネル起動構文などが使えるようになる．

### `device` 属性

GPU デバイスメモリ上に配置する変数には `device` 属性を付ける．

```fortran
real(dp), device :: x_gpu(nx)         ! 静的確保（GPU メモリ）
real(dp), allocatable, device :: y_gpu(:)  ! 動的確保
allocate(y_gpu(nx))
```

`device` 属性のない変数はすべてホスト (CPU) メモリ上にある．

### CPU ⇔ GPU メモリ転送

CUDA Fortran では代入演算子 `=` で CPU/GPU 間のメモリ転送が書ける（暗黙的転送）．

```fortran
real(dp) :: q_host(3, nx)          ! ホスト変数
real(dp), device :: q_gpu(3, nx)   ! デバイス変数

q_gpu = q_host   ! ホスト → デバイス転送 (cudaMemcpy H2D)
q_host = q_gpu   ! デバイス → ホスト転送 (cudaMemcpy D2H)
```

内部では `cudaMemcpy` が呼ばれる．サンプルコードの `main.f90` でこのパターンを使っている：

```fortran
! main.f90（AoS バリアント）
q_gpu = q    ! CPU → GPU
call time_step(nx, q_gpu)
q = q_gpu    ! GPU → CPU
```

---

## 3. CUF カーネルディレクティブ

### `!$cuf kernel do` ディレクティブ

既存の Fortran ループを GPU 上で並列実行する最も簡単な方法．ループの直前に書く．

```fortran
!$cuf kernel do(1)<<<grid_size, block_size>>>
do j = 1, nx
  ! このループ本体が GPU 上で並列実行される
  a(j) = b(j) + c(j)
enddo
```

- `do(1)`: 1 次元のループを並列化する
- `<<<grid_size, block_size>>>`: グリッドサイズとブロックサイズを指定

サンプルコードでは `<<<32, 4>>>` を使っている（合計 128 スレッド）：

```fortran
!$cuf kernel do(1)<<<32,4>>>
do j = 1, nx
  qq(1,j) = q(1,j)
  qq(2,j) = q(2,j) / qq(1,j)
  ...
enddo
```

### グリッド・ブロックサイズの選び方

グリッドサイズ × ブロックサイズ = 総スレッド数．総スレッド数がループの反復回数以上になるように設定する必要がある．

サンプルコードの `nx = 4096` に対して `<<<32, 4>>>` では 128 スレッドしかなく，各スレッドが複数の格子点を担当する（コンパイラが自動的にループを分割する）．

実用的なコードでは通常 `<<<(nx+255)/256, 256>>>` のように `nx` に応じて計算する．

### ディレクティブのメリット・デメリット

| | メリット | デメリット |
|---|---|---|
| **ディレクティブ** | 既存 Fortran コードへの変更が少ない | スレッド制御が粗い，最適化の自由度が低い |
| **明示的カーネル** | 完全なスレッド制御が可能 | 記述量が増える |

---

## 4. 明示的 CUDA カーネル関数

### `attributes(global)` サブルーチン

GPU 上で実行されるカーネル関数は `attributes(global)` 属性で宣言する：

```fortran
attributes(global) subroutine my_kernel(n, a, b, c)
  integer, intent(in), value :: n
  real(dp), intent(in),  device :: a(n), b(n)
  real(dp), intent(out), device :: c(n)
  integer :: j
  j = (blockIdx%x - 1) * blockDim%x + threadIdx%x
  if (j > n) return
  c(j) = a(j) + b(j)
end subroutine my_kernel
```

重要なポイント：
- **スカラー引数は `value` 属性が必要**（値渡し）．配列引数は参照渡しで `device` 属性をつける．
- **`threadIdx%x`, `blockIdx%x`, `blockDim%x`** などの CUDA 組み込み変数でグローバルスレッドインデックスを計算する．
- **境界チェック** (`if (j > n) return`) が必要．グリッドサイズの繰り上げで余分なスレッドが生まれるため．

### スレッドインデックスの計算

```
グローバルスレッドインデックス j:

  blockIdx%x  : 何番目のブロックか（1 始まり）
  blockDim%x  : 1 ブロックあたりのスレッド数
  threadIdx%x : ブロック内でのスレッド番号（1 始まり）

  j = (blockIdx%x - 1) * blockDim%x + threadIdx%x
```

例: 128 スレッド(=4ブロック×32スレッド)の場合

| blockIdx | threadIdx | j |
|---|---|---|
| 1 | 1 | 1 |
| 1 | 32 | 32 |
| 2 | 1 | 33 |
| 4 | 32 | 128 |

### カーネルの起動

```fortran
call my_kernel<<<grid_dim, block_dim>>>(n, a, b, c)
```

`grid_dim`, `block_dim` には整数またはネスト型 (`dim3`) を指定する：

```fortran
! 整数で指定（1次元）
call kernel<<<nx/128, 128>>>(nx, ...)

! dim3 で指定（多次元も可）
type(dim3) :: grid, block
grid  = dim3(nx/128, 1, 1)
block = dim3(128,    1, 1)
call kernel<<<grid, block>>>(nx, ...)
```

サンプルコードの `keep_kernel` バリアントでは：

```fortran
call kernel<<<nx/128,128>>>(nx, rho, u, p, T, mu, e, ev)
```

`nx = 4096` なので `4096/128 = 32` ブロック，各 128 スレッド，合計 4096 スレッドで 4095 個の界面を処理する．

---

## 5. コンパイル

### nvfortran のフラグ

サンプルコードの Makefile で使われているフラグ：

```makefile
FFLAGS = -cuda -acc -fast -mfma -Mchkptr -Minfo -gpu=ptxinfo,rdc,lto
```

| フラグ | 意味 |
|---|---|
| `-cuda` | CUDA Fortran を有効化（`device` 属性，カーネル起動構文など） |
| `-acc` | OpenACC を有効化（`!$cuf kernel do` はこれも必要） |
| `-fast` | 高レベル最適化（`-O2` + ベクトル化など） |
| `-mfma` | Fused Multiply-Add (FMA) 命令を使う |
| `-Mchkptr` | ポインタ境界チェック |
| `-Minfo` | 最適化情報を標準エラーに出力 |
| `-gpu=ptxinfo` | PTX（GPU 中間コード）の情報を出力 |
| `-gpu=rdc` | Relocatable Device Code（デバイスコードの分割コンパイル） |
| `-gpu=lto` | Link Time Optimization（リンク時最適化） |

### コンパイル順序

Fortran はモジュールの依存関係に従って順番にコンパイルする必要がある：

```
parameters.f90  →  solver.f90  →  main.f90
```

`solver.f90` は `parameters` モジュールを `use` し，`main.f90` は `parameters` と `solver` の両方を `use` するため，上記の順でコンパイルしなければならない．サンプルの Makefile はこの順序を保証している：

```makefile
SRC = parameters.f90 solver.f90 main.f90
OBJ = $(SRC:.f90=.o)
```
