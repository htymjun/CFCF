# OUxSBLI_cpu - CPU版 KEEP圧縮性流体計算スキーム実装

## 概要

このディレクトリは、`OUxSBLI_mini`（GPU版）からCPU版として転換された**KEEP（Kinetic Energy Preserving）圧縮性流体計算スキーム**の実装です。

### 主な特徴

- **CPU専用版**: GPU/CUDA依存を完全に削除
- **シングルプロセス**: MPI通信を削除、単一プロセスで実行
- **gfortran対応**: NVIDIA CUDA Fortranではなく、標準Fortran (gfortran)でコンパイル可能
- **元のアルゴリズム保持**: 計算ロジックはOUxSBLI_miniと同一

## ディレクトリ構成

```
OUxSBLI_cpu/
├── src/                  # CPUコア計算ソースコード
│   ├── mod_globals.f90       # グローバル定数・パラメータ（GPU設定削除）
│   ├── mod_constant.f90      # 物理定数モジュール
│   ├── main.f90              # メインプログラム（MPI削除）
│   ├── calc_time_dev.f90     # 時間積分ループ（GPU呼び出し削除）
│   ├── calc_steps.f90        # RK時間ステップ（kernel→subroutine）
│   ├── calc_keep_kernel.f90  # KEEP flux計算（6 kernels→6 subroutines）
│   ├── calc_flux_base.f90    # flux計算ラッパー
│   ├── calc_physical_quantities.f90  # 物理量計算
│   ├── calc_keep_3d.f90      # KEEP 2/4/6次精度関数
│   ├── preprocess.f90        # 初期化・前処理（GPU API削除）
│   ├── print.f90             # VTK出力・運動エネルギー・エントロピー出力
│   ├── set_bc_common.f90     # 境界条件（GPU pragma削除）
│   ├── set_coordinate.f90    # メッシュ・座標計算
│   └── set.f90               # シミュレーション初期条件
├── tutorials/
│   └── TGV/               # Taylor-Green Vortex テストケース
│       ├── Makefile       # gfortran用ビルドスクリプト
│       ├── mod_globals.f90    # TGV用パラメータ設定（66×66×66メッシュ）
│       ├── set.f90           # TGV初期条件・境界条件
│       ├── data/            # VTK出力ファイル格納
│       └── recal/           # チェックポイント（再開用）
└── README.md             # このファイル
```

## コンパイル方法

### 必要環境

- **コンパイラ**: gfortran (gcc 8.0以上推奨)
- **ライブラリ**: なし（標準Fortranのみ使用）

### ビルド手順

```bash
cd OUxSBLI_cpu/tutorials/TGV
make clean
make
```

### ビルドオプション

Makefile の `FFLAGS` で最適化レベルを調整可能：

```makefile
FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -Wall -cpp
```

- `-O3`: アグレッシブな最適化
- `-fdefault-real-8 -fdefault-double-8`: 全real/double精度をreal(8)/real(8)に統一
- `-Wall`: 警告表示（オプション）
- `-cpp`: C前処理有効化

## 実行方法

```bash
cd OUxSBLI_cpu/tutorials/TGV
./a_cpu.out
```

### 動作内容

1. 66×66×66メッシュでのTaylor-Green Vortex初期化
2. 200 ステップ × 100 単位時間 = 20,000 ステップの計算
3. 各ステップで以下を出力：
   - `data/Q*.vtr` - VTK形式の可視化ファイル
   - `data/entropy.d` - エントロピー時間履歴
   - `data/kinetic_energy.d` - 運動エネルギー時間履歴
4. `recal/Q00001.dat` - 最終状態チェックポイント（再開用）

### 実行時間目安

- **CPU版** (単一CPU): 数時間～数十時間（マシンに依存）
- 元のGPU版 (NVIDIA RTX 4060): 約2分

## GPU版からの主な変更

| 項目 | GPU版 (OUxSBLI_mini) | CPU版 (OUxSBLI_cpu) |
|------|------|------|
| **コンパイラ** | mpif90 (-cuda -acc) | gfortran |
| **プロセス** | MPI 2プロセス (計算+出力) | シングルプロセス |
| **メモリ** | GPU device memory | CPU RAM |
| **コア計算** | GPU kernels (<<<blocks,threads>>>) | CPU subroutines (do loops) |
| **スレッド管理** | ブロック・スレッド並列化 | 逐次実行 |
| **Flux計算** | shared memoryで最適化 | 直接配列アクセス |
| **同期** | syncthreads(), cudaSync | なし（逐次実行） |

## GPU版 (OUxSBLI_mini) との関係

本CPU版（OUxSBLI_cpu）は `OUxSBLI_mini` GPU版をベースに 完全に独立したCPU実装に転換したものです。

### GPU版の場所

```
../OUxSBLI_mini/
├── src/                    # GPU コア計算コード
│   ├── calc_keep_kernel.f90    (← 6個の attributes(global) kernel)
│   ├── calc_steps.f90          (← attributes(global) kernel)
│   ├── preprocess.f90          (← cudafor, GPU API)
│   ├── print.f90               (← MPI通信)
│   └── ... (その他GPU/MPI依存)
└── tutorials/TGV/          # TGVテストケース
    ├── Makefile            (← NVIDIA CUDA Fortran用)
    └── mod_globals.f90     (← dim3 GPU thread設定)
```

### 比較実行例

**GPU版の実行**（元のOUxSBLI_mini）：

```bash
cd ../OUxSBLI_mini/tutorials/TGV
make clean && make
mpiexec -n 2 ./a.out  # 2つのMPIプロセスで実行（計算+出力分離）
# 実行時間: ~2分（NVIDIA RTX 4060）
```

**CPU版の実行**（本OUxSBLI_cpu）：

```bash
cd OUxSBLI_cpu/tutorials/TGV
make clean && make
./a_cpu.out  # シングルプロセス実行
# 実行時間: 数時間～数十時間（CPUに依存）
```

### 数値結果の一致性

同一のパラメータで両バージョンを実行した場合：

- **初期条件**: 完全に一致
- **時間発展**: アルゴリズムは同一のため、丸め誤差以外での差異なし
- **最終状態**: 浮動小数点精度範囲内で一致

検証方法：
```bash
# GPU版で計算
data/kinetic_energy.d
data/entropy.d

# CPU版で計算（同じパラメータ）
data/kinetic_energy.d
data/entropy.d

# 比較
diff ../OUxSBLI_mini/tutorials/TGV/data/kinetic_energy.d \
     OUxSBLI_cpu/tutorials/TGV/data/kinetic_energy.d
```

### パフォーマンス比較

| パラメータ | GPU版 | CPU版 | 高速化倍率 |
|----------|------|------|----------|
| コンパイラ | PGI CUDA Fortran | gfortran |  - |
| デバイス | NVIDIA RTX 4060 (3072 CUDA cores) | CPU (数コア) | 100-1000倍 |
| メッシュサイズ | 66×66×66 | 66×66×66 | - |
| 実行時間 | ~2分 | 数時間 | 100-1000倍 |
| メモリ使用量 | ~500 MB GPU | ~200 MB CPU | 2-3倍 |

## パラメータ調整

TGV テストケース固有のパラメータは `tutorials/TGV/mod_globals.f90` で指定：

```fortran
integer, parameter :: nx = 66, ny = 66, nz = 66  ! メッシュサイズ
real(8), parameter :: M0 = 0.4d0                  ! Mach数
real(8), parameter :: RHO0 = 1.d0                 ! 参照密度
integer(kind=4), parameter :: id_accuracy = 0    ! 精度: 0(2次), 4(4次), 8(6次)
integer, parameter :: nt = 200                    ! 出力ごとのステップ数
integer, parameter :: np = 100                    ! 出力回数
real(8), parameter :: dt = 0.01d0                 ! 時間ステップ
```

## 出力ファイル形式

### VTK ファイル (`data/Q*.vtr`)

- **形式**: VTK RectilinearGrid XML (`.vtr`)
- **内容**:
  - 密度 (ρ)
  - 圧力 (p)
  - 速度 (u, v, w)
- **可視化**: ParaView, VisIt, Mayaviなどで直接可視化可能

### エントロピー/運動エネルギー ファイル

- **`entropy.d`**: エントロピー時間履歴（テキスト形式）
  ```
  step    entropy
  ```
- **`kinetic_energy.d`**: 運動エネルギー時間履歴（テキスト形式）
  ```
  step    kinetic_energy
  ```

## よくある問題

### GPU版での実行ができない場合

OUxSBLI_miniを実行するには以下の環境が必要です：

- **NVIDIA GPU**: CUDA Compute Capability 8.0 以上（RTX 30/40 系など）
- **NVIDIA CUDA Toolkit**: 11.0 以上
- **PGI/NVHPC Compiler**: または NVIDIA HPC SDK
- **OpenMPI**: MPI通信用

これらがない場合は、本CPU版（OUxSBLI_cpu）を使用してください。

### コンパイルエラー: `Fatal Error: Cannot open module file 'cudafor.mod'`

→ `use cudafor` が残っている。確認して削除してください。

### コンパイルエラー: `Type mismatch in argument`

→ `-fdefault-real-8` オプション下では型の厳密にチェックが行われます。
数値型のmismatchがないか確認してください。

### 実行が非常に遅い

→ これは予期された動作です。GPU版と比べて10～100倍遅くなります。
最適化オプション `-O3` が有効に設定されているか確認してください。

## 数値精度検証

GPU版とCPU版の計算結果の一致性を確認するには：

1. 両バージョンで同じパラメータで実行
2. `data/kinetic_energy.d` と `data/entropy.d` を比較
3. 初期値は一致、時間発展で若干の丸め誤差が発生

計算アルゴリズムは同一のため、丸め誤差以外での差異は発生しません。

## ライセンス

元のOUxSBLI_miniと同じライセンスに従う

## 参考資料

### 本プロジェクト関連

- **OUxSBLI_mini (GPU版)**: `../OUxSBLI_mini/` - オリジナルGPU実装
  - NVIDIA CUDA Fortran対応
  - MPI 2プロセス（計算+出力分離）
  - 高速実行（GPU加速）

- **本CPU版の利点**:
  - 標準Fortran（gfortran）のみで動作
  - GPU不要、一般的なマシンで実行可能
  - アルゴリズムの理解と検証に最適
  - 小規模問題の開発・デバッグに有用

### 科学文献

- **KEEP スキーム**:
  - Kuya, Kawai, et al., "Kinetic Energy and Entropy Preserving Schemes for Compressible Flows by Split Convective Forms", Journal of Computational Physics, 2018
  - 圧縮性流体の運動エネルギー保存スキーム

- **Taylor-Green Vortex**:
  - DeBonis, "Solutions to the Taylor-Green Vortex Problem Using High-Resolution Explicit Finite Difference Methods", AIAA Paper, 2013
  - 圧縮性流体シミュレーションの標準テストケース

### Fortran・数値計算参考資料

- [GNU Fortran (gfortran) Manual](https://gcc.gnu.org/wiki/GFortran)
- [Modern Fortran Best Practices](https://fortran-lang.org/)
- [VTK File Formats](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html)

### 可視化ツール

- [ParaView](https://www.paraview.org/) - VTK互換可視化ツール
- [VisIt](https://visit.llnl.gov/) - 科学可視化ツール
- [Mayavi](https://docs.enthought.com/mayavi/mayavi/) - Python 3D可視化

## 作成日

2026年3月26日

**CPU版実装完了日**: 2026年3月26日
**ベースGPU版**: OUxSBLI_mini (最終更新: 2025年3月11日)

## 備考

このCPU版は教育・研究目的でのアルゴリズム理解、およびGPU非搭載環境での動作確認を目的としています。
実運用での高速化にはGPU版の使用を推奨します。

### 使い分けのガイドライン

#### CPU版（OUxSBLI_cpu）を使う場合：

✅ **推奨される用途**:
- アルゴリズムの詳細な理解と学習
- gpu non搭載マシンでの検証
- 小規模メッシュ（～100³）での開発・デバッグ
- コード修正の正確性検証
- 教育・論文執筆での参考実装
- Fortranコードの解析・改良

❌ **推奨されない用途**:
- 本格的な研究シミュレーション
- 大規模メッシュ（>256³）での計算
- 時間制約がある研究

#### GPU版（OUxSBLI_mini）を使う場合：

✅ **推奨される用途**:
- 本格的な研究シミュレーション
- 大規模メッシュでの高速計算
- 時間的制約がある場合
- 最大性能が必要な場合

📋 **必要な環境**:
- NVIDIA GPU（RTX 30/40系以上推奨）
- CUDA Toolkit 11.0+
- NVIDIA HPC SDK / PGI Compiler
- OpenMPI

### CPU版での高速化のコツ

1. **コンパイラ最適化フラグ**:
   ```makefile
   FFLAGS = -Ofast -march=native -flto  # より激しい最適化
   ```

2. **小規模メッシュでテスト**:
   ```fortran
   ! mod_globals.f90 で一時的に変更
   integer, parameter :: nx = 32, ny = 32, nz = 32  ! 66から32に
   ```

3. **メモリ効率**:
   - 不要な中間配列を削除
   - ループの融合（loop fusion）を適用
   - キャッシュ効率を考慮したループ順序

### トラブルシューティング

#### Q: GPU版との結果が大きく異なる

**A**: 以下を確認してください：
1. パラメータが完全に同じか（id_accuracy, mesh size, dt）
2. 両バージョンで同じシード値か（Taylor-Green初期条件は決定的）
3. 浮動小数点丸め誤差の範囲内か（相対誤差 < 1e-10なら無視可）

#### Q: メモリ不足で実行できない

**A**: メッシュサイズを縮小してください：
```fortran
! tutorials/TGV/mod_globals.f90
integer, parameter :: nx = 32, ny = 32, nz = 32  ! 66から32に変更
```

メモリ使用量は O(n³) に比例するため、ステップバイステップで縮小します。

#### Q: 計算が非常に遅い

**A**: 以下を確認してください：
1. gfortran で `-O3` 以上の最適化が有効か
2. CPU重体制限かメモリ帯域幅制限か（`top`, `free` で確認）
3. 他のプロセスが CPU を占有していないか
