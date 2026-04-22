# CUDA Fortran for Compressible Flow (CFCF)

CUDA Fortranを用いたGPU加速CFDソルバのサンプルコード集．
KEEPスキームで粘性ソッド衝撃波管問題を解く3つの実装バリアントを収録している．

## 依存

- NVIDIA HPC SDK >= 24（`nvfortran` コンパイラ）
- CUDA 対応 GPU
- 環境モジュール: `nvhpc-openmpi3` をロードする

## サンプルコード

KEEPスキームを使用して粘性衝撃波管問題を解く3つのバリアント：

| ディレクトリ | GPU の使い方 | メモリレイアウト |
|---|---|---|
| `shock_tube_keep/keep_directive_AoS/` | `!$cuf kernel do` ディレクティブ | AoS: `q(3, nx)` |
| `shock_tube_keep/keep_directive_SoA/` | `!$cuf kernel do` ディレクティブ | SoA: `q(nx, 3)` |
| `shock_tube_keep/keep_kernel/` | 明示的 CUDA カーネル関数 | SoA: `q(nx, 3)` |

```bash
cd shock_tube_keep/keep_directive_AoS  # またはほかの2つ
make
./sod
```

実行すると `keep.dat`（無次元化された解）と計算時間が出力される．

## 勉強会資料

`./doc` に CUDA Fortran 学習者向けの解説資料を収録している：

| ファイル | 内容 |
|---|---|
| [doc/index.md](doc/index.md) | 目次・学習ガイド |
| [doc/01_cuda_fortran_intro.md](doc/01_cuda_fortran_intro.md) | CUDA Fortranの基礎（スレッド階層，`device`属性，カーネル起動） |
| [doc/02_physics.md](doc/02_physics.md) | 物理モデルと数値手法（KEEP スキーム，N-S 方程式） |
| [doc/03_code_explanation.md](doc/03_code_explanation.md) | コード詳細解説と3バリアントの比較 |

## CUDA Fortranの公式ドキュメント＆書籍
最新機能を使いたい場合(CUDA FortranはCUDA C++と比べて新機能の実装が遅れている)や生成AIに書かせたコードが動かない場合は公式ドキュメントを参照するのがよい．

https://docs.nvidia.com/hpc-sdk/compilers/fortran-cuda-interfaces/#

https://docs.nvidia.com/hpc-sdk/compilers/cuda-fortran-prog-guide/#

この書籍にも細かな情報が記載されているため有用．

@book{ruetsch2024cuda,
  title={CUDA Fortran for scientists and engineers: best practices for efficient CUDA Fortran programming},
  author={Ruetsch, Gregory and Fatica, Massimiliano},
  year={2024},
  publisher={Elsevier}
}
