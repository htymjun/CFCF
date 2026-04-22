# 圧縮性CFDのための CUDA Fortran 勉強会資料

本資料は，CUDA Fortran を用いた GPU 加速 CFD (Computational Fluid Dynamics) コードの作成・理解を目的とした勉強会資料です．
サンプルコードとして，KEEPスキームで粘性ソッド衝撃波管問題を解くコードを題材にしています．

## 対象読者

- Fortran の基本文法を知っている
- GPU 計算・CUDA Fortran を学びたい
- CFD コードへの GPU 適用方法を理解したい

## 資料一覧

| ファイル | 内容 |
|---|---|
| [01_cuda_fortran_intro.md](01_cuda_fortran_intro.md) | CUDA Fortranの基礎概念・文法 |
| [02_physics.md](02_physics.md) | 物理モデル・数値手法 (KEEP スキーム) |
| [03_code_explanation.md](03_code_explanation.md) | サンプルコードの詳細解説と3バリアントの比較 |

## サンプルコードの構成

```
shock_tube_keep/
├── keep_directive_AoS/   # CUF カーネルディレクティブ + AoS メモリレイアウト
├── keep_directive_SoA/   # CUF カーネルディレクティブ + SoA メモリレイアウト
└── keep_kernel/          # 明示的 CUDA カーネル関数 + SoA メモリレイアウト
```

3つのバリアントはすべて同じ物理問題・数値手法を解くが，GPU プログラミングのスタイルとデータ配置方法が異なる．
各バリアントの違いを理解することで，CUDA Fortran の主要なプログラミングパターンを一通り学べる．

## 推奨学習順序

1. **[02_physics.md](02_physics.md)** — まず解くべき物理問題と数値手法を理解する
2. **[01_cuda_fortran_intro.md](01_cuda_fortran_intro.md)** — CUDA Fortran の基礎文法・概念を学ぶ
3. **[03_code_explanation.md](03_code_explanation.md)** — サンプルコードを読み解き，3バリアントを比較する

## ビルド・実行

```bash
# 例: AoS バリアントをビルドして実行
cd shock_tube_keep/keep_directive_AoS
make
./sod
```

必要環境: NVIDIA HPC SDK >= 24 (`nvfortran` コンパイラ)，CUDA 対応 GPU
