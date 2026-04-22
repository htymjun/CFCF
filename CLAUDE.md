# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CFCF (CUDA Fortran for Compressible Flow) demonstrates GPU-accelerated 1D compressible flow solvers using CUDA Fortran. It solves the Sod shock tube problem with the KEEP (Kinetic Energy and Enthalpy Preserving) numerical scheme for convective fluxes.

## Build Commands

Each variant lives in its own subdirectory under `shock_tube_keep/` and has its own Makefile:

```bash
cd shock_tube_keep/keep_directive_AoS   # or keep_directive_SoA or keep_kernel
make          # build → produces executable 'sod'
make clean    # remove .o, .mod, and executable
./sod         # run the solver
```

**Requirement:** NVIDIA HPC SDK >= 24 (`nvfortran` compiler) with a CUDA-capable GPU. Load the environment module `nvhpc-openmpi3` before building.

**Compiler flags used:** `-cuda -acc -fast -mfma -Mchkptr -Minfo -gpu=ptxinfo,rdc,lto`

## Output

Running `./sod` produces:
- `keep.dat` — solution columns: x/Lx, ρ/ρ_l, u/a, p/(ρ_l·R·T_lr)
- Timing output (CPU→GPU transfer, computation, GPU→CPU transfer)

## Architecture

### Three Variants

| Directory | GPU Approach | Data Layout |
|---|---|---|
| `keep_directive_AoS/` | OpenACC `!$cuf kernel do(1)` directives | `q(3,nx)` — variables interleaved |
| `keep_directive_SoA/` | OpenACC `!$cuf kernel do(1)` directives | `q(nx,3)` — variables contiguous |
| `keep_kernel/` | Explicit `attributes(global)` CUDA kernels | `q(nx,3)` |

All three variants share the same module structure:

- **`parameters.f90`** — all physical/numerical constants, grid size (`nx=4096`), time steps (`nt=2500`), Re=25000, CFL=0.1
- **`solver.f90`** — `init()` sets Sod shock tube initial conditions; `eev()` (directive variants) or `kernel()` (kernel variant) computes KEEP convective + viscous fluxes; `time_step()` runs 1st-order Euler time integration
- **`main.f90`** — allocates GPU device arrays, transfers CPU→GPU, calls `time_step()`, transfers GPU→CPU, writes output

### Conservative Variable Convention

- `q(1)` = ρ (density)
- `q(2)` = ρu (momentum)
- `q(3)` = E (total energy)

Primitive variables (ρ, u, p, T) are derived on-the-fly. Viscosity follows Sutherland's law.

### KEEP Scheme (core of `eev()`/`kernel()`)

Fluxes are computed at cell interfaces `i+1/2`:
- Mass: `c = 0.25*(ρ_i + ρ_{i+1})*(u_i + u_{i+1})`
- Momentum: `0.5*c*(u_i + u_{i+1}) + 0.5*(p_i + p_{i+1})`
- Energy: `c*(ie + ke) + pressure_dilatation`

Convective (`e`) and viscous (`ev`) flux arrays are sized `(3, nx-1)`.
