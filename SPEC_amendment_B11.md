# Spec amendment B11 (verbatim)

Received from the user on 2026-09-21, quoted verbatim below; supplied by
another agent, to be incorporated like amendment set B (binding). Status:
AMENDMENT_B_STATUS.md.

---

### B11 — New Section 6.5: hot-loop performance rules (mandatory for src/blockMatrix and src/blockSolvers)

1. **No OpenFOAM field algebra inside solver kernels.** Amul, smoothers, restriction/prolongation, and all Krylov vector operations are hand-written loops over flat `solveScalar*` arrays obtained via `.data()`; never `Field` operators, never `tmp<>`, inside these loops.
2. **Fused Krylov kernels.** Implement the following as single-pass loops (one read of each operand): `axpy_dot` (y += a·x; return <y,y> and <y,z>), `update_residual_norm` (r −= a·q; return ‖r‖²), and the FGMRES modified-Gram–Schmidt step as a fused "dot then axpy" per basis vector. Reference: fused/pipelined Krylov kernels reduce memory traffic by 1.5–2× on bandwidth-bound hardware. Reductions still accumulate in double (6.0).
3. **Vectorisation.** Every kernel loop carries `#pragma omp simd` (compiled with `-fopenmp-simd`, no OpenMP runtime), pointers are declared `__restrict__`, arrays are allocated 64-byte aligned (`std::aligned_alloc` or `posix_memalign` in a small `alignedList<solveScalar>` wrapper). A 4×4 float block (64 B) is exactly one AVX-512 register / two AVX2 registers; the block4Ops matvec is written so the compiler sees a fixed-trip-count 4×4 loop (no runtime loop bounds).
4. **Data layout is fixed** as in 5.2 (block AoS): one 4×4 block contiguous per face/cell. Do not introduce SoA variants.
5. **Assembly writes directly into the flat arrays.** `coupledAssembler` computes face contributions in one pass over `mesh.owner()/neighbour()` and writes into `upper_/lower_/diag_` in place; native `fvm::` operators are used only to *obtain* scalar coefficient fields (laplacian, upwind div), never re-assembled through `fvMatrix` addition.
6. **Measure, don't guess.** Phase A gate adds `Test-kernelBandwidth`: run Amul and the fused axpy_dot on a 5 M-cell synthetic system, report achieved GB/s vs. STREAM triad on the same machine; pass: Amul ≥ 60 % of STREAM, axpy_dot ≥ 80 %. Results go to `results/gates/phase_A.json` and into the report (15.5).
7. Expression templates are explicitly out of scope (DECISIONS.md entry: hot loops are hand-written; field algebra is < 30 % of iteration time and lives in OpenFOAM core).
