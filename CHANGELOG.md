# Changelog

All notable changes. Code marked *uncompiled* has only passed
`scripts/syntax-check.sh` (g++ -fsyntax-only) and has not been built with
wmake or run yet.

## [Unreleased]

### Phase 0 / A - 2026-09-21
- Repository skeleton (spec Section 4), README, PLAN, DECISIONS (D-001..D-014),
  spec copy (`SPEC_coupledFoam.html` authoritative, `.md` text extraction).
- `src/blockMatrix`: `blockScalar.H` (float storage / double reductions,
  D-001), `block4Ops.H` (4x4 kernels, LU inversion with pivot guard),
  `doubleReduce` (double-accumulated reductions, MPI_DOUBLE),
  `blockLduInterface` (non-blocking processor exchange of 4-vectors),
  `blockLduMatrix4` (4x4-block LDU storage, Amul, row scaling). *uncompiled*
- `scripts/syntax-check.sh`: syntax-only compile with project warning flags.
