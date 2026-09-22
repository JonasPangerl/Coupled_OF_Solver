# Spec amendment set C (verbatim)

Received from the user on 2026-09-22 ~13:50 ("stelle sicher dass das gebaut
wird und wir es zur verfuegung haben"); written by another agent from an
external coupled-solver comparison. Binding like amendment sets B/B11.

Lead notes (not part of the verbatim text):
- C1 lowers nonOrthThreshold 70 -> 65. The user set 85 on 2026-09-22 (D-047,
  explicit user decision, later than the external text). Default stays 85;
  65 and 60 become benchmark variants. To be confirmed with the user.
- C4 overlaps the existing `agglomerationWeights pressure` (D-039: p-p block
  weights); `agglomerator algebraicPair` is added as the keyword for it.
- Status per item: AMENDMENT_B_STATUS.md, section "Amendment C".

---

## Amendment set C — additions from external coupled-solver comparison

### C1 — Section 8.1 extended: topological static criteria
Add to the static set, evaluated once at start-up (and after any mesh change):
- **wallStarved**: a cell with at least one wall face and ≤ `maxWallInternalFaces`
  (default 2) internal faces (3D cases only, i.e. `mesh.nSolutionD() == 3`).
  Rationale: nearly enclosed layer cells have a poorly conditioned row.
- **procAMI**: a cell that has both a `processor` patch face and a
  `cyclicAMI` patch face. Rationale: lagged AMI coupling (5.5) and
  processor-interface effects overlap; the set grows with rank count.
- **volumeJump**: |V_P − V_N| / max(V_P, V_N) > `volJumpThreshold` (default
  0.98, i.e. ratio > 50) across any internal face — complements
  `volRatioThreshold`.
Each criterion has its own count in the once-per-run block (12.2) and its own
bit in `remediationFlag` (extend: 4 = wallStarved, 8 = procAMI, 16 =
volumeJump; flags are OR-ed). Threshold `nonOrthThreshold` default lowered
from 70 to **65**; benchmark variant with 60 on T4b/T5 (14.2 configuration E).
Keywords: `remediation.static.maxWallInternalFaces 2; procAMI yes;
volJumpThreshold 0.98;`

### C2 — Section 5.3(d) replaced: tensorial Rhie–Chow diffusivity
`rhieChow.tensorial yes` (default yes). Per cell, let A_P be the 3×3
momentum diagonal block after step 5.3(a) (including PTC and implicit MRF
Coriolis, excluding the pressure column). Compute
  D_P = V_P · A_P^{-1}   (3×3, double, once per outer iteration)
with `block4Ops::invert3` (LU, partial pivoting). Guard: if
|det(A_P)| < 1e-12 · (trace(A_P)/3)^3, use the Moore–Penrose pseudoinverse
via 3×3 SVD (Jacobi, double) and count `nPseudoInverse`; a non-zero count is
a WARNING every 100 iterations and a test failure on T0–T3.
Face coefficient: D_f = n_f · interpolate(D)_f · n_f   (scalar, n_f = S_f/|S_f|),
used in the p–p Laplacian and in the flux update exactly where the scalar
D_f was used before. Interpolation of D is linear on the 3×3 entries.
With `tensorial no` the previous scalar definition applies (kept for the
benchmark comparison, 14.2 configuration E).

### C3 — New Section 7.6: Selective Frequency Damping (optional)
```
sfd { enabled no; chi 0.5; Delta 1.0; Lref 1.0; deactivateBelowR 1e-4; }
```
Purpose: suppress limit-cycle oscillations of the steady iteration
(non-dimensional χ, Δ scaled with U_ref (7.2) and Lref).
  χ* = chi · U_ref / Lref   [1/s],   Δ* = Delta · Lref / U_ref   [s]
Filtered velocity Ū (volVectorField, part of coupledState): initialised
Ū = U at activation; updated after every accepted outer step with the local
pseudo-time step of 5.4:
  Ū_P ← Ū_P + (Δt_P / (Δ* + Δt_P)) · (U_P − Ū_P)     (implicit Euler low-pass)
Momentum contribution (5.3a, before PTC): a_P += χ* V_P on the three
momentum diagonals; b_P += χ* V_P Ū_P. The forcing −χ*(U − Ū) vanishes at
the steady state, so the converged solution is unbiased.
Deactivation: when R_n < deactivateBelowR for nHold consecutive iterations,
set χ* = 0 (log event `SFD-off`) and keep Ū in state for restart.
Interaction: Ū is included in rollback (9.3) and flushed with the Anderson
history (7.5). Parameter guidance in the log header: χ ≈ 2× the
non-dimensional growth rate of the oscillation, Δ ≈ 1/(2 f_osc·Lref/U_ref);
if unknown, start with the defaults.
Test: T3-SST with `enabled yes` must converge to the same Cl/Cd as
`enabled no` within 1e-3 (unbiasedness check); report the limit-cycle
amplitude of Cl over the last 200 iterations for both (15.4).
Ref: Åkervik, Brandt, Henningson, Hœpffner, Marxen, Schlatter, *Steady
solutions of the Navier–Stokes equations by selective frequency damping*,
Phys. Fluids 18 (2006).

### C4 — Section 6.3.1 addition: agglomerator option
`agglomerator algebraicPair` selectable; coefficients = p–p block (row 3,
col 3) magnitudes, so agglomeration follows the anisotropy of the Rhie–Chow
Laplacian (stretched layer cells merge across their thin direction).
Benchmark variant on T4b/T5 vs faceAreaPair (14.2, E); the faster one on
T5 becomes the default and is recorded in DECISIONS.md. Coarsening-ratio
rule 6.3.1 applies unchanged.

### C5 — Section 8.3 extended: regex patch sets and built-in zones
`remediation.zonal.patchDistance[i].patches` accepts regular expressions
(native `wordRes`). Built-in named sets usable in `zones`: `_wallCells`
(all cells with a wall face), `_procCells` (cells with a processor face),
`_amiCells`, `_procAMICells`, `_wallStarved`, `_remediationStatic`.
A cell affected by several entries gets the **product** of the factors; the
effective per-cell multipliers are written as `volScalarField cflFactorEff`
and `betaEff` at every write time so the combined effect is always visible.

### C6 — Section 12 addition: diagnostic fields
At every write time (and on abort): `localDt` (Δt_P after all caps),
`localCFL` (Δt_P · (½Σ|φ_f|)/V_P), `cflFactorEff`, `betaEff`,
`remediationFlag`, and when SFD is active `USFD` (= Ū).

### C7 — Section 14.2 configuration E extended
Add variants: `rhieChow.tensorial no`, `nonOrthThreshold 60`,
`agglomerator algebraicPair`, `sfd.enabled yes` (T3 only), and
`etaMax 0.7` with `minIter 2` (tests the "solve loosely, iterate often"
design point against the default).

### C8 — DECISIONS.md entries (do not implement)
- No solve-on-variables form: the increment form is required for line
  search, adaptive tolerance and rollback (Sections 5.7, 7.2, 9.3).
- No multiplicative relaxation stack: local pseudo-time replaces it; the
  (2−r) approximation is numerically wrong for r < 0.8 and not adopted.
- Static-only remediation is insufficient: dynamic set (8.2) stays.
- Post-solve pressure shift for closed domains not adopted; row
  replacement (singular system avoided) stays.
