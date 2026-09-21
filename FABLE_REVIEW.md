# Open problems for a follow-up review (Fable)

Started 2026-09-21 ~23:10 at the user's request: everything that causes
problems overnight and is not solved cleanly is listed here with evidence
and a proposed fix, so a stronger model can pick it up in the morning.
Resolved items move to the "Resolved" section with the commit.

## Open

### 1. Round-off-sized change doubles T0 Re 1000 iterations (np1)

Evidence (debug agent, 2026-09-21, A/B on T0 Re 1000 np1):

| orthogonality test in nonOrthCorrection | FGMRES | outer iterations |
|---|---|---|
| exact `kMax == 0` (committed) | old | 58 |
| exact | mixed precision (committed) | 57 |
| `kMax <= 1e-10` | old | 123 (5 CFL cuts) |
| `kMax <= 1e-10` | mixed precision | 106 (6 CFL cuts) |

The T0 blockMesh cavity has max|k| = 7.1e-14 (round-off), so the only
difference is whether the non-orthogonal correction path runs with
round-off-sized vectors or is skipped. A correction of 1e-14 cannot
change the physics. Most likely `orthogonal_` also switches something
that is not just the k-correction: a different delta-coefficient field,
a different gradient or explicit-source path, or the deferred correction.
It may also be that Re 1000 is extremely sensitive at the start-up switch
(beta 0 -> 1, D-008).

Related: T0 Re 1000 at np4 needs 182-282 outer iterations (np1: 57),
with or without the tolerance. A coupled Newton-like solver should be
nearly independent of the decomposition; only the preconditioner
differs. This suggests the outer convergence at Re 1000 depends on
inner-solve quality (Eisenstat-Walker eta, K-cycle at np4) more than
expected.

To do:
- (a) Diff the code paths guarded by `orthogonal_` / `nonOrthCorrection`
  and find what else the flag controls.
- (b) Log per-iteration R, CFL, eta and cuts for np1 vs np4 at Re 1000
  and find where they diverge (likely around the start-up switch).
- (c) Decide the tolerance: the orthogonality tolerance is conceptually
  right (coupledDefaults::orthogonalityTolerance exists but is unused).
Committed state: exact test (f025546).
