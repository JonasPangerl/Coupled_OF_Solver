# F1 half-car, first coupledFoam run: three defects

> **Status 2026-09-30 afternoon: the defects below are fixed in the working
> tree** (DECISIONS.md D-074..D-079, CHANGELOG.md). With them the F1 case
> ran 300 iterations on the T4 template settings with 0 rollbacks. What
> remained was a small effective pseudo-time step (CFL*omega 1-4); D-075
> was narrowed to the remediation cells and D-079 (local step control)
> added for that. The run records are in the OpenFoam_cases repository,
> `f1_halfcar/runs/s5_ec2_coupledFoam/`. One correction to the text below
> is marked in section 4.

Written 2026-09-30 for the agent who will fix the solver. Everything below
was measured on this machine today; no number is quoted from memory or from
an earlier document. Read `F1_READINESS.md` first - it predicted two of the
three areas that broke.

**Bottom line:** the segregated reference run is complete. coupledFoam does
not survive iteration 18 on this case. The blocking defect is a local
pressure blow-up on the symmetry plane that gets WORSE as the step gets
smaller, which rules out a step-size problem and points at the assembly.

---

## 1. Where everything is

Machine: c7a.16xlarge, 64 physical cores, 123 GB RAM, OpenFOAM v2606
(system install), coupledFoam built at
`~/OpenFOAM/ubuntu-v2606/platforms/linux64GccDPInt32Opt/bin/coupledFoam`.
The machine has a watchdog that powers it off after 30 min without a
compute process (`/usr/local/bin/cfd-watchdog`, `status`, `verlaengern <h>`);
the disk survives that, all of the below was still there after yesterday's
shutdown.

    ~/f1/f1_halfcar_s5_it0000_startCase   shipped case, NEVER modified
    ~/f1/runs/s5_seg_np60                 simpleFoam reference, COMPLETE
    ~/f1/runs/s5_cpl_np60                 coupledFoam, 5 runs, all logged
    ~/f1/runs/potentialStart              the shared start state, untouched
    ~/f1/runs/README.md                   how the two cases were built

Mesh: 20,614,106 cells, decomposed once with scotch into 60 subdomains
(`SCOTCH_PTHREAD_NUMBER=1`) in the segregated case; the coupled case uses
the same decomposition through hard links, so both solvers see a
bit-identical decomposed mesh. Roughly 340 k cells per rank.

`potentialFoam -writephi` ran ONCE (7.5 s, continuity error 1.6e-4) and its
result is the start state for both solvers; `~/f1/runs/potentialStart` is
an untouched copy of it.

### Reproducer

    cd ~/f1/runs/s5_cpl_np60
    rm -rf processor*/[1-9]* processor*/*_lastValid
    CFENV_CWD=$PWD ~/bin/cfenv sys mpirun -np 60 coupledFoam -parallel \
        2>&1 | tee log.new

Dies at iteration 18, about 2 minutes in. The five logs of today's runs are
in that directory as `log.coupledFoam*`; the naming says what each one was.

---

## 2. The segregated reference (complete, use it as the target)

300 iterations, 1916 s wall (6.4 s/iteration) on 60 ranks.

| quantity | value |
|---|---|
| Cd at iteration 300 | **1.238505** |
| Cl at iteration 300 | **-3.447526** |
| Cd mean over 290-300 | 1.245900 |
| Cl mean over 290-300 | -3.446178 |

Not converged at 300 - Cd was still falling by ~0.0015 per iteration. The
budget was capped at 300 deliberately (both `runTimeControl` objects removed
from `system/controlDict`) so that both solvers get the same budget.

**Field extremes of that solution**, measured with
`simpleFoam -postProcess -func "fieldMinMax(fields=(p U),location=true)"`:

| quantity | value | Cp = p/pref | location |
|---|---|---|---|
| max \|U\| | 150.1 m/s | 3.0 U_inf | (3.437, -0.903, 0.001) rear wheel contact |
| min p | -14664.3 | **-11.7** | (3.558, -0.929, 0.006) rear wheel contact |
| max p | 1801.3 | +1.44 | (-0.081, -0.827, 0.0003) front wheel contact |

This is the yardstick: **the physical suction minimum on this geometry is
Cp = -11.7**, which matches the `remediation.dynamic.cp 15` setting in the
T4 template. Anything far beyond that is numerical.

---

## 3. DEFECT 1 (blocking): pressure blow-up on the symmetry plane

### What happens

With a correct `Uref` the run reaches iteration 13 healthily (R falls from
1 to 0.021), then the sentinel rejects every subsequent step because `min p`
runs away in a cluster of cells a few centimetres across at

    (3.40 .. 3.44, -0.007 .. -0.001, 0.427 .. 0.447)     all on rank 0

i.e. **on the symmetry plane (y ~ 0), at mid height, at the rear wing**.
After 5 consecutive rollbacks the run aborts. `max|U|` at that moment is
195-268 m/s, nowhere near the velocity limit - it is purely the pressure
criterion (`sentinel.H:17`, `max|p| > pFactor*pref`) that fires.

### The decisive observation

Two runs, identical except for `ptc.CFLmin`. Rollback reduces CFL by
`sentinel.cflFactor` on every rejection, so the two differ only in how small
the retry step may become:

CFL below is the value the rollback cut down TO, as reported on the
`sentinel rollback` line; min p is the value of the step that was rejected.

| iteration | min p, CFLmin 1 | CFL -> | min p, CFLmin 0.05 | CFL -> |
|---|---|---|---|---|
| 14 | -78383 | 1.50 | -78383 | 1.50 |
| 15 | -224401 | 1.00 | -224401 | 0.375 |
| 16 | -224401 | 1.00 | -224401 | 0.094 |
| 17 | -224401 | 1.00 | **-795833** | 0.05 |
| 18 | -224401 | 1.00 | -795772 | 0.05 |

(`log.coupledFoam.pSentinel_CFLmin1` and `log.coupledFoam` respectively.)

With `CFLmin 1` the cut saturates on the floor, so every retry repeats the
same step: R at iterations 16 and 17 is 0.022332884 and 0.022332841,
identical to six significant figures, and the rejected state's min p is the
same -224401.11 at iterations 15 through 18. The rollback machinery cannot
do anything at all. That alone is worth fixing: a floor the recovery path
reaches is a floor that disables the recovery path.

With the floor lowered, the step really does shrink (CFL 1.5 -> 0.05, and
the line search damps to its own floor `omega 0.1`, `cuts=2`) - **and the
pressure spike gets ten times worse instead of better.** min p goes from
-78383 to -795833, i.e. Cp -63 to Cp -637, against a physical minimum of
Cp -11.7.

A step-size instability improves when the step shrinks. This does the
opposite. So the pseudo-transient term is not what is holding this together,
and the problem is in what is assembled at those faces.

### Why the symmetry plane is the first suspect

- The cluster sits at y ~ -0.002, i.e. in the first cell layer off the
  symmetry plane, and nowhere else in 20.6 M cells.
- `F1_READINESS.md` lists symmetry/symmetryPlane as "dedicated zero-flux
  patches" - a separate code path from the generic BC linearisation.
- The FPE of defect 2 lands in the same neighbourhood of the code:
  `nonOrthCorrection::correctionFromGrad` <- `rhieChow::updateExplicit`
  <- `coupledAssembler::assembleContinuity`.
- simpleFoam solves the same cells on the same mesh with no trouble at all.

### Mesh quality is a co-factor, not an alibi

`checkMesh -constant` on this mesh **fails one check**
(`~/f1/runs/s5_cpl_np60/log.checkMesh`):

- max skewness 11.381773, **126 highly skew faces** (the solver's own
  `remediation.meshQuality.skewThreshold` is 6, so these are inside the
  static set)
- max non-orthogonality 74.92 deg, 1472 severely non-orthogonal faces
- 203 non-manifold points

coupledFoam already puts **4858 cells** under static remediation on this
mesh (`nStat=4858`, constant across every run) and blows up anyway.

**Not yet answered, and it is the first thing to check:** whether the
failing cluster at (3.40, 0, 0.43) is among those 126 skew faces. I did not
resolve it - `coupledStaticSet` only prints counts, not locations, and
getting coordinates out of the `skewFaces` faceSet needs a few lines of
code. If the cluster IS skew, the question becomes why full static
treatment does not hold it; if it is NOT, the assembler is wrong on clean
cells and that is a much more serious bug.

---

## 4. DEFECT 2: `potentialClip` dies with an FPE at this mesh size

`potentialInit yes; potentialClip 4;` (D-057, the mechanism that exists
precisely for singular potential-flow peaks) crashes in the **first**
continuity assembly, before iteration 1. Log: `log.coupledFoam.fpeClip`.

    coupledFoam: potentialClip 4: 6106 cells clipped to |U| = 200 (peak 961.03951)
    ...
    Foam::sigFpe::sigHandler(int)
    Foam::nonOrthCorrection::correctionFromGrad(...)
    Foam::rhieChow::updateExplicit(...)
    Foam::coupledAssembler::assembleContinuity(...)

The mechanism is visible in `applications/coupledFoam/coupledFoam.C:813-843`:
the clip keeps the total head by moving the kinetic energy into p,

    Ui[celli] *= Umax/mu;
    pi[celli] += 0.5*(sqr(mu) - sqr(Umax));

so the peak cell gets `p += 0.5*(961^2 - 200^2) ~ 4.4e5` while its untouched
neighbour stays near 0. Over 6106 scattered cells that is a pressure field
the non-orthogonality correction does not survive. **The clip trades a
velocity spike for a pressure spike, and at 20.6 M cells the pressure spike
is the worse of the two.** Whatever the fix is, it probably has to smooth
the head redistribution over a neighbourhood instead of dumping it in one
cell, or clip U without touching p at all.

**Correction (added later the same day):** an earlier version of this
section claimed that coupledFoam's own potential start peaks at 961 m/s
where potentialFoam's peaks at 679, i.e. that they are different fields.
That was wrong. `potentialInit` does not solve anything - it only declares
that the fields in `0/` came from potentialFoam (`coupledFoam.C`, "fresh
start (potentialInit: fields from potentialFoam expected)"). Both runs used
the same field, and `fieldMinMax` on it gives max|U| = 961.04 m/s at
(-0.510, -0.693, 0.250), p = 0 everywhere. The 679 m/s was the max|U| of a
REJECTED trial step reported on a sentinel rollback line, not of the start
field.

The measurement adds something important instead: **p is exactly 0 in the
whole start field.** potentialFoam without `-writep` writes U only. So the
clip's "keep the total head" had no head to keep - it put isolated spikes
on a flat p = 0.

---

## 5. DEFECT 3: `Uref boundary` is unusable on a case with slip/inletOutlet

Default `Uref boundary` on this case yields **Uref = 321.2 m/s** against a
free stream of 50 (`log.coupledFoam.Uref321_abgebrochen`).

`lineSearch.C` takes the maximum of `mag(U)` over all non-coupled patches.
On this case every patch with a prescribed value is at or below 50 m/s -
inlet 50, road 50, wheels `rotatingWallVelocity` 151.4876 rad/s x 0.323 m =
48.9, everything else `noSlip`. The only two patches that can exceed it are
`tunnelWall` (`slip`) and `outlet` (`inletOutlet`), and both inherit their
values from the interior - here from the peaks of the potential start.

This is not cosmetic. `Uref` scales the entire safety machinery:

| mechanism | with Uref 321 | with Uref 50 (correct) |
|---|---|---|
| `remediation.dynamic` at `cU*Uref` | 1284 m/s | 200 m/s |
| `localLimit` step at `fLoc*Uref` | 160 m/s per step | 25 m/s per step |
| `sentinel` at `UFactor*Uref` | 3212 m/s | 500 m/s |
| `sentinel` at `pFactor*pref` | 2.58e6 | 62500 |

With the inflated value the whole apparatus is inert (`nDyn=0` throughout)
and **defect 1 is invisible** - that run sailed past iteration 18 with a
falling residual while carrying a Cp of -180 it never reported. An
inflated reference does not merely mis-scale the diagnostics, it hides the
defects the diagnostics exist to catch.

`F1_READINESS.md` already had "slip tunnelWall ... on the review list to
double-check". This is that item coming due.

Suggested direction: `boundary` mode should consider only patches that
*prescribe* a value, or fall back when the boundary maximum exceeds the
field maximum by an implausible factor. Note `UrefFallbackFactor` already
implements the mirror-image guard (boundary far BELOW the field) - the
guard for the other direction is missing.

---

## 6. What I changed in the case, and how much to trust it

All of it is in `~/f1/runs/s5_cpl_np60/system/fvSolution`, each entry
commented in place with its measurement. Two are legitimate, two are not:

| change | verdict |
|---|---|
| `Uref 50` (was `boundary`) | **legitimate.** The default was measurably wrong; see defect 3. |
| `ptc.CFLmin 0.05` (was 1) | **legitimate.** The floor was disabling the recovery path; see the table in section 3. |
| `sentinel.UFactor 20` (was 10) | **workaround.** I widened a safety threshold to get past the start transient. Do not keep this. |
| `potentialClip 0`, `potentialInit no` | **avoidance.** I routed around defect 2 rather than fixing it. |

Also in that file, and unrelated to the defects: the convergence stop is
disabled (`forceCoeffsWindow 0`, `residualTol 0`, and
`forceCoeffsRmsWindow 100` because it must be >= 1 and defaults to
`forceCoeffsWindow`) so the run ends at exactly 300 iterations like the
segregated reference. `convergenceMonitor` accepts no "off" mode; that is
the only way to express it.

---

## 7. Suggested order of work

1. **Locate the failing cluster in the mesh.** Are the cells at
   (3.40, 0, 0.43) among the 126 skew faces or not? This single answer
   splits defect 1 into two very different bugs. Cheap: dump the face
   centres of the `skewFaces` set written by `checkMesh`.
2. **Build a small reproducer.** 20.6 M cells and 60 ranks is not a
   debugging loop. The ingredients look like: a symmetry plane, a sharp
   trailing edge on it, skew cells, kOmegaSST. If it reproduces at 100 k
   cells the rest is ordinary work.
3. **Instrument the assembly at those cells** - the Rhie-Chow non-orth
   correction is where both defect 1 and defect 2 surface.
4. Then fix defect 3 (`Uref boundary`) and defect 2 (`potentialClip`),
   which are independent of each other and of defect 1.
5. Only then re-run the F1 case and compare against Cd 1.2385 / Cl -3.4475.

## 8. What was NOT done

- No comparison result exists. coupledFoam has never completed an iteration
  budget on this case.
- The AMR loop was not exercised: the shipped mesh carries no `cellLevel` /
  `pointLevel`, so `refineHexMesh` refuses (`tools/amr_refine.py:340` in the
  cases repository). Enabling it means re-meshing with `./Allrun.remesh`,
  which produces a different mesh and therefore a different comparison.
  Separately, `tools/amr_loop.py` hardcodes `simpleFoam` in two places
  (`start_solver`, and `simpleFoam -postProcess` in `amr_refine.py:382`),
  so it would need a `--solver` flag before it could drive coupledFoam at
  all. coupledFoam itself uses `createMesh.H`, i.e. a static mesh, so live
  `dynamicRefineFvMesh` refinement is not an option for it either.
- Nothing was committed to git; the repository has unrelated uncommitted
  work in it from before today.
