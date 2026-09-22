"""Source gates of amendment D (D3, D4; DECISIONS.md D-064). No build needed.

- D3: `grep -nE '\\b(SMALL|VSMALL|ROOTVSMALL|GREAT|VGREAT)\\b'` over src/
  (and applications/coupledFoam) returns nothing: guards use the typed
  constants of src/blockMatrix/coupledConstants.H. applications/test/ is
  excluded (spec D3).
- D4: no native OpenFOAM reduction of field data (gSum, gSumMag, gAverage,
  sum(), sumMag(), ...) in the code (comments stripped); every accumulation
  goes through doubleReduce (double accumulator, MPI_DOUBLE). Calls
  qualified with doubleReduce:: or blockKernels:: are allowed; the
  definitions in doubleReduce*.{H,C} and blockKernels.H are exempt.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCAN = [REPO / "src", REPO / "applications" / "coupledFoam"]
SUFFIXES = {".C", ".H"}

D3_PATTERN = re.compile(r"\b(SMALL|VSMALL|ROOTVSMALL|GREAT|VGREAT)\b")

D4_NAMES = (
    "gSum", "gSumMag", "gSumSqr", "gSumProd", "gSumCmptMag", "gSumCmptProd",
    "gAverage", "gWeightedSum", "gWeightedAverage", "domainIntegrate",
    "weightedAverage", "sumProd", "sumCmptMag", "sumCmptProd",
    "sum", "sumMag", "sumSqr", "average",
)
D4_PATTERN = re.compile(
    r"(?P<pre>[\w:]*::|\.|->)?\b(?P<name>" + "|".join(D4_NAMES) + r")\s*\("
)
D4_ALLOWED_QUALIFIERS = ("doubleReduce::", "blockKernels::")
D4_EXEMPT = {"doubleReduce.H", "doubleReduce.C", "doubleReduceTemplates.H",
             "blockKernels.H"}


def _sources():
    for root in SCAN:
        for p in sorted(root.rglob("*")):
            if p.suffix in SUFFIXES and "lnInclude" not in p.parts \
                    and "Make" not in p.parts:
                yield p


def _strip_comments(text: str) -> str:
    """Remove /* */ and // comments (keeps line numbers)."""
    def block(m):
        return "\n" * m.group(0).count("\n")
    text = re.sub(r"/\*.*?\*/", block, text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def test_gate_D3_guard_constants():
    hits = []
    for p in _sources():
        for i, line in enumerate(p.read_text(errors="replace").splitlines(), 1):
            if D3_PATTERN.search(line):
                hits.append(f"{p.relative_to(REPO)}:{i}: {line.strip()}")
    assert not hits, "D3 gate: precision-dependent guard constants\n" \
        + "\n".join(hits)


def test_gate_D4_native_reductions():
    hits = []
    for p in _sources():
        if p.name in D4_EXEMPT:
            continue
        code = _strip_comments(p.read_text(errors="replace"))
        for i, line in enumerate(code.splitlines(), 1):
            for m in D4_PATTERN.finditer(line):
                pre = m.group("pre") or ""
                if pre in (".", "->"):
                    continue            # member function of another class
                if pre.endswith(D4_ALLOWED_QUALIFIERS):
                    continue
                hits.append(f"{p.relative_to(REPO)}:{i}: {line.strip()}")
    assert not hits, "D4 gate: native reduction of field data\n" \
        + "\n".join(hits)
