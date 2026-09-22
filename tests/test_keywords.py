"""Keyword reference completeness (D-066).

Every default constant in src/include/coupledDefaults.H must be documented
in docs/KEYWORDS.md (as `name` in backticks: in the constant column of a
keyword row, or in the compile-time table). A constant added without an
entry fails this test. No OpenFOAM environment needed.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULTS = REPO / "src" / "include" / "coupledDefaults.H"
KEYWORDS = REPO / "docs" / "KEYWORDS.md"

# constexpr <type> <name> = ...;   (types may contain spaces, *, ::)
_CONST = re.compile(r"^\s*constexpr\s+[\w:\s\*]+?\b(\w+)\s*=", re.M)


def _constants() -> list[str]:
    text = DEFAULTS.read_text()
    # drop comments so that commented-out constants do not count
    text = re.sub(r"//[^\n]*", "", text)
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return _CONST.findall(text)


def test_defaults_parsed():
    names = _constants()
    # sanity: the header has well over 100 constants
    assert len(names) > 100, names
    assert len(names) == len(set(names)), "duplicate constant names"


def test_every_constant_documented():
    doc = KEYWORDS.read_text()
    missing = [n for n in _constants() if f"`{n}`" not in doc]
    assert not missing, (
        "constants in coupledDefaults.H without an entry in docs/KEYWORDS.md "
        f"(add a keyword row or a compile-time row): {missing}")
