"""Minimal reader of OpenFOAM dictionaries (no OpenFOAM environment needed).

Enough for the harness checks of the case templates (fvSolution,
relaxation.*): nested sub-dictionaries, `key value;` entries (a value may
span several tokens and parenthesised lists), quoted keys, comments.
Directives (#include, #sinclude, ...) are skipped, $-macros are kept
verbatim. Values are strings; normalise() maps numbers and switches to
comparable values. Not a general parser: codeStreams, #calc and regex
semantics are not supported.
"""

from __future__ import annotations

import re
from pathlib import Path

_TOKEN = re.compile(r'"[^"]*"|[{}();]|[^\s{}();"]+')


def _strip(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    text = re.sub(r"//[^\n]*", " ", text)
    # directives: one line each
    return "\n".join(line for line in text.splitlines()
                     if not line.lstrip().startswith("#"))


def parse(text: str) -> dict:
    toks = _TOKEN.findall(_strip(text))
    pos = 0

    def block() -> dict:
        nonlocal pos
        out: dict = {}
        while pos < len(toks):
            t = toks[pos]
            if t == "}":
                pos += 1
                return out
            if t == ";":
                pos += 1
                continue
            key = t.strip('"')
            pos += 1
            if pos < len(toks) and toks[pos] == "{":
                pos += 1
                out[key] = block()
                continue
            val, depth = [], 0
            while pos < len(toks):
                v = toks[pos]
                if v == "(":
                    depth += 1
                elif v == ")":
                    depth -= 1
                elif v == ";" and depth <= 0:
                    pos += 1
                    break
                elif v == "}" and depth <= 0:
                    break
                val.append(v)
                pos += 1
            out[key] = " ".join(val)
        return out

    return block()


def read(path: str | Path) -> dict:
    return parse(Path(path).read_text(errors="replace"))


_SWITCH = {"yes": True, "on": True, "true": True, "y": True,
           "no": False, "off": False, "false": False, "n": False}


def normalise(v):
    """Comparable value: bools for switches, floats for numbers, nested
    dicts recursively, other strings stripped of quotes."""
    if isinstance(v, dict):
        return {k: normalise(x) for k, x in v.items()}
    s = str(v).strip().strip('"')
    if s.lower() in _SWITCH:
        return _SWITCH[s.lower()]
    try:
        return float(s)
    except ValueError:
        return s


def get(d: dict, dotted: str, default=None):
    cur = d
    for p in dotted.split("."):
        if not isinstance(cur, dict) or p not in cur:
            return default
        cur = cur[p]
    return cur


def set_dotted(d: dict, dotted: str, value) -> None:
    """Set a dotted entry, creating missing sub-dictionaries (as
    cflib.case.set_entry does with foamDictionary)."""
    parts = dotted.split(".")
    cur = d
    for p in parts[:-1]:
        if not isinstance(cur.get(p), dict):
            cur[p] = {}
        cur = cur[p]
    cur[parts[-1]] = value
