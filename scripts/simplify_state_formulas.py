#!/usr/bin/env python3
"""Conservative display-only cleanup of TeX coefficients in Dicke expansions.

This is a small structural formatter, not a symbolic algebra engine.  It only
unwraps a complete round coefficient immediately before a Dicke ket, at a
state-term boundary.  Sums and complex coefficients remain grouped.  A leading
minus in a sum may be factored out by reversing its *top-level* signs; signs in
fractions, radicals, indices and exponents are never distributed through.

Import ``simplify_state_formula`` to format one expression.  Executing this
file prints a read-only JSON report; it does not rewrite any catalogue.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any, Iterator


@dataclass(frozen=True)
class _Token:
    start: int
    end: int
    value: str
    kind: str = "atom"
    bracket: str = ""


_SIZERS = {"left", "right", "big", "Big", "bigg", "Bigg", "bigl", "bigr", "Bigl", "Bigr", "biggl", "biggr", "Biggl", "Biggr"}
_SPACES = {"quad", "qquad", "thinspace", "enspace", "hspace", ",", ";", "!", " ", ":"}
_RELATIONS = {"=", "&", r"\propto", r"\simeq", r"\approx", r"\equiv", r"\\"}
_PRODUCTS = {"*", "/", r"\cdot", r"\times", r"\otimes", "^", "_"}
_KET = re.compile(r"(?:(?:\\(?:left|bigl|Bigl|biggl|Biggl)\s*)?(?:\||\\(?:lvert|vert)\b)\s*D\s*[_^]|\\ket\s*\{\s*D|(?<!\\)D\s*_)")


def _tokens(tex: str) -> list[_Token]:
    out: list[_Token] = []
    i = 0
    while i < len(tex):
        start = i
        char = tex[i]
        if char.isspace() or char == "~":
            i += 1
            while i < len(tex) and (tex[i].isspace() or tex[i] == "~"):
                i += 1
            out.append(_Token(start, i, tex[start:i], "space"))
            continue
        if char == "\\":
            i += 1
            if i < len(tex) and tex[i].isalpha():
                while i < len(tex) and tex[i].isalpha():
                    i += 1
            elif i < len(tex):
                i += 1
            command = tex[start + 1:i]
            if command in _SIZERS:
                j = i
                while j < len(tex) and tex[j].isspace():
                    j += 1
                if j < len(tex) and tex[j] in "()[]":
                    char = tex[j]
                    i = j + 1
                    out.append(_Token(start, i, tex[start:i], "open" if char in "([" else "close", "(" if char in "()" else "["))
                    continue
            out.append(_Token(start, i, tex[start:i], "space" if command in _SPACES else "atom"))
            continue
        i += 1
        if char in "({[":
            out.append(_Token(start, i, char, "open", char))
        elif char in ")}]":
            out.append(_Token(start, i, char, "close", {")": "(", "}": "{", "]": "["}[char]))
        else:
            out.append(_Token(start, i, char))
    return out


def _structure(tex: str):
    tokens = _tokens(tex)
    stack: list[int] = []
    pairs: dict[int, int] = {}
    parents: dict[int, int | None] = {}
    for i, token in enumerate(tokens):
        if token.kind == "close":
            if not stack or tokens[stack[-1]].bracket != token.bracket:
                return None
            opening = stack.pop()
            pairs[opening] = i
            parents[i] = stack[-1] if stack else None
        else:
            parents[i] = stack[-1] if stack else None
            if token.kind == "open":
                stack.append(i)
    return None if stack else (tokens, pairs, parents)


def _top_indices(tokens, pairs, parent=None) -> list[int]:
    start = 0 if parent is None else parent + 1
    end = len(tokens) if parent is None else pairs[parent]
    result = []
    i = start
    while i < end:
        if tokens[i].kind != "space":
            result.append(i)
        i = pairs[i] + 1 if tokens[i].kind == "open" else i + 1
    return result


def _top_terms(tex: str):
    """Return signed shallow terms; nested groups are indivisible atoms."""
    parsed = _structure(tex)
    if parsed is None:
        return None
    tokens, pairs, _ = parsed
    top = _top_indices(tokens, pairs)
    if any(tokens[i].value in {r"\pm", r"\mp", r"\begin", r"\end"} for i in top):
        return None
    terms = []
    sign = 1
    start = 0
    have_body = False
    previous = None
    for i in top:
        token = tokens[i]
        if token.value in {"+", "-"} and previous not in _PRODUCTS:
            if have_body:
                body = tex[start:token.start].strip()
                if not body:
                    return None
                terms.append((sign, body))
                sign = -1 if token.value == "-" else 1
                have_body = False
            elif token.value == "-":
                sign *= -1
            start = token.end
        else:
            have_body = True
        previous = token.value
    body = tex[start:].strip()
    if not body or not have_body:
        return None
    terms.append((sign, body))
    return terms


def _skip_spaces(tex: str, start: int) -> int:
    for token in _tokens(tex[start:]):
        if token.kind != "space":
            return start + token.start
    return len(tex)


def _at_term_start(tokens, pairs, parents, opening: int) -> bool:
    parent = parents[opening]
    siblings = _top_indices(tokens, pairs, parent)
    before = [i for i in siblings if i < opening]
    if not before:
        return True
    previous = tokens[before[-1]].value
    return previous in {"+", "-"} | _RELATIONS


def _unwrap_coefficients(tex: str) -> str:
    parsed = _structure(tex)
    if parsed is None:
        return tex
    tokens, pairs, parents = parsed
    edits = []
    for opening, closing in pairs.items():
        left, right = tokens[opening], tokens[closing]
        if left.bracket != "(" or not _at_term_start(tokens, pairs, parents, opening):
            continue
        following = _skip_spaces(tex, right.end)
        if _KET.match(tex, following) is None:
            continue
        inner = tex[left.end:right.start]
        if _KET.search(inner) or any(mark in inner for mark in (r"\rangle", r"\langle", "|", "=", "&")):
            continue
        terms = _top_terms(inner)
        if terms is None:
            continue
        if len(terms) == 1:
            sign, body = terms[0]
            # Do not reintroduce the redundant unit coefficients removed earlier.
            if body == "1":
                body = ""
            replacement = ("-" if sign < 0 else "") + body
        elif terms[0][0] < 0:
            # -a+b becomes -(a-b), without entering a radical or fraction.
            positive = "".join(("" if i == 0 else "+" if sign < 0 else "-") + body for i, (sign, body) in enumerate(terms))
            replacement = r"-\left(" + positive + r"\right)"
        else:
            continue  # A genuine positive-leading sum still needs its grouping.
        if replacement != tex[left.start:right.end]:
            edits.append((left.start, right.end, replacement))
    # Candidates are sibling coefficient atoms; nevertheless reject overlap.
    last_start = len(tex) + 1
    for start, end, replacement in sorted(edits, reverse=True):
        if end <= last_start:
            tex = tex[:start] + replacement + tex[end:]
            last_start = start
    return tex


def _coefficient_signs(tex: str) -> str:
    parsed = _structure(tex)
    if parsed is None:
        return tex
    tokens, pairs, parents = parsed
    edits = {}
    for ket in _KET.finditer(tex):
        containing = next((i for i, token in enumerate(tokens) if token.start <= ket.start() < token.end), None)
        if containing is None:
            continue
        parent = parents[containing]
        if parent is not None:
            outer = _top_indices(tokens, pairs, parents[parent])
            before = [i for i in outer if i < parent]
            if before and tokens[before[-1]].value in {"^", "_"}:
                continue  # Never edit superscript/subscript contents.
        top = [i for i in _top_indices(tokens, pairs, parent) if tokens[i].start < ket.start()]
        i = 0
        while i < len(top):
            j = i
            while j < len(top) and tokens[top[j]].value in {"+", "-"}:
                j += 1
            if j - i >= 2:
                previous = tokens[top[i - 1]].value if i else None
                if previous not in {"^", "_"}:
                    minus = sum(tokens[top[k]].value == "-" for k in range(i, j)) % 2
                    symbol = "-" if minus else "+"
                    if symbol == "+" and (previous is None or previous in _RELATIONS or previous in _PRODUCTS):
                        symbol = ""
                    edits[(tokens[top[i]].start, tokens[top[j - 1]].end)] = symbol
            i = max(j, i + 1)
    for (start, end), replacement in sorted(edits.items(), reverse=True):
        tex = tex[:start] + replacement + tex[end:]
    return tex


def simplify_state_formula(tex: str) -> str:
    """Simplify redundant Dicke-coefficient wrappers without altering data.

    Unbalanced input is returned unchanged.  Only ``|D...>``/``D_k`` state
    expressions are considered; no general algebraic rewriting is attempted.
    """
    if not isinstance(tex, str):
        raise TypeError("A TeX state expression must be a string.")
    if _KET.search(tex) is None or _structure(tex) is None:
        return tex
    # Extra nested scalar wrappers become removable on the following pass.
    for _ in range(len(tex) + 1):
        simplified = _coefficient_signs(_unwrap_coefficients(tex))
        if simplified == tex:
            return tex
        tex = simplified
    raise RuntimeError("State-expression formatting did not converge.")


def _state_fields(value: Any, path: str = "$") -> Iterator[tuple[str, str]]:
    if isinstance(value, dict):
        for key, child in value.items():
            if key in {"stateLatex", "state_latex", "stateTex", "stateFormula"} and isinstance(child, str):
                yield path + "." + key, child
            else:
                yield from _state_fields(child, path + "." + key)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from _state_fields(child, path + f"[{index}]")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    paths = [root / folder / "data/catalogue.json" for folder in ("symmetric-sector-catalogue", "anticoherent-state-catalogue", "ac-subspace-catalogue", "majorana-degeneracy-catalogue")]
    paths += sorted((root / "sector-length-sets/data").glob("*.json"))
    report = []
    for path in paths:
        fields = list(_state_fields(json.loads(path.read_text())))
        changed = [{"path": key, "before": formula, "after": simplify_state_formula(formula)} for key, formula in fields if simplify_state_formula(formula) != formula]
        if fields:
            report.append({"file": str(path.relative_to(root)), "stateFields": len(fields), "changedFields": len(changed), "examples": changed[:3]})
    print(json.dumps({"readOnly": True, "files": report}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
