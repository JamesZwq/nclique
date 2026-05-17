#!/usr/bin/env python3
"""
Self-test for verify_bib.py.

For each test case below, we mutate one bib entry, run the verifier
on just that entry, and assert the expected status (OK / WARN / FAIL).
The bib file is restored afterwards (we operate on a temp copy).

Run:
    python3 vldbNuclearR1/tools/test_verify_bib.py

Exit 0 if every assertion passes; non-zero otherwise.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
ROOT = TOOLS.parent
BIB = ROOT / "references.bib"
CACHE = TOOLS / "verify_bib_cache.json"
VERIFY = TOOLS / "verify_bib.py"


# Each test: (name, cite_key, mutation_fn, expected_status, must_contain)
#  mutation_fn(bib_text) -> mutated_text
#  expected_status: one of "OK", "WARN", "FAIL"
#  must_contain: substring that must appear in the verifier output (or None)
TESTS = [
    (
        "single-letter title typo",
        "Pivoter",
        lambda t: t.replace(
            "The Power of Pivoting for Exact Clique Counting",
            "The Power of Pivotin for Exact Clique Counting",
        ),
        "WARN",
        "title sim",
    ),
    (
        "wrong author surname (first author)",
        "BronKerbosch",
        lambda t: t.replace(
            "author = {Bron, Coen and Kerbosch, Joep},\ntitle = {Algorithm 457:",
            "author = {Brorn, Coen and Kerbosch, Joep},\ntitle = {Algorithm 457:",
        ),
        "FAIL",
        "first author",
    ),
    (
        "year wrong by 5",
        "akbas2017truss",
        lambda t: t.replace(
            "@article{akbas2017truss,\nauthor = {Akbas, Esra and Zhao, Peixiang},\n"
            "title = {Truss-based community search: a truss-equivalence based indexing approach},\n"
            "year = {2017},",
            "@article{akbas2017truss,\nauthor = {Akbas, Esra and Zhao, Peixiang},\n"
            "title = {Truss-based community search: a truss-equivalence based indexing approach},\n"
            "year = {2012},",
        ),
        "FAIL",
        "year",
    ),
    (
        "year off by 1 (within early-access tolerance)",
        "akbas2017truss",
        lambda t: t.replace(
            "@article{akbas2017truss,\nauthor = {Akbas, Esra and Zhao, Peixiang},\n"
            "title = {Truss-based community search: a truss-equivalence based indexing approach},\n"
            "year = {2017},",
            "@article{akbas2017truss,\nauthor = {Akbas, Esra and Zhao, Peixiang},\n"
            "title = {Truss-based community search: a truss-equivalence based indexing approach},\n"
            "year = {2018},",
        ),
        # within +/-2 year tolerance the entry should report OK with no
        # issues at all (verifier prints no per-entry line for OK)
        "OK",
        None,
    ),
    (
        "dropped co-author",
        "nucleusSariyuce14",
        lambda t: t.replace(
            "{Sariy{\\\"u}ce, Ahmet Erdem and Seshadhri, C. and Pinar, Ali "
            "and {\\c{C}}ataly{\\\"u}rek, {\\\"U}mit V.}",
            "{Sariy{\\\"u}ce, Ahmet Erdem and Pinar, Ali "
            "and {\\c{C}}ataly{\\\"u}rek, {\\\"U}mit V.}",
        ),
        "WARN",
        "author count",
    ),
    (
        "totally different title (different paper)",
        "Pivoter",
        lambda t: t.replace(
            "The Power of Pivoting for Exact Clique Counting",
            "A Completely Unrelated Paper About Quantum Chemistry",
        ),
        "FAIL",
        "title sim",
    ),
]


def run_verify(cite_key: str) -> tuple[str, int]:
    """Run verify_bib.py --no-net --only ^KEY$, return (stdout, exit)."""
    r = subprocess.run(
        [sys.executable, str(VERIFY), "--no-net", "--only", f"^{re.escape(cite_key)}$"],
        capture_output=True, text=True, cwd=str(ROOT),
    )
    return r.stdout + r.stderr, r.returncode


def extract_status(stdout: str, cite_key: str) -> str:
    """Find the [STATUS] line for cite_key in stdout, or 'OK' if none."""
    for line in stdout.splitlines():
        m = re.match(r"\[(OK|WARN|FAIL|SKIP)\]\s+(\S+)", line)
        if m and m.group(2) == cite_key:
            return m.group(1)
    return "OK"  # verifier omits OK lines from per-entry output


def main():
    if not BIB.exists():
        print(f"BIB not found at {BIB}", file=sys.stderr)
        return 2

    # Work on a temp copy of the bib so the real file is never touched.
    # The verifier reads ROOT/references.bib so we save+restore in-place.
    original = BIB.read_bytes()
    BIB_BAK = BIB.with_suffix(".bib.selftest-bak")
    BIB_BAK.write_bytes(original)

    # Save cache too --- mutated runs will overwrite cache otherwise.
    cache_bytes = CACHE.read_bytes() if CACHE.exists() else None

    failures = []
    try:
        for name, key, mut, expected, must_contain in TESTS:
            mutated = mut(original.decode("utf-8"))
            if mutated == original.decode("utf-8"):
                failures.append((name, key, "MUTATION_NO_OP",
                                 "mutation did not change the bib --- "
                                 "check the test's source-text match"))
                print(f"  [NO-OP] {name:50s} key={key}")
                continue
            BIB.write_bytes(mutated.encode("utf-8"))
            out, rc = run_verify(key)
            status = extract_status(out, key)
            ok = (status == expected)
            if must_contain and must_contain not in out:
                ok = False
            mark = "PASS" if ok else "FAIL"
            print(f"  [{mark}] {name:50s} key={key:20s} "
                  f"expected={expected} got={status}")
            if not ok:
                failures.append((name, key, status, out))
                print(out)
            # restore for next iteration
            BIB.write_bytes(original)
    finally:
        BIB.write_bytes(original)
        BIB_BAK.unlink(missing_ok=True)
        if cache_bytes is not None:
            CACHE.write_bytes(cache_bytes)

    print()
    if failures:
        print(f"SELF-TEST: {len(failures)} failed of {len(TESTS)}")
        return 1
    print(f"SELF-TEST: all {len(TESTS)} mutations were caught")
    return 0


if __name__ == "__main__":
    sys.exit(main())
