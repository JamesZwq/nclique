#!/usr/bin/env python3
"""
External bibliography verifier for vldbNuclearR1/.

For every cite key actually USED in the paper, looks up the
authoritative metadata via:

  1. Crossref DOI resolution (if the bib entry has a `doi` field), or
  2. DBLP title search (fallback for entries without a DOI).

Compares author surnames, year, title similarity, and venue against
the bib body and prints a per-entry verdict.

Network results are cached to tools/verify_bib_cache.json so that
re-runs are instant and offline.

Usage:
    python3 vldbNuclearR1/tools/verify_bib.py             # full verify
    python3 vldbNuclearR1/tools/verify_bib.py --no-net    # cache only
    python3 vldbNuclearR1/tools/verify_bib.py --clear     # wipe cache

Exit code is 0 if every used entry is either OK or SKIP (anonymous
/ no-DOI-and-no-DBLP-match); non-zero if any entry has a hard FAIL
(wrong author surname, wrong year, or title similarity < 0.55).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from difflib import SequenceMatcher
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIB = ROOT / "references.bib"
TEX_FILES = [ROOT / "main.tex"] + sorted((ROOT / "sections").glob("*.tex"))
CACHE = Path(__file__).resolve().parent / "verify_bib_cache.json"

USER_AGENT = (
    "vldbNuclearR1-bib-verifier/1.0 "
    "(mailto:zhangwenqian6915@gmail.com)"
)

TITLE_SIM_FAIL = 0.55
TITLE_SIM_WARN = 0.80


# ---------------------------------------------------------------------- bib
def parse_bib(path: Path):
    """Same brace-counting parser as check_bib.py."""
    text = path.read_text()
    out = []
    i, n = 0, len(text)
    while i < n:
        m = re.search(r"@(\w+)\s*\{\s*([^,\s]+)\s*,", text[i:])
        if not m:
            break
        etype = m.group(1).lower()
        key = m.group(2)
        body_start = i + m.end()
        depth = 1
        j = body_start
        while j < n and depth > 0:
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
            j += 1
        body = text[body_start : j - 1]
        fields = {}
        k = 0
        while k < len(body):
            fm = re.search(r"(\w+)\s*=\s*", body[k:])
            if not fm:
                break
            fname = fm.group(1).lower()
            val_start = k + fm.end()
            if val_start >= len(body):
                break
            ch = body[val_start]
            if ch == "{":
                depth2 = 1
                p = val_start + 1
                while p < len(body) and depth2 > 0:
                    if body[p] == "{":
                        depth2 += 1
                    elif body[p] == "}":
                        depth2 -= 1
                    p += 1
                fields[fname] = body[val_start + 1 : p - 1]
                k = p
            elif ch == '"':
                p = val_start + 1
                while p < len(body) and body[p] != '"':
                    p += 1
                fields[fname] = body[val_start + 1 : p]
                k = p + 1
            else:
                p = val_start
                while p < len(body) and body[p] not in ",\n":
                    p += 1
                fields[fname] = body[val_start:p].strip()
                k = p
            while k < len(body) and body[k] in ", \n\r\t":
                k += 1
        out.append((key, etype, fields))
        i = j
    return out


def collect_cites():
    cites = set()
    cite_re = re.compile(r"\\(?:cite|citep|citet|citeauthor|citeyear)\*?\{([^}]+)\}")
    for f in TEX_FILES:
        if not f.exists():
            continue
        for m in cite_re.finditer(f.read_text()):
            for key in m.group(1).split(","):
                key = key.strip()
                if key:
                    cites.add(key)
    return cites


# ---------------------------------------------------------------------- net
def http_json(url: str, timeout: float = 15.0):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                                "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def crossref_lookup(doi: str):
    doi = doi.strip().lower()
    doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", doi)
    if not doi:
        return None
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi, safe='/')}"
    try:
        j = http_json(url)
    except Exception as e:
        return {"_error": f"crossref: {e}"}
    msg = j.get("message", {})
    authors = []
    for a in msg.get("author", []) or []:
        family = a.get("family", "")
        given = a.get("given", "")
        if family:
            authors.append((given, family))
    title = (msg.get("title") or [""])[0]
    container = (msg.get("container-title") or [""])[0]
    year = None
    for k in ("published-print", "published-online", "published",
              "issued", "created"):
        dp = msg.get(k, {}).get("date-parts", [])
        if dp and dp[0]:
            year = str(dp[0][0])
            break
    return {"source": "crossref", "doi": doi, "authors": authors,
            "title": title, "venue": container, "year": year}


def dblp_lookup(title: str, year: str | None):
    q = re.sub(r"[^\w\s]", " ", title).strip()
    if not q:
        return None
    url = ("https://dblp.org/search/publ/api?h=5&format=json&q=" +
           urllib.parse.quote(q))
    try:
        j = http_json(url)
    except Exception as e:
        return {"_error": f"dblp: {e}"}
    hits = j.get("result", {}).get("hits", {}).get("hit", []) or []
    if not hits:
        return None
    target = normalize_title(title)
    best = None
    best_score = 0.0
    for h in hits:
        info = h.get("info", {})
        cand_title = info.get("title", "")
        score = SequenceMatcher(None, target, normalize_title(cand_title)).ratio()
        # Heavily penalise candidates whose year is far from the bib year
        # to avoid "similar-titled but different paper" false matches.
        if year:
            try:
                ydiff = abs(int(info.get("year", "0")) - int(year))
                if ydiff > 3:
                    score -= 0.40
                elif ydiff <= 1:
                    score += 0.05
            except ValueError:
                pass
        if score > best_score:
            best_score = score
            best = info
    if not best or best_score < 0.70:
        return None
    a_field = best.get("authors", {}).get("author", [])
    if isinstance(a_field, dict):
        a_field = [a_field]
    authors = []
    for a in a_field:
        if isinstance(a, dict):
            full = a.get("text", "")
        else:
            full = str(a)
        # DBLP names: "Linyuan Lu 0002" or "Tao Zhou"
        full = re.sub(r"\s+\d{4}$", "", full).strip()
        parts = full.split()
        if parts:
            authors.append((" ".join(parts[:-1]), parts[-1]))
    return {"source": "dblp", "score": round(best_score, 3),
            "authors": authors,
            "title": best.get("title", ""),
            "venue": best.get("venue", ""),
            "year": best.get("year", "")}


# ---------------------------------------------------------- normalization
_ACCENT_MAP = str.maketrans({
    "á": "a", "à": "a", "â": "a", "ä": "a", "ã": "a", "å": "a",
    "é": "e", "è": "e", "ê": "e", "ë": "e",
    "í": "i", "ì": "i", "î": "i", "ï": "i",
    "ó": "o", "ò": "o", "ô": "o", "ö": "o", "õ": "o", "ø": "o",
    "ú": "u", "ù": "u", "û": "u", "ü": "u",
    "ñ": "n", "ç": "c", "ş": "s", "ğ": "g", "ı": "i", "ş": "s",
    "ý": "y", "ÿ": "y",
    "Á": "a", "À": "a", "Â": "a", "Ä": "a", "Ã": "a",
    "É": "e", "È": "e", "Ê": "e", "Ë": "e",
    "Í": "i", "Ì": "i", "Î": "i", "Ï": "i",
    "Ó": "o", "Ò": "o", "Ô": "o", "Ö": "o",
    "Ú": "u", "Ù": "u", "Û": "u", "Ü": "u",
    "Ñ": "n", "Ç": "c", "Ş": "s", "Ğ": "g", "İ": "i",
    "Ý": "y",
})


def normalize_title(t: str) -> str:
    if not t:
        return ""
    t = t.lower()
    # Strip embedded MathML (Crossref sometimes returns it).
    t = re.sub(r"<\s*mml:[^>]*>", " ", t)
    t = re.sub(r"<\s*/\s*mml:[^>]*>", " ", t)
    t = re.sub(r"<[^>]+>", " ", t)
    # Strip TeX math markers.
    t = re.sub(r"\$[^$]*\$", lambda m: m.group(0)[1:-1], t)
    # Strip TeX accent commands like \'{i}, \"u, \v{c}, \u{g} etc.
    t = re.sub(r"\\[a-zA-Z'`\"^~=]+\{?([a-zA-Z])\}?", r"\1", t)
    t = re.sub(r"[\{\}\\]", "", t)
    # Strip Unicode accents.
    t = t.translate(_ACCENT_MAP)
    t = re.sub(r"[^\w\s]", " ", t)
    t = re.sub(r"\s+", " ", t).strip()
    stop = {"a", "an", "the", "of", "for", "in", "on", "and", "with",
            "to", "by", "via"}
    return " ".join(w for w in t.split() if w not in stop)


def normalize_surname(s: str) -> str:
    s = s.lower()
    # TeX dotless-letter macros (\i, \j) keep the letter.
    s = re.sub(r"\\([ij])\b", r"\1", s)
    # TeX accent commands without curly braces: \'a -> a, \"u -> u.
    s = re.sub(r"\\['\"`^~=]([a-z])", r"\1", s)
    # TeX accent commands with curly braces: \'{a} -> a, \v{c} -> c.
    s = re.sub(r"\\[a-z'\"`^~=]+\{([a-z])\}", r"\1", s)
    # Strip remaining TeX commands (\textbf, \textit, etc.).
    s = re.sub(r"\\[a-z]+", "", s)
    # Strip stray backslashes, braces, quotes.
    s = re.sub(r"[\\\{\}\"']", "", s)
    # Unicode accents.
    s = s.translate(_ACCENT_MAP)
    s = re.sub(r"[^a-z]", "", s)
    return s.strip()


def extract_surnames(author_field: str) -> list[str]:
    if not author_field:
        return []
    out = []
    # BibTeX " and " separator is case-insensitive in practice; some
    # entries use "AND".
    for a in re.split(r"\s+and\s+", author_field, flags=re.IGNORECASE):
        a = a.strip()
        if not a:
            continue
        if "," in a:
            sur = a.split(",")[0].strip()
        else:
            sur = a.split()[-1] if a.split() else ""
        sur = normalize_surname(sur)
        if sur:
            out.append(sur)
    return out


# ---------------------------------------------------------------- verify
def verdict(bib_fields: dict, auth: dict | None):
    """Return (status, [issue_lines])."""
    if auth is None:
        return ("SKIP", ["no DOI and no plausible DBLP match"])
    if "_error" in auth:
        return ("SKIP", [auth["_error"]])

    issues = []
    # year (tolerate +/-2 for IEEE/ACM early-access vs journal-issue year)
    by = (bib_fields.get("year") or "").strip()
    ay = (auth.get("year") or "").strip()
    if by and ay and by != ay:
        try:
            if abs(int(by) - int(ay)) > 2:
                issues.append(f"year: bib={by} authoritative={ay} (gap>2)")
        except ValueError:
            issues.append(f"year: bib={by} authoritative={ay}")
    # surnames
    bib_sur = extract_surnames(bib_fields.get("author", ""))
    auth_sur = [normalize_surname(s) for _, s in auth.get("authors", [])]
    auth_sur = [s for s in auth_sur if s]
    if bib_sur and auth_sur:
        # need substring overlap; allow accent / spelling variants
        def matches(a, b):
            return a == b or a[:5] == b[:5] or a in b or b in a
        if not any(matches(bib_sur[0], a) for a in auth_sur[:3]):
            issues.append(
                f"first author: bib={bib_sur[0]!r} authoritative={auth_sur[:3]}"
            )
        else:
            # check author count plausibility
            if len(bib_sur) < len(auth_sur) - 1:
                issues.append(
                    f"author count: bib={len(bib_sur)} authoritative={len(auth_sur)}"
                )
    # title similarity
    bt = normalize_title(bib_fields.get("title", ""))
    at = normalize_title(auth.get("title", ""))
    sim = SequenceMatcher(None, bt, at).ratio() if (bt and at) else 0
    if bt and at:
        if sim < TITLE_SIM_FAIL:
            issues.append(
                f"title sim={sim:.2f}: bib={bt[:60]!r} authoritative={at[:60]!r}"
            )
        elif sim < TITLE_SIM_WARN:
            issues.append(
                f"title sim={sim:.2f} (warning, may be rephrased): "
                f"bib={bt[:60]!r} authoritative={at[:60]!r}"
            )
    # venue
    bv = normalize_title(
        bib_fields.get("journal") or bib_fields.get("booktitle") or ""
    )
    av = normalize_title(auth.get("venue", ""))
    if bv and av:
        v_sim = SequenceMatcher(None, bv, av).ratio()
        # accept any short-substring overlap (e.g. "pvldb" inside "pvldbendowment")
        bv_short = "".join(w[0] for w in bv.split())
        av_short = "".join(w[0] for w in av.split())
        if (v_sim < 0.40 and bv_short != av_short and
                bv not in av and av not in bv):
            issues.append(f"venue: bib={bv!r} authoritative={av!r}")

    if any("first author" in i or "year" in i or "title sim" in i
           and "warning" not in i for i in issues):
        # downgrade title-only warnings to WARN; only hard issues fail
        hard = [i for i in issues if (
            i.startswith("first author") or i.startswith("year:") or
            (i.startswith("title sim") and "warning" not in i)
        )]
        if hard:
            return ("FAIL", issues)
    if issues:
        return ("WARN", issues)
    return ("OK", [])


# ----------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-net", action="store_true", help="cache only")
    ap.add_argument("--clear", action="store_true", help="wipe cache")
    ap.add_argument("--delay", type=float, default=1.0,
                    help="seconds between live HTTP requests")
    ap.add_argument("--only", default=None,
                    help="only verify keys matching this regex")
    args = ap.parse_args()

    if args.clear and CACHE.exists():
        CACHE.unlink()

    cache = {}
    if CACHE.exists():
        try:
            cache = json.loads(CACHE.read_text())
        except Exception:
            cache = {}

    bib = {k: (etype, fields) for k, etype, fields in parse_bib(BIB)}
    used = sorted(collect_cites() & set(bib))
    if args.only:
        used = [k for k in used if re.search(args.only, k)]

    print(f"Verifying {len(used)} used cite keys ...\n")

    results = {}  # key -> (status, issues, auth)
    n_live = 0
    for k in used:
        etype, fields = bib[k]
        # short-circuit: anonymous / under-review
        au = fields.get("author", "").lower()
        note = fields.get("note", "").lower()
        if "anonymous" in au or "under review" in note:
            results[k] = ("SKIP", ["anonymous / under review"], None)
            continue

        auth = cache.get(k)
        if auth is None and not args.no_net:
            # try crossref first
            doi = fields.get("doi", "")
            if doi:
                auth = crossref_lookup(doi)
                n_live += 1
                time.sleep(args.delay)
            if not auth or auth.get("_error"):
                title = fields.get("title", "")
                year = fields.get("year", "")
                if title:
                    auth = dblp_lookup(title, year) or auth
                    n_live += 1
                    time.sleep(args.delay)
            cache[k] = auth
            # persist incrementally
            CACHE.write_text(json.dumps(cache, indent=2,
                                        ensure_ascii=False))

        status, issues = verdict(fields, auth)
        results[k] = (status, issues, auth)

    if not args.no_net and n_live:
        print(f"({n_live} live HTTP lookups; cache written to {CACHE.name})\n")

    counts = {"OK": 0, "WARN": 0, "FAIL": 0, "SKIP": 0}
    for k in used:
        status, issues, auth = results[k]
        counts[status] += 1
        if status == "OK":
            continue
        print(f"[{status}] {k}")
        for line in issues:
            print(f"   {line}")

    print()
    print("Summary:")
    for s in ("OK", "WARN", "FAIL", "SKIP"):
        print(f"  {s:4s}: {counts[s]}")

    return 1 if counts["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
