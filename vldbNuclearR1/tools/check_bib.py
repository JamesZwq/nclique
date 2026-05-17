#!/usr/bin/env python3
"""
Bibliography linter for vldbNuclearR1/.

Run from anywhere:
    python3 vldbNuclearR1/scripts/check_bib.py

Checks performed (each is a separate pass):

  1. BROKEN cite key used in .tex but no entry in references.bib
  2. ORPHAN entry defined but never cited (informational)
  3. DUPLICATE same paper under two keys (detected by DOI, by
     normalized-title similarity, or by author+year+venue match)
  4. MISSING required field for the entry type
  5. PLACEHOLDER or TODO-style entries
  6. KEY/AUTHOR mismatch heuristic --- the cite key implies a first
     author (e.g. "cohen2008trusses") that does not appear in the bib
     body's author list
  7. KEY/TITLE topic mismatch heuristic --- the cite key implies a
     topic word (e.g. "truss") that does not appear in the bib title

Exit code is 0 if no BROKEN / MISSING / PLACEHOLDER / KEY-MISMATCH
issues are found; non-zero otherwise.  Orphans are not failures.
"""

import re
import sys
from collections import defaultdict
from difflib import SequenceMatcher
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIB = ROOT / "references.bib"
TEX_FILES = [ROOT / "main.tex"] + sorted((ROOT / "sections").glob("*.tex"))

REQUIRED_FIELDS = {
    "article":       {"author", "title", "journal", "year"},
    "inproceedings": {"author", "title", "booktitle", "year"},
    "conference":    {"author", "title", "booktitle", "year"},
    "book":          {"author", "title", "publisher", "year"},
    "incollection":  {"author", "title", "booktitle", "publisher", "year"},
    "phdthesis":     {"author", "title", "school", "year"},
    "techreport":    {"author", "title", "institution", "year"},
    "misc":          {"author", "title"},
    "unpublished":   {"author", "title", "note"},
    "manual":        {"title"},
}

PLACEHOLDER_PATTERNS = [
    r"\bPlaceholder\b",
    r"\bTODO\b",
    r"\btbd\b",
    r"\bXXX\b",
    r"\bFIXME\b",
    r"\bAuthor,\s*Author\b",
]


def parse_bib(path: Path):
    """Return [(key, etype, {field: value}, raw_block)] for each entry."""
    text = path.read_text()
    entries = []
    i = 0
    n = len(text)
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
        raw_block = text[i + m.start() : j]
        fields = {}
        # parse "field = {value}" or "field = \"value\""
        k = 0
        while k < len(body):
            fm = re.search(r"(\w+)\s*=\s*", body[k:])
            if not fm:
                break
            fname = fm.group(1).lower()
            val_start = k + fm.end()
            if val_start >= len(body):
                break
            if body[val_start] == "{":
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
            elif body[val_start] == '"':
                p = val_start + 1
                while p < len(body) and body[p] != '"':
                    p += 1
                fields[fname] = body[val_start + 1 : p]
                k = p + 1
            else:
                # bare token
                p = val_start
                while p < len(body) and body[p] not in ",\n":
                    p += 1
                fields[fname] = body[val_start:p].strip()
                k = p
            # skip past comma
            while k < len(body) and body[k] in ", \n\r\t":
                k += 1
        entries.append((key, etype, fields, raw_block))
        i = j
    return entries


def collect_cites():
    """Return {cite_key: [(file, line), ...]} across all .tex files."""
    out = defaultdict(list)
    cite_re = re.compile(r"\\(?:cite|citep|citet|citeauthor|citeyear)\*?\{([^}]+)\}")
    for f in TEX_FILES:
        if not f.exists():
            continue
        for lineno, line in enumerate(f.read_text().splitlines(), 1):
            for m in cite_re.finditer(line):
                for key in m.group(1).split(","):
                    key = key.strip()
                    if key:
                        out[key].append((f.name, lineno))
    return out


def normalize_title(t: str) -> str:
    t = re.sub(r"\{|\}|\\[A-Za-z]+|[\W_]+", " ", t).lower()
    t = re.sub(r"\s+", " ", t).strip()
    # drop common filler words for fuzzy match
    filler = {"a", "an", "the", "of", "for", "in", "on", "and", "with", "to"}
    return " ".join(w for w in t.split() if w not in filler)


def normalize_doi(d: str) -> str:
    d = d.strip().lower()
    d = re.sub(r"^https?://(dx\.)?doi\.org/", "", d)
    return d


def first_author_surname(author_field: str) -> str | None:
    if not author_field:
        return None
    first = author_field.split(" and ")[0].strip()
    if "," in first:
        return first.split(",")[0].strip().lower()
    parts = first.split()
    return parts[-1].lower() if parts else None


def first_year_in_key(key: str) -> str | None:
    m = re.search(r"(19|20)\d{2}", key)
    return m.group(0) if m else None


def key_alpha_prefix(key: str) -> str:
    """Letters of the key before the first digit, lowercased."""
    return re.match(r"^[A-Za-z\-_]*", key).group(0).strip("-_").lower()


def main():
    print(f"Linting {BIB}")
    print(f"Across {len(TEX_FILES)} .tex files: "
          f"{', '.join(p.name for p in TEX_FILES)}\n")

    entries = parse_bib(BIB)
    by_key = {k: (etype, fields, raw) for k, etype, fields, raw in entries}
    cite_uses = collect_cites()

    fail = False
    sep = lambda title: print(f"\n{'='*5} {title} {'='*5}")

    # 1. Broken cite keys ----------------------------------------------------
    sep("1. BROKEN cite keys (used but not defined)")
    broken = sorted(set(cite_uses) - set(by_key))
    if broken:
        fail = True
        for k in broken:
            usages = ", ".join(f"{f}:{ln}" for f, ln in cite_uses[k][:3])
            print(f"  {k:35s}  used at {usages}")
    else:
        print("  (none)")

    # 2. Orphan entries ------------------------------------------------------
    sep("2. ORPHAN entries (defined but never cited)  [informational]")
    orphans = sorted(set(by_key) - set(cite_uses))
    if orphans:
        print(f"  {len(orphans)} unused entries (top 10 shown):")
        for k in orphans[:10]:
            print(f"  {k}")
    else:
        print("  (none)")

    # 3. Duplicates: by DOI / normalized title / author+year+venue -----------
    sep("3. DUPLICATE entries (same paper under multiple keys)")
    used_keys = [k for k in by_key if k in cite_uses]
    dups_by_doi = defaultdict(list)
    dups_by_title = defaultdict(list)
    for k in used_keys:
        _, fields, _ = by_key[k]
        if doi := fields.get("doi"):
            dups_by_doi[normalize_doi(doi)].append(k)
        if title := fields.get("title"):
            nt = normalize_title(title)
            if len(nt) > 12:  # avoid 1-word collisions
                dups_by_title[nt].append(k)
    seen_pairs = set()
    found_dups = False
    for k_doi, keys in dups_by_doi.items():
        if len(keys) > 1:
            found_dups = True
            pair = tuple(sorted(keys))
            seen_pairs.add(pair)
            print(f"  DOI {k_doi}: {', '.join(keys)}")
    for nt, keys in dups_by_title.items():
        if len(keys) > 1 and tuple(sorted(keys)) not in seen_pairs:
            found_dups = True
            print(f"  title \"{nt[:60]}...\": {', '.join(keys)}")
    # also fuzzy author+year+venue check among used entries
    fingerprint = defaultdict(list)
    for k in used_keys:
        _, fields, _ = by_key[k]
        au = first_author_surname(fields.get("author", ""))
        yr = fields.get("year", "").strip()
        venue = (fields.get("booktitle") or fields.get("journal") or "").lower()
        venue_short = re.sub(r"[^a-z]", "", venue)[:15]
        if au and yr:
            fp = (au, yr, venue_short)
            fingerprint[fp].append(k)
    for fp, keys in fingerprint.items():
        pair = tuple(sorted(keys))
        if len(keys) > 1 and pair not in seen_pairs:
            found_dups = True
            seen_pairs.add(pair)
            print(f"  author={fp[0]}, year={fp[1]}, venue~={fp[2]}: {', '.join(keys)}")
    if found_dups:
        fail = True
    else:
        print("  (none)")

    # 4. Missing required fields --------------------------------------------
    sep("4. MISSING required fields")
    any_missing = False
    for k in used_keys:
        etype, fields, _ = by_key[k]
        req = REQUIRED_FIELDS.get(etype, set())
        missing = [f for f in req if not fields.get(f, "").strip()]
        if missing:
            any_missing = True
            fail = True
            print(f"  {k}  (@{etype}): missing {', '.join(missing)}")
    if not any_missing:
        print("  (none)")

    # 5. Placeholder / TODO entries -----------------------------------------
    sep("5. PLACEHOLDER / TODO-style entries")
    any_ph = False
    for k in used_keys:
        _, _, raw = by_key[k]
        for pat in PLACEHOLDER_PATTERNS:
            if re.search(pat, raw, re.IGNORECASE):
                any_ph = True
                fail = True
                print(f"  {k}: matched pattern /{pat}/")
                break
    if not any_ph:
        print("  (none)")

    # 6. Key/author mismatch ------------------------------------------------
    # Heuristic: only flag SURNAME-STYLE keys, i.e. /^[a-z]+\d{4}\w*$/.
    # Descriptive keys like "Pivoter", "communityFinder", "kcore_3",
    # "communities_dec1" are intentional and skipped.
    sep("6. KEY/AUTHOR mismatch (surname-style keys only)")
    surname_style = re.compile(r"^([a-z]{4,})(\d{4})[a-z]*$")
    any_mismatch = False
    for k in used_keys:
        m = surname_style.match(k)
        if not m:
            continue
        prefix, year = m.group(1), m.group(2)
        _, fields, _ = by_key[k]
        au = first_author_surname(fields.get("author", ""))
        bib_year = fields.get("year", "").strip()
        # author check
        if au and not (prefix.startswith(au[:4]) or au[:4] in prefix):
            any_mismatch = True
            fail = True
            print(f"  key={k!r:35s} prefix \"{prefix}\" mismatches bib author \"{au}\"")
            continue
        # year check
        if bib_year and bib_year != year:
            any_mismatch = True
            fail = True
            print(f"  key={k!r:35s} year {year} mismatches bib year {bib_year}")
    if not any_mismatch:
        print("  (none)")

    # 7. Key/title topic mismatch ------------------------------------------
    # Heuristic: same surname-style key set, then check that any trailing
    # topic word (after the year) appears in the title.
    sep("7. KEY/TITLE topic mismatch (surname-style keys only)")
    surname_topic = re.compile(r"^[a-z]+\d{4}([a-z]+)$")
    any_tmismatch = False
    for k in used_keys:
        m = surname_topic.match(k)
        if not m:
            continue
        topic = m.group(1)
        if len(topic) < 4:
            continue  # too short to be a content word
        _, fields, _ = by_key[k]
        raw_title = fields.get("title", "").lower()
        # Normalize: drop math markup ($...$), TeX braces, hyphens, dots.
        title = re.sub(r"\$[^$]*\$", lambda m: m.group(0)[1:-1], raw_title)
        title = re.sub(r"[\{\}\\]", "", title)
        title = re.sub(r"[-.]", "", title)
        title = re.sub(r"\s+", " ", title)
        # accept if topic appears anywhere in title (substring is fine)
        if topic in title:
            continue
        # accept simple morphology: drop trailing 's' or 'es'
        if topic.rstrip("s") in title or topic.rstrip("es") in title:
            continue
        any_tmismatch = True
        fail = True
        print(f"  key={k!r:35s} suggests topic \"{topic}\" but title is "
              f"\"{raw_title[:60]}...\"")
    if not any_tmismatch:
        print("  (none)")

    # 8. Anonymous / "under review" entries (informational) -----------------
    sep("8. ANONYMOUS / under-review entries  [informational, expected for "
        "blind submissions]")
    for k in used_keys:
        _, fields, _ = by_key[k]
        au = fields.get("author", "")
        note = fields.get("note", "").lower()
        if "anonymous" in au.lower() or "under review" in note:
            print(f"  {k}: author={au!r}, note={fields.get('note','')!r}")

    print()
    if fail:
        print("RESULT: FAIL --- see issues above")
        return 1
    print("RESULT: OK --- no broken/missing/duplicate/mismatch issues")
    return 0


if __name__ == "__main__":
    sys.exit(main())
