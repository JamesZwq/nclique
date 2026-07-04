#!/usr/bin/env python3
"""
v5 zero-change certificate — decisive kill/confirm experiment.

Claim (relaxation lemma): after inserting e=(u,v), if EVERY seed
x in {u,v} ∪ W (W = common neighbors) has
    rd'(x) := #{ s-cliques K ∋ x in G'  with all other members ranked
                 AFTER x in the base peel order O(G-e) }   <=  core_{G-e}(x),
then NO core changes.  rd'(x) is capped at core(x)+1 (only the <= decision
matters), so each seed costs one small capped clique count over x's
higher-ranked neighborhood.

Per sampled edge e of G:
  base = G - e; run V3 on base dumping core[] + pop-rank[] (original-id
  space via PIVOTER_DUMP_MAPPING); the C++ helper `v5cert` computes the
  certificate; ground-truth change = (core(base) != core(G)) per vertex.

Metrics:
  absorption      = fraction of edges where the certificate FIRES
  soundness       = among fired edges, fraction that truly did NOT change
                    (MUST be 100%; a single violation kills the lemma)
  zero-coverage   = among no-change edges, fraction the certificate fires
  seed Δ stats    = per-edge #seeds and max (rd'(x)-core(x)) clamp

Usage:
  python3 bench_v5_certificate.py --graph soc-Epinions1.edges --s 3 \
      --edges 300 --workers 6 --tmpdir <scratch>
"""
from __future__ import annotations
import argparse, os, random, subprocess, tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from bench_dynamic_locality import read_edge_file, run_v3

CERT = os.environ.get("V5CERT", "./build/bin/v5cert")

_W = {}


def _init(graph_path, s, workdir, full_cores):
    n, lines = read_edge_file(graph_path)
    _W.update(n=n, lines=lines, s=s, workdir=workdir, full=full_cores)


def run_v3_order(graph_path, s, tag, workdir):
    """V3 -> (core dict, rank dict) in original-id space."""
    cf = os.path.join(workdir, f"c_{tag}.tsv")
    rf = os.path.join(workdir, f"r_{tag}.tsv")
    mf = os.path.join(workdir, f"m_{tag}.tsv")
    env = os.environ.copy()
    env.update({"PIVOTER_RUN_ST_V3": "1", "PIVOTER_DUMP_CORE": cf,
                "PIVOTER_DUMP_ORDER": rf, "PIVOTER_DUMP_MAPPING": mf,
                "OMP_NUM_THREADS": "1"})
    p = subprocess.run(["./build/bin/degeneracy_cliques", graph_path, "1", str(s)],
                       env=env, capture_output=True, text=True, timeout=1800)
    if p.returncode != 0:
        raise RuntimeError(p.stderr[-400:])
    new2orig = {}
    for l in open(mf):
        if l[0] == "#":
            continue
        a, b = l.split(); new2orig[int(a)] = int(b)
    core, rank = {}, {}
    for l in open(cf):
        if l[0] == "#":
            continue
        a, b = l.split(); core[new2orig[int(a)]] = float(b)
    for l in open(rf):
        if l[0] == "#":
            continue
        a, b = l.split(); rank[new2orig[int(a)]] = int(b)
    for f in (cf, rf, mf):
        os.unlink(f)
    return core, rank


def one_edge(idx):
    lines, n, s, wd = _W["lines"], _W["n"], _W["s"], _W["workdir"]
    p = lines[idx].split(); u, v = int(p[0]), int(p[1])
    tag = f"e{idx}_{os.getpid()}"
    gf = os.path.join(wd, f"g_{tag}.edges")
    with open(gf, "w") as f:
        f.write(f"{n} {len(lines)-1}\n")
        for i, ln in enumerate(lines):
            if i != idx:
                f.write(ln + "\n")
    try:
        base_core, base_rank = run_v3_order(gf, s, tag, wd)
    finally:
        os.unlink(gf)

    # ground truth: did any core change?  (base vs full G)
    full = _W["full"]
    changed = 0
    for x in base_core.keys() | full.keys():
        if base_core.get(x, 0.0) != full.get(x, 0.0):
            changed += 1

    # write base graph edges + core + rank for the C++ certificate helper
    cf = os.path.join(wd, f"cc_{tag}.tsv")
    with open(cf, "w") as f:
        for x, c in base_core.items():
            f.write(f"{x}\t{c:.0f}\t{base_rank.get(x, -1)}\n")
    try:
        # v5cert loads base graph (G-e), inserts (u,v), tests the certificate.
        # It re-reads the base edge file; reuse gf via a fresh write.
        with open(gf, "w") as f:
            f.write(f"{n} {len(lines)-1}\n")
            for i, ln in enumerate(lines):
                if i != idx:
                    f.write(ln + "\n")
        r = subprocess.run([CERT, gf, str(s), cf, str(u), str(v)],
                           capture_output=True, text=True, timeout=600)
    finally:
        for f in (gf, cf):
            if os.path.exists(f):
                os.unlink(f)
    if r.returncode != 0:
        return dict(u=u, v=v, error=r.stderr[-200:])
    fires = None
    nseeds = maxover = cert_us = -1
    for line in r.stdout.splitlines():
        if line.startswith("CERT"):
            kv = dict(t.split("=") for t in line.split()[1:])
            fires = kv["fires"] == "1"
            nseeds = int(kv["nseeds"]); maxover = int(kv["maxover"])
            cert_us = float(kv["cert_us"])
    return dict(u=u, v=v, changed=changed, fires=int(fires),
                nseeds=nseeds, maxover=maxover, cert_us=cert_us)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--edges", type=int, default=300)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tmpdir", default=None)
    args = ap.parse_args()
    gname = Path(args.graph).stem
    wd = args.tmpdir or tempfile.mkdtemp(prefix="v5c_")
    Path(wd).mkdir(parents=True, exist_ok=True)
    n, lines = read_edge_file(args.graph)
    print(f"[{gname} s={args.s}] n={n} m={len(lines)}", flush=True)
    full = run_v3(args.graph, args.s, "full", wd)
    rng = random.Random(args.seed)
    picks = rng.sample(range(len(lines)), args.edges)
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_init,
                             initargs=(args.graph, args.s, wd, full)) as ex:
        for i, r in enumerate(ex.map(one_edge, picks)):
            rows.append(r)
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(picks)}", flush=True)

    ok = [r for r in rows if "error" not in r]
    errs = len(rows) - len(ok)
    fired = [r for r in ok if r["fires"]]
    nochange = [r for r in ok if r["changed"] == 0]
    change = [r for r in ok if r["changed"] > 0]
    # soundness: fired but actually changed = LEMMA VIOLATION
    violations = [r for r in fired if r["changed"] > 0]
    cov = sum(1 for r in nochange if r["fires"])
    ct = sorted(r["cert_us"] for r in ok)
    def pct(a, q): return a[min(len(a)-1, int(q*len(a)))] if a else 0
    print(f"\n===== V5 CERTIFICATE  {gname} s={args.s}  ({len(ok)} edges, {errs} err) =====")
    print(f"  no-change edges: {len(nochange)}   change edges: {len(change)}")
    print(f"  ABSORPTION (cert fires): {len(fired)}/{len(ok)} "
          f"({100.0*len(fired)/max(1,len(ok)):.1f}%)")
    print(f"  SOUNDNESS (fired & truly unchanged): "
          f"{len(fired)-len(violations)}/{len(fired)} "
          f"-- VIOLATIONS={len(violations)} (MUST be 0)")
    if violations:
        for r in violations[:5]:
            print(f"    !! VIOLATION u={r['u']} v={r['v']} changed={r['changed']} "
                  f"maxover={r['maxover']}")
    print(f"  ZERO-COVERAGE (fires | no-change): {cov}/{len(nochange)} "
          f"({100.0*cov/max(1,len(nochange)):.1f}%)")
    print(f"  cert cost us: median={pct(ct,0.5):.0f} p90={pct(ct,0.9):.0f} "
          f"max={ct[-1] if ct else 0:.0f}")
    seeds = sorted(r["nseeds"] for r in ok)
    print(f"  #seeds: median={pct(seeds,0.5)} p90={pct(seeds,0.9)} max={seeds[-1] if seeds else 0}")


if __name__ == "__main__":
    main()
