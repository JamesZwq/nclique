#!/usr/bin/env python3
"""Build and validate a multi-r NSI2 index on ca-GrQc (configurable).

The gate deliberately uses two independent correctness checks:

* exact positive-core distributions from fresh fixed-(r,s) native runs; and
* exact per-r-clique values dumped by the repository's reference decomposition.

All generated indexes, dumps, query files, command logs, and benchmark output
live under --work-dir (under /tmp by default).
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import os
from pathlib import Path
import platform
import random
import re
import shlex
import subprocess
import sys
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
INT_RE = re.compile(r"[+-]?\d+")
CORE_RE = re.compile(r"^core=([+-]?\d+) count=(\d+)$")
PLANE_CELL_RE = re.compile(r"^\[nsi-plane-cell\] r=(\d+) s=(\d+)\b")

Clique = Tuple[int, ...]
Cell = Tuple[int, int]
Distribution = Dict[int, int]


class GateFailure(RuntimeError):
    pass


def clean_environment() -> Dict[str, str]:
    """Remove algorithm-selection state inherited from the caller."""
    env = dict(os.environ)
    for key in list(env):
        if key.startswith("SCT_") or key.startswith("PIVOTER_"):
            del env[key]
    env.pop("MEM_DBG", None)
    env.update({
        "LC_ALL": "C",
        "LANG": "C",
        "OMP_NUM_THREADS": "1",
        "OMP_DYNAMIC": "FALSE",
    })
    return env


def tail(text: str, lines: int = 80) -> str:
    xs = text.splitlines()
    return "\n".join(xs[-lines:])


class Runner:
    def __init__(self, work_dir: Path, timeout: float):
        self.work_dir = work_dir
        self.timeout = timeout
        self.command_log = work_dir / "commands.log"
        self.command_log.write_text(
            f"# nsi2_validate.py\n# host={platform.platform()}\n"
            "# child env: LC_ALL=C LANG=C OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE; "
            "inherited SCT_*/PIVOTER_* and MEM_DBG removed\n",
            encoding="utf-8",
        )

    def run(
        self,
        label: str,
        command: Sequence[object],
        env: Mapping[str, str],
        env_overrides: Optional[Mapping[str, object]] = None,
    ) -> str:
        cmd = [str(x) for x in command]
        shown_env = ""
        if env_overrides:
            shown_env = "env " + " ".join(
                f"{k}={shlex.quote(str(v))}" for k, v in sorted(env_overrides.items())
            ) + " "
        shown = shown_env + shlex.join(cmd)
        print(f"$ {shown}", flush=True)
        with self.command_log.open("a", encoding="utf-8") as log:
            log.write(f"\n# {label}\n$ {shown}\n")
        try:
            proc = subprocess.run(
                cmd,
                cwd=REPO_ROOT,
                env=dict(env),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=self.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise GateFailure(
                f"{label} timed out after {self.timeout:.0f}s: {shown}"
            ) from exc

        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)
        (self.work_dir / f"{safe}.stdout").write_text(proc.stdout, encoding="utf-8")
        (self.work_dir / f"{safe}.stderr").write_text(proc.stderr, encoding="utf-8")
        if proc.returncode != 0:
            raise GateFailure(
                f"{label} failed with exit {proc.returncode}\n"
                f"command: {shown}\n"
                f"--- stdout tail ---\n{tail(proc.stdout)}\n"
                f"--- stderr tail ---\n{tail(proc.stderr)}"
            )
        return proc.stdout


def require_file(path: Path, executable: bool = False) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise GateFailure(f"required file does not exist: {path}")
    if executable and not os.access(path, os.X_OK):
        raise GateFailure(f"required binary is not executable: {path}")
    return path


def expected_cells(r_min: int, r_max: int, s_max: int) -> List[Cell]:
    return [(r, s) for r in range(r_min, r_max + 1)
            for s in range(r + 1, s_max + 1)]


def parse_plane_distributions(stdout: str, cells: Iterable[Cell]) -> Dict[Cell, Distribution]:
    wanted = set(cells)
    result: Dict[Cell, Distribution] = {}
    current: Optional[Cell] = None
    for line_no, line in enumerate(stdout.splitlines(), 1):
        m = PLANE_CELL_RE.match(line)
        if m:
            current = (int(m.group(1)), int(m.group(2)))
            if current in result:
                raise GateFailure(f"duplicate plane cell {current} at stdout line {line_no}")
            result[current] = {}
            continue
        m = CORE_RE.match(line)
        if m and current is not None:
            core, count = int(m.group(1)), int(m.group(2))
            if core in result[current]:
                raise GateFailure(f"duplicate core {core} in plane cell {current}")
            result[current][core] = count
    missing = wanted - set(result)
    extra = set(result) - wanted
    if missing or extra:
        raise GateFailure(f"plane cell set mismatch: missing={sorted(missing)} extra={sorted(extra)}")
    return result


def parse_fixed_distribution(stdout: str) -> Distribution:
    result: Distribution = {}
    for line in stdout.splitlines():
        m = CORE_RE.match(line)
        if not m:
            continue
        core, count = int(m.group(1)), int(m.group(2))
        if core in result:
            raise GateFailure(f"duplicate fixed-run core line for core={core}")
        result[core] = count
    return result


def positive_distribution(dist: Mapping[int, int]) -> Distribution:
    return {core: count for core, count in dist.items() if core != 0 and count != 0}


def compare_distributions(label: str, got: Mapping[int, int], expected: Mapping[int, int]) -> None:
    if dict(got) == dict(expected):
        return
    keys = sorted(set(got) | set(expected))
    diffs = [(k, got.get(k), expected.get(k)) for k in keys if got.get(k) != expected.get(k)]
    raise GateFailure(
        f"{label} distribution mismatch ({len(diffs)} differing levels); first 12: {diffs[:12]}"
    )


def parse_reference_record(line: str, r: int, path: Path, line_no: int) -> Optional[Tuple[Clique, int]]:
    stripped = line.rstrip("\r\n")
    if not stripped or stripped.startswith("#"):
        return None
    fields = stripped.split("\t")
    if len(fields) != 2:
        raise GateFailure(f"{path}:{line_no}: expected one tab separating clique and core")
    vertex_tokens = fields[0].split()
    if len(vertex_tokens) != r or any(INT_RE.fullmatch(x) is None for x in vertex_tokens):
        raise GateFailure(f"{path}:{line_no}: malformed {r}-clique")
    if INT_RE.fullmatch(fields[1].strip()) is None:
        raise GateFailure(f"{path}:{line_no}: core is not an exact integer: {fields[1]!r}")
    clique = tuple(sorted(int(x) for x in vertex_tokens))
    if len(set(clique)) != r:
        raise GateFailure(f"{path}:{line_no}: repeated vertex in clique {clique}")
    core = int(fields[1])
    if core < 0:
        raise GateFailure(f"{path}:{line_no}: negative reference core {core}")
    return clique, core


@dataclass
class ReservoirResult:
    cliques: List[Clique]
    total: int
    zero_total: int
    forced_zero: bool


def reservoir_sample_reference(path: Path, r: int, sample_size: int, seed: int,
                               require_zero: bool) -> ReservoirResult:
    sample_rng = random.Random(seed + 1_000_003 * r)
    zero_rng = random.Random((seed ^ 0x5A17C0DE) + 1_000_033 * r)
    reservoir: List[Clique] = []
    reservoir_core: List[int] = []
    zero_choice: Optional[Clique] = None
    total = zero_total = 0
    with path.open("r", encoding="utf-8") as src:
        for line_no, line in enumerate(src, 1):
            rec = parse_reference_record(line, r, path, line_no)
            if rec is None:
                continue
            clique, core = rec
            total += 1
            if len(reservoir) < sample_size:
                reservoir.append(clique)
                reservoir_core.append(core)
            else:
                j = sample_rng.randrange(total)
                if j < sample_size:
                    reservoir[j] = clique
                    reservoir_core[j] = core
            if core == 0:
                zero_total += 1
                if zero_rng.randrange(zero_total) == 0:
                    zero_choice = clique

    if total == 0:
        raise GateFailure(f"reference dump contains no {r}-cliques: {path}")
    if require_zero and zero_choice is None:
        raise GateFailure(
            f"reference boundary dump has no zero-core {r}-clique; "
            "use --allow-no-zero only for graphs where none exists"
        )
    forced = False
    if zero_choice is not None and not any(core == 0 for core in reservoir_core):
        at = sample_rng.randrange(len(reservoir))
        reservoir[at] = zero_choice
        reservoir_core[at] = 0
        forced = True
    if len(set(reservoir)) != len(reservoir):
        raise GateFailure(f"reference dump yielded duplicate sampled {r}-cliques")
    reservoir.sort()
    return ReservoirResult(reservoir, total, zero_total, forced)


@dataclass
class ReferenceScan:
    answers: Dict[Clique, int]
    distribution: Distribution
    total: int
    zero_total: int


def scan_reference(path: Path, r: int, wanted: Set[Clique]) -> ReferenceScan:
    answers: Dict[Clique, int] = {}
    dist: Counter[int] = Counter()
    total = zero_total = 0
    with path.open("r", encoding="utf-8") as src:
        for line_no, line in enumerate(src, 1):
            rec = parse_reference_record(line, r, path, line_no)
            if rec is None:
                continue
            clique, core = rec
            total += 1
            if core == 0:
                zero_total += 1
            else:
                dist[core] += 1
            if clique in wanted:
                if clique in answers:
                    raise GateFailure(f"duplicate sampled clique in {path}: {clique}")
                answers[clique] = core
    missing = wanted - set(answers)
    if missing:
        raise GateFailure(f"{path} is missing {len(missing)} sampled cliques; first: {sorted(missing)[:3]}")
    return ReferenceScan(answers, dict(dist), total, zero_total)


def parse_integer_matrix(stdout: str, rows: int, columns: int, label: str) -> List[List[int]]:
    matrix: List[List[int]] = []
    for line_no, line in enumerate(stdout.splitlines(), 1):
        if not line.strip():
            continue
        tokens = line.split()
        if len(tokens) != columns or any(INT_RE.fullmatch(x) is None for x in tokens):
            raise GateFailure(
                f"{label}: malformed answer line {line_no}; expected {columns} exact integers: {line!r}"
            )
        matrix.append([int(x) for x in tokens])
    if len(matrix) != rows:
        raise GateFailure(f"{label}: answer count {len(matrix)} != expected {rows}")
    return matrix


def compare_answer_matrix(label: str, cliques: Sequence[Clique], got: Sequence[Sequence[int]],
                          expected: Sequence[Sequence[int]], s_values: Sequence[int]) -> None:
    if got == expected:
        return
    mismatches = []
    for i, (g_row, e_row) in enumerate(zip(got, expected)):
        for j, (g, e) in enumerate(zip(g_row, e_row)):
            if g != e:
                mismatches.append((cliques[i], s_values[j], g, e))
                if len(mismatches) == 10:
                    break
        if len(mismatches) == 10:
            break
    raise GateFailure(f"{label}: bit-exact mismatch; first entries (clique,s,got,ref): {mismatches}")


def write_clique_file(path: Path, cliques: Sequence[Clique]) -> None:
    with path.open("w", encoding="utf-8") as out:
        for clique in cliques:
            out.write(" ".join(str(v) for v in clique) + "\n")


def write_benchmark_file(path: Path, samples: Mapping[int, Sequence[Clique]], s_max: int,
                         seed: int) -> int:
    # One point cell per distinct sampled clique.  The workload is balanced by
    # r and cycles uniformly through that column's admissible s values.  A
    # round-robin cell prefix makes nsi_query's first 1000 cold samples balanced.
    strata: Dict[Cell, List[Tuple[int, int, Clique]]] = {}
    rng = random.Random(seed ^ 0x4E534932)
    for r, cliques in sorted(samples.items()):
        width = s_max - r
        for j, clique in enumerate(cliques):
            s = r + 1 + (j % width)
            strata.setdefault((r, s), []).append((r, s, clique))
    for records in strata.values():
        rng.shuffle(records)

    keys = sorted(strata)
    positions = {key: 0 for key in keys}
    prefix: List[Tuple[int, int, Clique]] = []
    prefix_size = min(1000, sum(len(x) for x in strata.values()))
    while len(prefix) < prefix_size:
        progressed = False
        for key in keys:
            at = positions[key]
            if at < len(strata[key]) and len(prefix) < prefix_size:
                prefix.append(strata[key][at])
                positions[key] = at + 1
                progressed = True
        if not progressed:
            break
    remainder: List[Tuple[int, int, Clique]] = []
    for key in keys:
        remainder.extend(strata[key][positions[key]:])
    rng.shuffle(remainder)
    records = prefix + remainder
    with path.open("w", encoding="utf-8") as out:
        out.write(f"# seed={seed} schema: r s v1 ... vr\n")
        for r, s, clique in records:
            out.write(f"{r} {s} " + " ".join(str(v) for v in clique) + "\n")
    return len(records)


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Build and bit-exactly validate a serialized multi-r NSI2 plane index."
    )
    p.add_argument("--engine", type=Path, default=SCRIPT_DIR / "region_native_sct_peel")
    p.add_argument("--query", type=Path, default=SCRIPT_DIR / "nsi_query")
    p.add_argument("--reference", type=Path, default=REPO_ROOT / "build/bin/degeneracy_cliques")
    p.add_argument("--graph", type=Path, default=REPO_ROOT / "data/ca-GrQc.edges")
    p.add_argument("--r-min", type=int, default=3)
    p.add_argument("--r-max", type=int, default=5)
    p.add_argument("--s-max", type=int, default=6)
    p.add_argument("--sample-per-r", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=20_260_710)
    p.add_argument("--work-dir", type=Path, default=Path("/tmp/pivoter-nsi2-validate"))
    p.add_argument("--index", type=Path, default=None,
                   help="index output (default: WORK_DIR/plane.nsi2)")
    p.add_argument("--mce-budget", type=float, default=120.0)
    p.add_argument("--timeout", type=float, default=600.0,
                   help="timeout in seconds for each child process")
    p.add_argument("--cold-mib", type=int, default=128)
    p.add_argument("--warm-reps", type=int, default=1,
                   help="calls per warm timed interval (default: 1; values >1 measure repeated-query locality)")
    p.add_argument("--reuse-index", action="store_true",
                   help="validate an existing --index instead of rebuilding it")
    p.add_argument("--allow-no-zero", action="store_true",
                   help="do not require each r sample to exercise an unsupported zero-core clique")
    p.add_argument("--skip-benchmark", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = make_parser().parse_args(argv)
    if args.r_min < 1 or args.r_max < args.r_min or args.r_max > 255 or args.s_max < args.r_max + 1:
        raise GateFailure("require 1 <= r-min <= r-max <= 255 and s-max >= r-max+1")
    if args.sample_per_r < 1:
        raise GateFailure("--sample-per-r must be positive")
    if args.timeout <= 0 or args.mce_budget <= 0:
        raise GateFailure("--timeout and --mce-budget must be positive")
    if not 1 <= args.cold_mib <= 4096 or not 1 <= args.warm_reps <= 100_000:
        raise GateFailure("invalid --cold-mib or --warm-reps")

    engine = require_file(args.engine, executable=True)
    query = require_file(args.query, executable=True)
    reference = require_file(args.reference, executable=True)
    graph = require_file(args.graph)
    work_dir = args.work_dir.expanduser().resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    index = (args.index.expanduser().resolve() if args.index else work_dir / "plane.nsi2")
    runner = Runner(work_dir, args.timeout)
    base_env = clean_environment()
    cells = expected_cells(args.r_min, args.r_max, args.s_max)
    print(
        f"NSI2 gate: graph={graph} plane=r{args.r_min}..{args.r_max},s<={args.s_max} "
        f"sample={args.sample_per_r}/r seed={args.seed} work={work_dir}",
        flush=True,
    )

    # 1. Build the serialized plane and retain stdout for the aggregate gate.
    if args.reuse_index:
        index = require_file(index)
        plane_stdout_path = work_dir / "plane.stdout"
        if not plane_stdout_path.is_file():
            raise GateFailure("--reuse-index also requires WORK_DIR/plane.stdout for aggregate validation")
        plane_stdout = plane_stdout_path.read_text(encoding="utf-8")
        print(f"reusing index {index}", flush=True)
    else:
        if index.exists():
            index.unlink()
        plane_overrides = {
            "SCT_RSWEEP": "1",
            "SCT_RMIN": args.r_min,
            "SCT_RMAX": args.r_max,
            "SCT_SMAX": args.s_max,
            "SCT_INDEX_OUT": index,
        }
        plane_env = dict(base_env)
        plane_env.update({k: str(v) for k, v in plane_overrides.items()})
        plane_stdout = runner.run(
            "plane",
            [engine, graph, args.r_min, args.r_min + 1, "--mce-budget", args.mce_budget],
            plane_env,
            plane_overrides,
        )
        if not index.is_file() or index.stat().st_size < 4:
            raise GateFailure(f"plane build did not create index: {index}")
    with index.open("rb") as src:
        if src.read(4) != b"NSI2":
            raise GateFailure(f"index does not have NSI2 magic: {index}")
    plane_dists = parse_plane_distributions(plane_stdout, cells)
    print(f"plane build OK: {index.stat().st_size} serialized bytes, {len(cells)} cells", flush=True)

    # 2. Load-time and precise byte accounting from the real NSI2 loader.
    stats_stdout = runner.run("stats", [query, index, "stats"], base_env)
    if "format=NSI2" not in stats_stdout or "byte-accounting" not in stats_stdout:
        raise GateFailure("NSI2 stats output is missing format or byte-accounting confirmation")
    print("--- NSI2 stats ---")
    print(stats_stdout.rstrip())

    # 3. Fresh standalone fixed-(r,s) aggregate distributions.
    fixed_dists: Dict[Cell, Distribution] = {}
    for r, s in cells:
        stdout = runner.run(f"fixed_r{r}_s{s}", [engine, graph, r, s, "--mce-budget", args.mce_budget], base_env)
        fixed = positive_distribution(parse_fixed_distribution(stdout))
        plane = positive_distribution(plane_dists[(r, s)])
        compare_distributions(f"plane vs fixed r={r},s={s}", plane, fixed)
        fixed_dists[(r, s)] = fixed
        print(
            f"aggregate r={r} s={s}: PASS ({len(fixed)} core levels, "
            f"{sum(fixed.values())} positive-core r-cliques)",
            flush=True,
        )

    # 4. Independent per-clique REF dumps and deterministic reservoir samples.
    dump_paths: Dict[Cell, Path] = {}
    for r, s in cells:
        dump = work_dir / f"ref_r{r}_s{s}.core"
        if dump.exists():
            dump.unlink()
        ref_overrides = {"PIVOTER_RUN_REF": "1", "PIVOTER_DUMP_CORE": dump}
        ref_env = dict(base_env)
        ref_env.update({k: str(v) for k, v in ref_overrides.items()})
        runner.run(
            f"ref_r{r}_s{s}",
            [reference, graph, r, s, "default"],
            ref_env,
            ref_overrides,
        )
        if not dump.is_file() or dump.stat().st_size == 0:
            raise GateFailure(f"reference did not create dump: {dump}")
        dump_paths[(r, s)] = dump

    samples: Dict[int, List[Clique]] = {}
    expected: Dict[Cell, Dict[Clique, int]] = {}
    total_sampled_cells = 0
    for r in range(args.r_min, args.r_max + 1):
        boundary = r + 1
        reservoir = reservoir_sample_reference(
            dump_paths[(r, boundary)], r, args.sample_per_r, args.seed,
            require_zero=not args.allow_no_zero,
        )
        samples[r] = reservoir.cliques
        wanted = set(reservoir.cliques)
        boundary_total: Optional[int] = None
        for s in range(boundary, args.s_max + 1):
            scan = scan_reference(dump_paths[(r, s)], r, wanted)
            if boundary_total is None:
                boundary_total = scan.total
            elif scan.total != boundary_total:
                raise GateFailure(
                    f"reference r-clique universe changed for r={r}: "
                    f"boundary={boundary_total}, s={s} has {scan.total}"
                )
            compare_distributions(f"REF vs fixed r={r},s={s}", scan.distribution, fixed_dists[(r, s)])
            expected[(r, s)] = scan.answers
            total_sampled_cells += len(reservoir.cliques)
        print(
            f"sample r={r}: {len(reservoir.cliques)}/{reservoir.total} cliques, "
            f"boundary-zero={reservoir.zero_total}, forced-zero={reservoir.forced_zero}",
            flush=True,
        )

    # 5. Kernel-only and adjacency-validated point and row APIs.
    point_checks = row_checks = 0
    for r in range(args.r_min, args.r_max + 1):
        cliques = samples[r]
        qfile = work_dir / f"queries_r{r}.txt"
        write_clique_file(qfile, cliques)
        s_values = list(range(r + 1, args.s_max + 1))
        expected_rows = [[expected[(r, s)][clique] for s in s_values] for clique in cliques]

        for j, s in enumerate(s_values):
            want = [[row[j]] for row in expected_rows]
            for validated in (False, True):
                mode = "pointfile-validated" if validated else "pointfile"
                command: List[object] = [query, index, mode]
                if validated:
                    command.append(graph)
                command.extend([r, s, qfile])
                label = f"point_{'validated' if validated else 'kernel'}_r{r}_s{s}"
                stdout = runner.run(label, command, base_env)
                got = parse_integer_matrix(stdout, len(cliques), 1, label)
                compare_answer_matrix(label, cliques, got, want, [s])
                point_checks += len(cliques)
            print(f"point r={r} s={s}: kernel+validated PASS ({len(cliques)} each)", flush=True)

        for validated in (False, True):
            mode = "rowfile-validated" if validated else "rowfile"
            command = [query, index, mode]
            if validated:
                command.append(graph)
            command.extend([r, qfile])
            label = f"row_{'validated' if validated else 'kernel'}_r{r}"
            stdout = runner.run(label, command, base_env)
            got = parse_integer_matrix(stdout, len(cliques), len(s_values), label)
            compare_answer_matrix(label, cliques, got, expected_rows, s_values)
            row_checks += len(cliques) * len(s_values)
        print(f"row r={r}: kernel+validated PASS ({len(cliques)} rows each)", flush=True)

    # 6. Mixed-r latency benchmark.  nsi_query reports warm/cache-cold,
    # point/row, and kernel/validation-inclusive median and p95.
    benchmark_stdout = ""
    benchmark_count = 0
    if not args.skip_benchmark:
        bench_file = work_dir / "benchmark_queries.txt"
        benchmark_count = write_benchmark_file(bench_file, samples, args.s_max, args.seed)
        if benchmark_count < 1000:
            raise GateFailure(
                f"benchmark has only {benchmark_count} queries; increase --sample-per-r or use --skip-benchmark"
            )
        benchmark_stdout = runner.run(
            "benchmark",
            [query, index, "bench", graph, bench_file,
             "--cold-mib", args.cold_mib, "--warm-reps", args.warm_reps],
            base_env,
        )
        print("--- latency benchmark ---")
        print(benchmark_stdout.rstrip())

    distinct_samples = sum(len(x) for x in samples.values())
    print("=== NSI2 VALIDATION PASS ===")
    print(f"aggregate distributions: {len(cells)}/{len(cells)} exact")
    print(
        f"random sample: {distinct_samples} distinct r-cliques; "
        f"{total_sampled_cells} distinct (r,s,R) answers"
    )
    print(f"pointfile answers: {point_checks} exact (kernel + validated)")
    print(f"rowfile cells: {row_checks} exact (kernel + validated)")
    if benchmark_stdout:
        print(f"latency workload: {benchmark_count} mixed-r queries; median+p95 table above")
    print(f"artifacts and exact commands: {work_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except GateFailure as exc:
        print(f"NSI2 GATE FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
