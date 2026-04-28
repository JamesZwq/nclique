"""
Shared parallel-scheduling harness for paper r=1 experiments.

Each experiment script imports from here:
    from bench_lib import (
        ServerConfig, ParallelRunner, parse_phase_timings,
        load_done, raise_stack
    )

Design:
  * ProcessPoolExecutor-style scheduler with memory gating, CPU load gating,
    resume support, and per-job log capture.
  * Each script supplies (a) the job iterator (yielding job dicts) and
    (b) the runner function (takes job dict, returns result dict).
  * Output is one CSV row per job, schema chosen by the script.
"""

from __future__ import annotations
import csv, os, re, signal, socket, subprocess, sys, time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional, Any

# ---------- stack limit ----------
def raise_stack(target_bytes: int = 1 << 30) -> None:
    try:
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
        if hard == resource.RLIM_INFINITY:
            t = target_bytes
        elif hard >= target_bytes:
            t = target_bytes
        else:
            t = hard
        if t != soft:
            try:
                resource.setrlimit(resource.RLIMIT_STACK, (t, hard))
            except (ValueError, OSError):
                resource.setrlimit(resource.RLIMIT_STACK, (max(soft, t - 4096), hard))
        new_soft = resource.getrlimit(resource.RLIMIT_STACK)[0]
        def _fmt(x):
            return "unlimited" if x == resource.RLIM_INFINITY else f"{x/1024/1024:.0f}MB"
        print(f"[stack] {_fmt(soft)} -> {_fmt(new_soft)} (hard={_fmt(hard)})", flush=True)
    except Exception as e:
        print(f"[stack] WARN: {e}", flush=True)


# ---------- per-server config ----------
@dataclass
class ServerConfig:
    name: str
    max_workers: int
    cpu_load_target: Optional[float]   # None disables load gate
    mem_limit_gb: int
    mem_kill_gb: int
    per_proc_mem_gb: int
    datadir: str = "/data/wenqianz"

    @classmethod
    def detect(cls, table: dict[str, "ServerConfig"]) -> "ServerConfig":
        # 1) CLI arg
        if len(sys.argv) > 1 and sys.argv[1] in table:
            return table[sys.argv[1]]
        # 2) Hostname
        host = socket.gethostname().lower()
        for key, cfg in table.items():
            if key in host:
                return cfg
        # 3) Default = first entry
        first_key = next(iter(table))
        print(f"[server] could not detect — defaulting to {first_key}", flush=True)
        return table[first_key]


DEFAULT_SERVERS = {
    "tods1": ServerConfig("tods1", max_workers=60, cpu_load_target=0.85,
                          mem_limit_gb=300, mem_kill_gb=450, per_proc_mem_gb=250),
    "tods2": ServerConfig("tods2", max_workers=80, cpu_load_target=None,
                          mem_limit_gb=300, mem_kill_gb=450, per_proc_mem_gb=250),
}


# ---------- mem / cpu helpers ----------
def get_used_mem_gb() -> float:
    try:
        with open("/proc/meminfo") as f:
            info = {}
            for line in f:
                p = line.split()
                if len(p) >= 2:
                    info[p[0].rstrip(":")] = int(p[1])
        total, avail = info.get("MemTotal", 0), info.get("MemAvailable", 0)
        return (total - avail) / 1024 / 1024
    except Exception:
        return 0.0

def get_proc_rss_gb(pid: int) -> float:
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024 / 1024
    except Exception:
        pass
    return 0.0

def get_load_avg_1min() -> float:
    try:
        return os.getloadavg()[0]
    except Exception:
        return 0.0

def cpu_has_headroom(cfg: ServerConfig, running_count: int) -> bool:
    if running_count >= cfg.max_workers:
        return False
    if cfg.cpu_load_target is None:
        return True
    import multiprocessing as mp
    target = mp.cpu_count() * cfg.cpu_load_target
    return get_load_avg_1min() < target


# ---------- output parsers ----------
_RE_TOTAL = re.compile(r'NucleusCoreDecomposition took:\s*([\d.]+)')
# Peel timing labels emitted by the binary, in priority order:
#   * "ST_V[23] r=1 (peel) took: <ms>"
#   * "ST r=1 took: <ms>"            (immutable-tree ST: peel + nCr-delta)
#   * "Peeling time: <ms>"            (mutable-tree REF and r>=3 variants)
_RE_PEEL_V23  = re.compile(r'ST_V[23]\s+r=1\s+\(peel\)\s+took:\s*([\d.]+)\s*ms')
_RE_PEEL_ST   = re.compile(r'ST\s+r=1\s+took:\s*([\d.]+)\s*ms')
_RE_PEEL_GEN  = re.compile(r'Peeling time:\s*([\d.]+)\s*ms')
_RE_BUILD = re.compile(r'(?:SDCT_Fused|SDCT_MaxClique|ST_V[23]\s+Build)[^\d]*([\d.]+)\s*ms')
_RE_MEM   = re.compile(r'\[Memory-\w+\]\s*Final Memory:\s*([\d.]+)\s*kB')
_RE_LOCAL_ITER = re.compile(r'Local H-index V[2-4][^\d]*converged in (\d+) iterations')

def parse_phase_timings(txt: str) -> dict[str, float]:
    """Pull total / peel / build / mem / iter from binary stdout. -1.0 if absent."""
    out = {"total_ms": -1.0, "peel_ms": -1.0, "build_ms": -1.0,
           "mem_kB": -1.0, "iter_count": -1}
    if (m := _RE_TOTAL.search(txt)): out["total_ms"] = float(m.group(1))
    # peel: try V2/V3 explicit label first, fall back to ST, then generic.
    if (m := _RE_PEEL_V23.search(txt)): out["peel_ms"] = float(m.group(1))
    elif (m := _RE_PEEL_ST.search(txt)):  out["peel_ms"] = float(m.group(1))
    elif (m := _RE_PEEL_GEN.search(txt)): out["peel_ms"] = float(m.group(1))
    if (m := _RE_BUILD.search(txt)): out["build_ms"] = float(m.group(1))
    if (m := _RE_MEM.search(txt)):   out["mem_kB"]   = float(m.group(1))
    if (m := _RE_LOCAL_ITER.search(txt)): out["iter_count"] = int(m.group(1))
    return out


# ---------- CSV resume ----------
def load_done(csv_path: Path, key_fields: tuple[str, ...]) -> set[tuple]:
    """Return {tuple(row[k] for k in key_fields)} for rows where status==OK."""
    done = set()
    if not csv_path.exists():
        return done
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            if row.get("status") != "OK":
                continue
            try:
                done.add(tuple(row[k] for k in key_fields))
            except KeyError:
                pass
    return done


# ---------- parallel runner ----------
@dataclass
class Job:
    key: tuple                                # (graph, s, algo, ...) — composes done-set key
    cmd: list[str]                            # subprocess argv
    env: dict[str, str]                       # env overrides (added on top of os.environ)
    log_path: Optional[Path] = None           # where to write stdout/stderr
    timeout: int = 3600
    extra: dict[str, Any] = None              # passed through to result row


class ParallelRunner:
    """Lightweight scheduler: poll loop with mem/cpu gates."""
    def __init__(self, cfg: ServerConfig, csv_path: Path, fieldnames: list[str],
                 settle_sec: float = 0.1, poll_sec: float = 3.0):
        self.cfg = cfg
        self.csv_path = csv_path
        self.fieldnames = fieldnames
        self.settle_sec = settle_sec
        self.poll_sec = poll_sec
        self.running: list[tuple] = []   # (proc, job, t_start)
        self.shutdown = False

        signal.signal(signal.SIGINT,  self._handle_signal)
        signal.signal(signal.SIGTERM, self._handle_signal)

        # init csv
        if not self.csv_path.exists() or self.csv_path.stat().st_size == 0:
            self.csv_path.parent.mkdir(parents=True, exist_ok=True)
            with self.csv_path.open("w", newline="") as f:
                csv.DictWriter(f, fieldnames=fieldnames).writeheader()

    def _handle_signal(self, sig, frame):
        print("\n[shutdown] killing children...", flush=True)
        self.shutdown = True
        for proc, *_ in self.running:
            try: proc.kill()
            except: pass

    def append_row(self, row: dict):
        with self.csv_path.open("a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=self.fieldnames)
            w.writerow({k: row.get(k, "") for k in self.fieldnames})

    def _launch(self, job: Job) -> None:
        env = os.environ.copy(); env.update(job.env)
        log_f = job.log_path.open("w") if job.log_path else subprocess.PIPE
        proc = subprocess.Popen(
            job.cmd, env=env,
            stdout=log_f, stderr=subprocess.STDOUT,
            preexec_fn=os.setsid,
        )
        proc._log_handle = log_f if job.log_path else None
        self.running.append((proc, job, time.time()))

    def _reap(self, on_finish: Callable[[Job, str, str, dict], None]) -> None:
        """on_finish(job, status, log_text, parsed) appends rows."""
        still = []
        for proc, job, t0 in self.running:
            if proc.poll() is None:
                still.append((proc, job, t0)); continue
            ret = proc.returncode
            # close log handle
            if hasattr(proc, "_log_handle") and proc._log_handle:
                try: proc._log_handle.close()
                except: pass
            # read log
            log_text = ""
            if job.log_path and job.log_path.exists():
                log_text = job.log_path.read_text(errors="replace")
            else:
                try:
                    out, _ = proc.communicate(timeout=5)
                    log_text = out.decode(errors="replace") if isinstance(out, bytes) else (out or "")
                except: pass

            if ret == -9 or ret == 137:           status = "OOM"
            elif ret == 0:                         status = "OK"
            elif ret in (-11, -10, -8, -6, -4):    status = f"SIGNAL({ret})"
            else:                                  status = f"ERROR({ret})"
            wall_ms = (time.time() - t0) * 1000.0
            parsed = parse_phase_timings(log_text)
            parsed["wall_ms"] = wall_ms
            on_finish(job, status, log_text, parsed)
        self.running = still

    def _check_timeouts(self, on_finish):
        now = time.time()
        for i in range(len(self.running) - 1, -1, -1):
            proc, job, t0 = self.running[i]
            if now - t0 > job.timeout:
                try:
                    proc.kill(); proc.wait(timeout=10)
                except: pass
                wall_ms = (time.time() - t0) * 1000.0
                if hasattr(proc, "_log_handle") and proc._log_handle:
                    try: proc._log_handle.close()
                    except: pass
                log_text = job.log_path.read_text(errors="replace") if job.log_path and job.log_path.exists() else ""
                parsed = parse_phase_timings(log_text)
                parsed["wall_ms"] = wall_ms
                on_finish(job, "TIMEOUT", log_text, parsed)
                self.running.pop(i)

    def _check_per_proc_mem(self, on_finish):
        for i in range(len(self.running) - 1, -1, -1):
            proc, job, t0 = self.running[i]
            rss = get_proc_rss_gb(proc.pid)
            if rss > self.cfg.per_proc_mem_gb:
                try:
                    proc.kill(); proc.wait(timeout=10)
                except: pass
                if hasattr(proc, "_log_handle") and proc._log_handle:
                    try: proc._log_handle.close()
                    except: pass
                wall_ms = (time.time() - t0) * 1000.0
                parsed = {"total_ms": -1.0, "peel_ms": -1.0, "build_ms": -1.0,
                          "mem_kB": -1.0, "iter_count": -1, "wall_ms": wall_ms}
                on_finish(job, "OOM", f"OOM RSS={rss:.0f}GB", parsed)
                self.running.pop(i)

    def _kill_newest(self, on_finish):
        if not self.running: return
        proc, job, t0 = self.running.pop()
        try:
            proc.kill(); proc.wait(timeout=10)
        except: pass
        if hasattr(proc, "_log_handle") and proc._log_handle:
            try: proc._log_handle.close()
            except: pass
        wall_ms = (time.time() - t0) * 1000.0
        parsed = {"total_ms": -1.0, "peel_ms": -1.0, "build_ms": -1.0,
                  "mem_kB": -1.0, "iter_count": -1, "wall_ms": wall_ms}
        on_finish(job, "OOM_KILL_NEWEST", "kill newest for OOM", parsed)

    def run(self, jobs: Iterable[Job],
            on_finish: Callable[[Job, str, str, dict], None],
            announce_done: Callable[[Job, str, dict], None] = None) -> None:
        """Main loop. on_finish fires for every completed job (OK / TIMEOUT / OOM /...)."""
        def _wrapped(job, status, log_text, parsed):
            on_finish(job, status, log_text, parsed)
            if announce_done:
                announce_done(job, status, parsed)

        job_iter = iter(jobs)
        exhausted = False

        while not self.shutdown:
            # 1) reap finished
            self._reap(_wrapped)
            self._check_timeouts(_wrapped)
            self._check_per_proc_mem(_wrapped)

            # 2) memory protection
            mem = get_used_mem_gb()
            if mem > self.cfg.mem_kill_gb and self.running:
                print(f"  [mem-emergency] {mem:.0f}GB > kill {self.cfg.mem_kill_gb}GB -> kill newest", flush=True)
                self._kill_newest(_wrapped)
                continue

            # 3) launch new
            if not exhausted and cpu_has_headroom(self.cfg, len(self.running)) and mem < self.cfg.mem_limit_gb:
                try:
                    job = next(job_iter)
                    self._launch(job)
                    time.sleep(self.settle_sec)
                    continue
                except StopIteration:
                    exhausted = True

            # 4) done?
            if exhausted and not self.running:
                break

            # 5) idle
            time.sleep(self.poll_sec)


# ---------- graph helpers ----------
def link_graphs(graphs: list[str], cfg: ServerConfig) -> list[str]:
    os.makedirs("graphs", exist_ok=True)
    missing = []
    for g in graphs:
        f = f"graphs/{g}.edges"
        src = f"{cfg.datadir}/{g}.edges"
        if os.path.exists(f):
            continue
        if os.path.exists(src):
            try:
                os.symlink(src, f)
            except FileExistsError:
                pass
        else:
            missing.append(g)
    if missing:
        print(f"[graphs] missing on {cfg.name}: {missing}", flush=True)
    return [g for g in graphs if os.path.exists(f"graphs/{g}.edges")]


def build_main(cmake_targets: list[str] = None) -> None:
    """Run cmake build with -j 12 cap (project rule)."""
    print("[build] cmake configure...", flush=True)
    subprocess.run("cmake -S . -B build -DCMAKE_BUILD_TYPE=Release".split(),
                   capture_output=True)
    if cmake_targets is None:
        cmake_targets = ["degeneracy_cliques"]
    cmd = ["cmake", "--build", "build", "-j", "12"]
    for t in cmake_targets:
        cmd += ["--target", t]
    print(f"[build] {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[build] FAILED:\n{r.stderr[-1000:]}", flush=True)
        sys.exit(1)
    print("[build] OK", flush=True)
