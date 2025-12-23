# -*- coding: utf-8 -*-
"""
Aggregate bazel run logs and make time & memory plots.

Assumptions:
- Logs are in an extracted folder like:
    /Users/zhangwenqian/UNSW/pivoter/python/experiment/data/bazel_runs_20250829_153227
  containing many *.log (or .out/.txt) files produced by:
    nohup /bin/time -v bazel run :NucleusDecomposition_main -- -s -rounds 1 -r <R> -ss <S> ... <dataset.adj> > runRS.log 2>&1 &

- We parse:
   * dataset name from the command line or helper lines
   * r from "-r N"
   * s from "-ss N"
   * total wall time from "/bin/time -v" (Elapsed (wall clock) time ...)
   * Maximum resident set size (kbytes)

Output:
- Per-dataset, three figures per r in {1,2,3}:
    <dataset>_r<r>_time_overall.png   (total wall-clock seconds)
    <dataset>_r<r>_mem.png            (Max RSS in MB)
- An index CSV: _plots_index.csv
- A small _preview.html with inline images.
"""

import os
import re
import tarfile
import matplotlib.pyplot as plt
import pandas as pd

# ---- Helpers to gather candidate log directories ----
def _candidate_dirs(base_dir, base_root_dir):
    """Return a sorted list of folders to scan for logs.
    Includes `base_dir` (if exists) and all subfolders of `base_root_dir`
    whose names start with 'bazel_runs_'.
    """
    dirs = []
    if base_dir and os.path.isdir(base_dir):
        dirs.append(base_dir)
    if base_root_dir and os.path.isdir(base_root_dir):
        # include the root itself for standalone logs
        dirs.append(base_root_dir)
        for name in os.listdir(base_root_dir):
            full = os.path.join(base_root_dir, name)
            if os.path.isdir(full) and name.startswith('bazel_runs_'):
                dirs.append(full)
    return sorted(set(dirs))

# ---- Configuration ----
# On your Mac, set BASE_DIR to the extracted folder and leave BASE_TAR as None.
BASE_ROOT_DIR = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data"
BASE_DIR = None
BASE_TAR = None  # e.g., "/Users/zhangwenqian/Downloads/bazel_runs_20250829_153227.tar.gz"
OUT_DIR  = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/image_bazel"

os.makedirs(OUT_DIR, exist_ok=True)

def parse_elapsed_to_seconds(text):
    s = text.strip()
    m = re.search(r'(\d+):(\d{2}):(\d{2}(?:\.\d+)?)$', s)   # h:mm:ss(.ms)
    if m:
        h = int(m.group(1)); mi = int(m.group(2)); sec = float(m.group(3))
        return h*3600 + mi*60 + sec
    m = re.search(r'(\d+):(\d{2}(?:\.\d+)?)$', s)            # m:ss(.ms)
    if m:
        mi = int(m.group(1)); sec = float(m.group(2))
        return mi*60 + sec
    m = re.search(r'(\d+(?:\.\d+)?)$', s)                     # plain seconds
    if m:
        return float(m.group(1))
    return None

def dataset_from_path(path):
    if not path:
        return None
    return os.path.basename(path.strip().strip('"').strip("'"))

def dataset_nice_name(base):
    if base is None:
        return None
    name = os.path.basename(base)
    known_exts = {".edges", ".adj", ".txt", ".mtx", ".csv", ".tsv", ".bin", ".gz"}
    while True:
        root, ext = os.path.splitext(name)
        if ext and ext.lower() in known_exts and root:
            name = root
            continue
        break
    return name

def scan_filename_for_dataset_r_s(path):
    """Fallback: parse r, s, dataset from filename patterns like
    nucleus_r2_s4_soc-pokec-relationships.log
    or any ...r<r>_s< s or ss><s>_<dataset>.<ext>
    """
    name = os.path.basename(path)
    r_val = None
    s_val = None
    dataset = None
    # Try r<r>_s< s or ss><s>_<dataset>
    m = re.search(r"r(\d+)[-_]s(?:s)?(\d+)[-_](.+?)(?:\.[^.]+)?$", name)
    if m:
        try:
            r_val = int(m.group(1))
            s_val = int(m.group(2))
        except Exception:
            r_val = None; s_val = None
        ds_token = m.group(3)
        # If token already carries an extension we keep it; else assume .adj
        if re.search(r"\.(?:edges|adj)$", ds_token, flags=re.IGNORECASE):
            dataset = ds_token
        else:
            dataset = ds_token + ".adj"
    return dataset, r_val, s_val

def scan_text_for_dataset_r_s(lines):
    r_val = None
    s_val = None
    dataset = None

    # 1) Dataset from explicit "dataset :" line
    for ln in lines:
        m = re.search(r'^\s*dataset\s*:\s*(\S+)', ln)
        if m:
            dataset = dataset_from_path(m.group(1))
            break

    # 2) Dataset from helper line
    if dataset is None:
        for ln in lines:
            m = re.search(r'about to call runAndPrint for dataset\s+(\S+)', ln)
            if m:
                dataset = dataset_from_path(m.group(1))
                break

    # 3) Try to parse from command lines (handle -r, -s/-ss, and dataset *.edges|*.adj anywhere)
    for ln in lines:
        if ('Command being timed' in ln) or ('bazel run' in ln) or (' -ss ' in ln) or (' -s ' in ln) or (' -r ' in ln):
            if r_val is None:
                m = re.search(r'(?:(?:^)|\s)-r\s+(\d+)', ln)
                if m:
                    r_val = int(m.group(1))
            if s_val is None:
                # accept -ss N or -s N
                m = re.search(r'(?:(?:^)|\s)-(?:ss|s)\s+(\d+)', ln)
                if m:
                    s_val = int(m.group(1))
            if dataset is None:
                # accept absolute or relative paths ending with .edges or .adj
                m = re.search(r'(\S+\.(?:edges|adj))', ln)
                if m:
                    dataset = dataset_from_path(m.group(1))

    # 4) Fallback: look for "r : N" and "s : N" summary lines
    if r_val is None:
        for ln in lines:
            m = re.search(r'^\s*r\s*:\s*(\d+)', ln)
            if m:
                r_val = int(m.group(1))
                break
    if s_val is None:
        for ln in lines:
            m = re.search(r'^\s*s\s*:\s*(\d+)', ln)
            if m:
                s_val = int(m.group(1))
                break

    # 5) Final fallback: dataset anywhere in the file
    if dataset is None:
        for ln in lines[::-1]:
            m = re.search(r'(\S+\.(?:edges|adj))', ln)
            if m:
                dataset = dataset_from_path(m.group(1))
                break

    return dataset, r_val, s_val

def scan_text_for_time_and_mem(lines):
     elapsed = None
     maxrss_kb = None
     exit_status = None
     terminated_by_signal = False
     for ln in lines:
         # if 'Elapsed (wall clock) time' in ln:
         #     val = ln.split(':', maxsplit=1)[-1].strip()
         #     e = parse_elapsed_to_seconds(val)
         #     if e is not None:
         #         elapsed = e
         if '### Running Time' in ln:
         #    ### Running Time: 1508.4400189999999
            m = re.search(r'###\s*Running Time\s*:\s*([\d\.]+)', ln)
            if m:
                e = parse_elapsed_to_seconds(m.group(1))
                if e is not None:
                    elapsed = e
         elif 'Maximum resident set size (kbytes)' in ln:
             m = re.search(r'(\d+)', ln)
             if m:
                 maxrss_kb = int(m.group(1))
         elif ln.strip().startswith('Exit status:'):
             m = re.search(r'Exit status:\s*(\d+)', ln)
             if m:
                 exit_status = int(m.group(1))
         # 新增：信号终止检测
         if 'Command terminated by signal' in ln:
             terminated_by_signal = True
     # 只要被信号终止，就强制置为非零
     if terminated_by_signal:
         exit_status = -1
     maxrss_mb = (maxrss_kb / 1024.0) if maxrss_kb is not None else None
     return elapsed, maxrss_mb, exit_status

# ---- Parse thread count from log ----
def scan_text_for_threads(lines):
    """Parse a thread count line like '### Threads: 16' (or 'Threads: 16')."""
    for ln in lines:
        m = re.search(r"(?:###\s*)?Threads:\s*(\d+)", ln)
        if m:
            try:
                return int(m.group(1))
            except Exception:
                pass
    return None

def read_lines(path):
    with open(path, 'rb') as fh:
        data = fh.read()
    try:
        text = data.decode('utf-8', errors='replace')
    except Exception:
        text = str(data)
    return text.splitlines()

def best_pick(group: pd.DataFrame) -> pd.Series:
    g_ok = group[group['status']=='ok']
    if len(g_ok) > 0:
        g_ok_time = g_ok[g_ok['total_sec'].notna()]
        if len(g_ok_time) > 0:
            return g_ok_time.sort_values('total_sec').iloc[0]
        return g_ok.iloc[0]
    with_info = group[(group['total_sec'].notna()) | (group['max_rss_mb'].notna())]
    if len(with_info)>0:
        return with_info.iloc[0]
    return group.iloc[0]

def make_time_plot(subdf: pd.DataFrame, dataset_name: str, r_val: int, out_dir: str):
    d = subdf[(subdf['dataset_name']==dataset_name) & (subdf['r']==r_val)].copy()
    d = d.dropna(subset=['s','total_sec']).sort_values('s')
    if d.empty: return None
    xs = d['s'].astype(int).tolist()
    ys = d['total_sec'].astype(float).tolist()
    fig = plt.figure(figsize=(8,4.8))
    plt.bar([str(x) for x in xs], ys)
    plt.xlabel("s")
    plt.ylabel("Total time (s)")
    plt.title(f"{dataset_name}   r={r_val}")
    plt.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.5)
    plt.tight_layout()
    fname = f"{dataset_name}_r{r_val}_time_overall.png"
    path = os.path.join(out_dir, fname)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path

def make_mem_plot(subdf: pd.DataFrame, dataset_name: str, r_val: int, out_dir: str):
    d = subdf[(subdf['dataset_name']==dataset_name) & (subdf['r']==r_val)].copy()
    d = d.dropna(subset=['s','max_rss_mb']).sort_values('s')
    if d.empty: return None
    xs = d['s'].astype(int).tolist()
    ys = d['max_rss_mb'].astype(float).tolist()
    fig = plt.figure(figsize=(8,4.8))
    plt.bar([str(x) for x in xs], ys)
    plt.xlabel("s")
    plt.ylabel("Max RSS (MB)")
    plt.title(f"{dataset_name}   r={r_val}")
    plt.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.5)
    plt.tight_layout()
    fname = f"{dataset_name}_r{r_val}_mem.png"
    path = os.path.join(out_dir, fname)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path

def arcData(base_root_dir=BASE_ROOT_DIR) -> pd.DataFrame:
    # Gather all candidate folders: all bazel_runs_* under base_root_dir
    folders = _candidate_dirs(None, base_root_dir)

    # Find candidate log files
    log_files = []
    for folder in folders:
        for root, _, files in os.walk(folder):
            for f in files:
                if f.lower().endswith((".log", ".out", ".txt")):
                    log_files.append(os.path.join(root, f))

    # Parse each log
    records = []
    for p in sorted(log_files):
        lines = read_lines(p)
        dataset, r_val, s_val = scan_text_for_dataset_r_s(lines)
        if (dataset is None) or (r_val is None) or (s_val is None):
            f_ds, f_r, f_s = scan_filename_for_dataset_r_s(p)
            dataset = dataset or f_ds
            r_val = r_val if r_val is not None else f_r
            s_val = s_val if s_val is not None else f_s
        elapsed, maxrss_mb, exit_status = scan_text_for_time_and_mem(lines)
        threads = scan_text_for_threads(lines)
        status = 'ok' if (exit_status == 0 or exit_status is None) else 'fail'
        if (dataset is None) or (r_val is None) or (s_val is None):
            continue
        records.append({
            'log_path': p,
            'dataset_file': dataset,
            'dataset_name': dataset_nice_name(dataset),
            'r': r_val,
            's': s_val,
            'total_sec': elapsed,
            'max_rss_mb': maxrss_mb,
            'exit_status': exit_status,
            'status': status,
            'threads': threads,
        })

    df = pd.DataFrame.from_records(records)
    if df.empty:
        return df
    df = df.sort_values(['dataset_name','r','s'])

    # Deduplicate multiple attempts for the same (dataset,r,s)
    group_keys = ['dataset_name', 'r', 's', 'threads'] if 'threads' in df.columns else ['dataset_name', 'r', 's']
    dedup = df.groupby(group_keys, as_index=False).apply(best_pick).reset_index(drop=True)
    # keep only successful exits if exit_status is present
    dedup = dedup[(dedup['exit_status'].isna()) | (dedup['exit_status'].isin([0]))].reset_index(drop=True)
    return dedup

def main():
    # Scan BASE_DIR and all 'bazel_runs_*' folders under BASE_ROOT_DIR
    dedup = arcData()
    print(dedup.columns.tolist())
    if dedup.empty:
        print(f"No runs parsed under {BASE_ROOT_DIR}. Please check parsing patterns (dataset/r/s) or filename pattern r*_s*_<dataset>.")
        return

    rows = []
    for dataset_name in sorted(dedup['dataset_name'].unique()):
        for r_val in [1,2,3,4]:
            p_time = make_time_plot(dedup, dataset_name, r_val, OUT_DIR)
            p_mem  = make_mem_plot(dedup, dataset_name, r_val, OUT_DIR)
            rows.append({'dataset': dataset_name, 'r': r_val,
                         'plot_time': p_time or '', 'plot_mem': p_mem or ''})

    idx = pd.DataFrame(rows)
    idx.to_csv(os.path.join(OUT_DIR, "_plots_index.csv"), index=False)

    with open(os.path.join(OUT_DIR, "_preview.html"), "w", encoding="utf-8") as f:
        f.write("<html><body><h2>Bazel Runs: Time & Memory Plots</h2>\n")
        for _, row in idx.iterrows():
            if row['plot_time']:
                f.write(f"<div><p><b>{row['dataset']} (r={row['r']}) time</b></p>"
                        f"<img src='{os.path.basename(row['plot_time'])}' width='650'></div>\n")
            if row['plot_mem']:
                f.write(f"<div><p><b>{row['dataset']} (r={row['r']}) memory</b></p>"
                        f"<img src='{os.path.basename(row['plot_mem'])}' width='650'></div>\n")
        f.write("</body></html>")

    print("Done. Index:", os.path.join(OUT_DIR, "_plots_index.csv"))
    print("Preview:", os.path.join(OUT_DIR, "_preview.html"))

if __name__ == "__main__":
    main()