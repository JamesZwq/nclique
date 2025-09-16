#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_nucleus_logs.py
解析 run_nucleus.sh 生成的 nucleus_experiment.log，
输出与 make_plots_bazel.py 风格相容的 DataFrame：
- dataset_name：仅文件名且去扩展名
- r, s：由 k 拆出（34->r=3,s=4；114->r=1,s=14）
"""
from __future__ import annotations
import re
from pathlib import Path
from typing import Dict, Any, List
import pandas as pd
import numpy as np
from dateutil import parser as date_parser

KNOWN_EXTS = {".edges", ".adj", ".txt", ".mtx", ".csv", ".tsv", ".bin", ".gz"}

def _wall_time_to_sec(s: str) -> float:
    s = s.strip()
    if not s:
        return np.nan
    parts = s.split(":")
    try:
        if len(parts) == 3:
            h, m, sec = parts
            return int(h) * 3600 + int(m) * 60 + float(sec)
        elif len(parts) == 2:
            m, sec = parts
            return int(m) * 60 + float(sec)
        else:
            return float(parts[0])
    except Exception:
        return np.nan

def _dataset_nice_name(path: str | None) -> str | None:
    if not path:
        return None
    name = Path(path.strip().strip('"').strip("'")).name
    while True:
        stem, suffix = Path(name).stem, Path(name).suffix
        if suffix and suffix.lower() in KNOWN_EXTS and stem:
            name = stem
            continue
        if suffix.lower() == ".gz":
            name = Path(name).with_suffix("").name
            continue
        break
    return name

def _extract_time_stats(block: str) -> Dict[str, Any]:
    metrics: Dict[str, Any] = {}
    patterns = [
        (r"^Elapsed \(wall clock\) time.*?:\s*(.+)$", "wall_time_sec", _wall_time_to_sec),
        (r"^User time \(seconds\):\s*(.+)$",          "user_time_sec", float),
        # (r"^System time \(seconds\):\s*(.+)$",        "sys_time_sec", float),
        # (r"^Percent of CPU this job got:\s*([\d\.]+)\s*%$", "cpu_percent", float),
        (r"^Maximum resident set size \(kbytes\):\s*(\d+)$", "max_rss_mb", lambda s: int(s)/1024.0),
        (r"^Exit status:\s*(\d+)$",                   "exit_status_inner", int),
    ]
    for line in block.splitlines():
        line = line.strip()
        for pat, col, conv in patterns:
            m = re.match(pat, line)
            if m:
                try:
                    metrics[col] = conv(m.group(1))
                except Exception:
                    metrics[col] = np.nan
    return metrics

def _split_rs_from_k(k_val) -> tuple[int | None, int | None]:
    """34->(3,4)；114->(1,14)。无法解析则返回(None,None)。"""
    if k_val is None:
        return None, None
    try:
        s = str(int(k_val))
    except Exception:
        s = str(k_val).strip()
        if not s.isdigit():
            return None, None
    if len(s) >= 2:
        return int(s[0]), int(s[1:])
    return None, None

def arcData(log_path: str | Path = "nucleus_experiment.log") -> pd.DataFrame:
    log_path = Path(log_path)
    text = log_path.read_text(encoding="utf-8", errors="ignore")

    blocks = text.split("================== RUN ==================")
    records: List[Dict[str, Any]] = []

    for blk in blocks[1:]:
        blk_main = blk.split("================ END OF RUN ==================", 1)[0] if "================ END OF RUN ==================" in blk else blk

        def _grab(pattern: str, cast=lambda x: x.strip(), default=None):
            m = re.search(pattern, blk_main)
            return cast(m.group(1)) if m else default

        ts_str = _grab(r"timestamp\s*:\s*(.+)")
        try:
            timestamp = date_parser.parse(ts_str) if ts_str else None
        except Exception:
            timestamp = ts_str

        dataset_path = _grab(r"dataset\s*:\s*(.+)")
        dataset_name = _dataset_nice_name(dataset_path)
        k = _grab(r"\nk\s*:\s*(\d+)", cast=int)
        option = _grab(r"option\s*:\s*(YES|NO)")
        command = _grab(r"command\s*:\s*(.+)")

        m_tb = re.search(r"-{4,}\s*/bin/time\s+-v\s*-{4,}\s*(.*?)\-{4,}\s+STATUS\s+-{4,}", blk_main, flags=re.S|re.M)
        if not m_tb:
            m_tb = re.search(r"/bin/time\s+-v.*?\n(.*?)\n-+\s+STATUS", blk_main, flags=re.S|re.M)
        time_block = m_tb.group(1).strip() if m_tb else None

        status_shell = _grab(r"exit_status:\s*(\d+)", cast=int)
        timed_out = bool(re.search(r"exit_status:\s*\d+\s*\(TIMEOUT\)", blk_main))

        r_val, s_val = _split_rs_from_k(k)

        rec = {
            # "timestamp": timestamp,
            "dataset_name": dataset_name,
            "r": r_val,
            "s": s_val,
            "option": option,
            # "command": command,
            "exit_status": status_shell,
            "timed_out": timed_out,
        }
        if time_block:
            rec.update(_extract_time_stats(time_block))
        records.append(rec)

    if not records:
        raise ValueError("No runs parsed — 请检查日志格式。")

    df = pd.DataFrame(records)
    # 列顺序 & 排序
    preferred = ["timestamp","dataset_name","r","s","option","command",
                 "wall_time_sec","user_time_sec","sys_time_sec","cpu_percent","max_rss_mb",
                 "exit_status","exit_status_inner","timed_out"]
    df = df[[c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]]
    df = df.sort_values([c for c in ["dataset_name","r","s","option","timestamp"] if c in df.columns], kind="stable")
    for c in ["wall_time_sec","user_time_sec","sys_time_sec","cpu_percent","max_rss_mb","exit_status","exit_status_inner","r","s"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

if __name__ == "__main__":
    # import sys
    # path = sys.argv[1] if len(sys.argv) > 1 else "nucleus_experiment.log"
    d = arcData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/nucleus_experiment.log")
    out = Path("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/nucleus_experiment.log").with_suffix(".with_rs.csv")
    d.to_csv(out, index=False)
    print(f"[Saved] {out}")