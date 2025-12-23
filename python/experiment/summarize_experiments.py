#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re, sys, pathlib, pandas as pd

START = re.compile(r"^=+\s*RUN\s*=+$", re.M)
END   = re.compile(r"^=+\s*END OF RUN\s*=+$", re.M)

FIELDS = {
    "timestamp": r"^timestamp\s*:\s*(.+)$",
    "dataset":   r"^dataset\s*:\s*(.+)$",
    "s":         r"^s\s*:\s*(\d+)\s*$",
    "r":         r"^r\s*:\s*(\d+)\s*$",
    "cmd":       r"^cmd\s*:\s*(.+)$",
}
RE = {k: re.compile(v, re.M|re.I) for k,v in FIELDS.items()}

METRICS = {
    "tree_build_ms":    re.compile(r"Tree Build took:\s*([\d.]+)\s*ms", re.I),
    "prog_time_ms":     re.compile(r"\btime\s*:\s*([\d.]+)\s*ms\b", re.I),
    "num_leaves":       re.compile(r"(?:nu[mn])\s*Leaf\s*:\s*([0-9]+)", re.I),
    "clique_index_ms":  re.compile(r"clique Index build took:\s*([\d.]+)\s*ms", re.I),
    "max_rss_kb":       re.compile(r"Maximum resident set size \(kbytes\):\s*([0-9]+)", re.I),
    "exit_status":      re.compile(r"Exit status:\s*([0-9]+)", re.I),
}

def norm_cmd(s: str) -> str:
    # 压缩多空白，去首尾空白
    return re.sub(r"\s+", " ", s.strip())

def parse_block(txt: str) -> dict:
    row = {k: None for k in FIELDS}
    for k,pat in RE.items():
        m = pat.search(txt)
        if m:
            row[k] = m.group(1).strip()
    if row["cmd"]:
        row["cmd"] = norm_cmd(row["cmd"])

    for k,pat in METRICS.items():
        m = pat.search(txt)
        if not m:
            row[k] = None
        else:
            val = m.group(1)
            if k in ("num_leaves","max_rss_kb","exit_status"):
                row[k] = int(val)
            else:
                row[k] = float(val)
    return row

def info_score(row: pd.Series) -> int:
    keys = ["tree_build_ms","prog_time_ms","num_leaves","clique_index_ms","max_rss_kb","exit_status"]
    return int(sum(pd.notna(row.get(k)) for k in keys))

def main(in_path: str, out_dir: str):
    text = pathlib.Path(in_path).read_text(errors="replace")

    # 只取成对的 RUN…END 块
    starts = [m.end() for m in START.finditer(text)]
    ends   = [m.start() for m in END.finditer(text)]
    blocks = []
    j = 0
    for s in starts:
        while j < len(ends) and ends[j] < s:
            j += 1
        if j >= len(ends):
            break
        e = ends[j]
        blocks.append(text[s:e])
        j += 1

    rows = []
    for b in blocks:
        r = parse_block(b)
        rows.append(r)

    if not rows:
        print("No completed RUN blocks found.", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(rows)

    # 类型整理
    for col in ("s","r"):
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    # 去重键：dataset, s, r, cmd
    for key in ("dataset","cmd"):
        df[key] = df[key].fillna("").map(lambda x: norm_cmd(x))

    # 选择每组最佳
    def pick(g: pd.DataFrame) -> pd.Series:
        g = g.copy()
        # 分数：信息多者优先
        g["__score"] = g.apply(info_score, axis=1)
        # 其次有无 exit_status（无也行），最后用时间戳（字典序）选最新
        g = g.sort_values(["__score","timestamp"], ascending=[False, True])
        return g.iloc[0]

    dedup = df.groupby(["dataset","s","r","cmd"], dropna=False, as_index=False).apply(pick).reset_index(drop=True)
    dedup = dedup.drop(columns=["__score"], errors="ignore")

    outp = pathlib.Path(out_dir); outp.mkdir(parents=True, exist_ok=True)
    cols = ["dataset","s","r","tree_build_ms","prog_time_ms","num_leaves","clique_index_ms","max_rss_kb","exit_status","cmd","timestamp"]
    for c in cols:
        if c not in dedup.columns: dedup[c] = pd.NA
    dedup = dedup[cols].sort_values(["dataset","s","r","timestamp"]).reset_index(drop=True)
    dedup.to_csv(outp/"summary_all_runs_sorted.csv", index=False)
    # r=3 的子表
    dedup[dedup["r"]==3][["dataset","s","r","tree_build_ms","clique_index_ms","prog_time_ms","num_leaves","max_rss_kb","exit_status","cmd","timestamp"]] \
        .to_csv(outp/"summary_r_eq_3.csv", index=False)

    # 简单重复检查（应为 0）
    dup_counts = dedup.groupby(["dataset","s","r","cmd"]).size().reset_index(name="n")
    still_dup = dup_counts[dup_counts["n"]>1]
    if not still_dup.empty:
        print("WARNING: still duplicated keys:\n", still_dup.to_string(index=False))

    print(outp/"summary_all_runs_sorted.csv")
    print(outp/"summary_r_eq_3.csv")

if __name__ == "__main__":
    main("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdata", "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/outdir")