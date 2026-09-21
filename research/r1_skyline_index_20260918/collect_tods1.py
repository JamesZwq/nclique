#!/usr/bin/env python3
"""After the tods1 queue has pushed its records (git pull first): copy the tods1 CND sweep records
(src-r1index/scripts/prior_original_<g>.json whose binary path is on tods1) into prior/tods1/, keeping the three
records reconstructed from the log, then regenerate the paper's tables and figures."""
import glob, json, os, shutil, subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
dst = HERE / 'prior' / 'tods1'; dst.mkdir(exist_ok=True)
for f in sorted(glob.glob(str(ROOT / 'src-r1index' / 'scripts' / 'prior_original_*.json'))):
    r = json.load(open(f))
    if r.get('binary', '').startswith('/data/wenqianz') and 'total_wall_s' in r and '_p' not in os.path.basename(f):
        shutil.copy(f, dst / os.path.basename(f)); print('copied', os.path.basename(f), r['s_max'], r['sizes_ok'], round(r['total_wall_s']))
paper = ROOT / 'Sigmod2027ChainIndex'
subprocess.run(['python3', str(paper / 'make_tables.py')], check=True)
subprocess.run(['python3', str(paper / 'make_figures.py')], check=True)
print('tables and figures regenerated; recompile the paper and check tables/prior_stats.tex and tables/query_stats.tex')
