#!/usr/bin/env python3
"""General-purpose compression baseline (bytes only, no timing): the per-vertex S trees of every graph whose index file is in
cx/ are materialized by strees_dump, once with the degeneracy labels of the input (the labels the decomposition works with)
and once with the index's aligned labels, and compressed with zstd -19, zstd -19 --long=31 and xz -9e; the chain index file
itself is compressed the same way.  Writes compress.json; refuses to overwrite.  Needs build/strees_dump, zstd and xz."""
import datetime, json, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent; OUT = HERE / 'compress.json'
TOOLS = {'zstd19': ['zstd', '-19', '-T1', '-q', '-c'], 'zstd19long': ['zstd', '-19', '--long=31', '-T1', '-q', '-c'], 'xz9e': ['xz', '-9e', '-T1', '-c']}

def compressed_sizes(path):
    return {name: len(subprocess.run(cmd + [str(path)], capture_output=True, check=True).stdout) for name, cmd in TOOLS.items()}

def main():
    if OUT.exists(): raise SystemExit('refusing to overwrite compress.json')
    recorded = {}
    for f in ('final.json', 'more.json'):
        for r in json.loads((HERE / f).read_text())['runs']:
            if 'result' in r: recorded[Path(r['graph']).stem] = r['result']
    rec = {'started': datetime.datetime.now().astimezone().isoformat(), 'tools': {k: ' '.join(v) for k, v in TOOLS.items()}, 'graphs': {}}
    tmp = Path('/tmp/strees_dump.bin')
    for cx in sorted((HERE / 'cx').glob('*.cx')):
        g = cx.stem
        if g not in recorded: print('skip', g, '(no recorded run)'); continue
        print(g, flush=True)
        entry = {'index_file_bytes': cx.stat().st_size, 'index_compressed': compressed_sizes(cx), 'index_bytes_total': recorded[g]['bytes_total'],
                 'baseline_vertex_bytes': recorded[g]['baseline_vertex_bytes'], 'count_bits': recorded[g]['count_bits']}
        for labels, extra in (('degeneracy', [str(cx) + '.perm']), ('aligned', [])):
            dump = json.loads(subprocess.run([str(HERE / 'build' / 'strees_dump'), str(cx), str(tmp)] + extra, capture_output=True, text=True, check=True).stdout)
            assert dump['labels'] == labels
            entry['dump'] = dump; entry[f'strees_compressed_{labels}'] = compressed_sizes(tmp)
        rec['graphs'][g] = entry; OUT.write_text(json.dumps(rec, indent=1) + '\n')
        d, a = entry['strees_compressed_degeneracy'], entry['strees_compressed_aligned']
        print(f"  S trees {dump['strees_bytes']:,} -> xz {d['xz9e']:,} (degeneracy labels) / {a['xz9e']:,} (aligned labels); index {entry['index_file_bytes']:,} -> xz {entry['index_compressed']['xz9e']:,}", flush=True)
    tmp.unlink(missing_ok=True); rec['finished'] = datetime.datetime.now().astimezone().isoformat(); OUT.write_text(json.dumps(rec, indent=1) + '\n')

if __name__ == '__main__':
    main()
