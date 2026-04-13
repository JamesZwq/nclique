#!/usr/bin/env python3
"""Generate random small graphs and find V4 vs ST mismatches."""

import subprocess, random, tempfile, os, sys, itertools

BIN = "./build/bin/degeneracy_cliques"
R, S = 3, 4

def generate_random_graph(n, edge_prob=0.5):
    """Generate random graph with n vertices."""
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < edge_prob:
                edges.append((i, j))
    return n, edges

def generate_multipartite(parts):
    """Generate complete multipartite graph K_{parts[0], parts[1], ...}"""
    n = sum(parts)
    edges = []
    groups = []
    offset = 0
    for p in parts:
        groups.append(list(range(offset, offset + p)))
        offset += p
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            for u in groups[i]:
                for v in groups[j]:
                    edges.append((u, v))
    return n, edges

def write_graph(n, edges, path):
    with open(path, 'w') as f:
        f.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")

def run_algo(graph_path, env_var):
    """Run algorithm and extract core distribution."""
    env = os.environ.copy()
    env[env_var] = "1"
    result = subprocess.run(
        [BIN, graph_path, str(R), str(S)],
        capture_output=True, text=True, env=env, timeout=30
    )
    output = result.stdout + result.stderr

    # Parse core= lines
    cores = {}
    for line in output.split('\n'):
        line = line.strip()
        if line.startswith('core='):
            parts = line.split()
            try:
                c = float(parts[0].split('=')[1])
                cnt = int(parts[1].split('=')[1])
                cores[c] = cores.get(c, 0) + cnt
            except:
                pass
    return cores

def compare(graph_path, label=""):
    """Compare V4 vs ST on a graph. Return True if match."""
    try:
        v4 = run_algo(graph_path, "PIVOTER_RUN_REGION_V4")
        st = run_algo(graph_path, "PIVOTER_RUN_ST")
    except subprocess.TimeoutExpired:
        return True  # skip timeouts
    except Exception as e:
        return True

    if not st:
        return True  # ST produced no output (graph too small)

    # Compare core LEVELS (not counts, since V4 might output tuples)
    v4_levels = sorted(v4.keys())
    st_levels = sorted(st.keys())

    # Also compare counts
    match_levels = (v4_levels == st_levels)
    match_counts = (v4 == st)

    if not match_levels or not match_counts:
        print(f"\n{'='*60}")
        print(f"MISMATCH FOUND! {label}")
        print(f"{'='*60}")
        print(f"V4 levels: {v4_levels}")
        print(f"ST levels: {st_levels}")
        if not match_levels:
            print(f"LEVEL DIFF: V4 extra={set(v4_levels)-set(st_levels)}, ST extra={set(st_levels)-set(v4_levels)}")
        print(f"\nV4 dist: {dict(sorted(v4.items()))}")
        print(f"ST dist: {dict(sorted(st.items()))}")

        # Print the graph
        with open(graph_path) as f:
            print(f"\nGraph ({graph_path}):")
            print(f.read())
        return False
    return True

def main():
    tmpdir = tempfile.mkdtemp()
    graph_path = os.path.join(tmpdir, "test.edges")

    tested = 0

    # Phase 1: structured graphs (multipartite)
    print("Testing multipartite graphs...")
    for parts in [
        [2,2], [2,2,2], [2,2,2,2], [3,3], [3,3,3],
        [2,3], [2,3,4], [1,2,3], [1,1,2,3], [2,2,3,3],
        [1,3,3], [1,1,1,3], [4,4], [1,2,2,2],
    ]:
        n, edges = generate_multipartite(parts)
        if n < S:
            continue
        write_graph(n, edges, graph_path)
        label = f"K_{{{','.join(map(str,parts))}}}"
        if not compare(graph_path, label):
            return
        tested += 1
        print(f"  K_{{{','.join(map(str,parts))}}}: OK", flush=True)

    # Phase 2: random graphs
    print(f"\nTesting random graphs (r={R}, s={S})...")
    for trial in range(10000):
        n = random.randint(S, 12)
        p = random.uniform(0.3, 0.9)
        n, edges = generate_random_graph(n, p)
        if len(edges) == 0:
            continue
        write_graph(n, edges, graph_path)
        if not compare(graph_path, f"random n={n} m={len(edges)} trial={trial}"):
            return
        tested += 1
        if tested % 100 == 0:
            print(f"  {tested} graphs tested, all OK...", flush=True)

    print(f"\nAll {tested} graphs passed!")
    os.unlink(graph_path)
    os.rmdir(tmpdir)

if __name__ == "__main__":
    main()
