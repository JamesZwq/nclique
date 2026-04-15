
import subprocess, os, random, tempfile

BIN = './build/bin/degeneracy_cliques'
R, S = 3, 4
tmpf = tempfile.NamedTemporaryFile(mode='w', suffix='.edges', delete=False)

def run(path, env_var):
    env = os.environ.copy()
    env[env_var] = '1'
    out = subprocess.run([BIN, path, str(R), str(S)], capture_output=True, text=True, env=env, timeout=60)
    txt = out.stdout + out.stderr
    cores = {}
    found_first = False
    for line in txt.split('\n'):
        s = line.strip()
        if s.startswith('core=') and 'count=' in s:
            p = s.split()
            c = float(p[0].split('=')[1])
            cnt = int(float(p[1].split('=')[1]))
            if c in cores:
                break  # duplicate core level = second block, stop
            cores[c] = cnt
            found_first = True
        elif found_first and not s.startswith('core='):
            break
    return cores

for trial in range(5000):
    n = random.randint(8, 14)
    p = random.uniform(0.35, 0.85)
    edges = [(i,j) for i in range(n) for j in range(i+1,n) if random.random() < p]
    if not edges: continue
    with open(tmpf.name, 'w') as f:
        f.write(f'{n} {len(edges)}\n')
        for u,v in edges: f.write(f'{u} {v}\n')
    v2f = run(tmpf.name, 'PIVOTER_RUN_REGION_V2F')
    st = run(tmpf.name, 'PIVOTER_RUN_ST')
    if v2f != st:
        print(f'MISMATCH at trial={trial} n={n} m={len(edges)}')
        print(f'V2F: {dict(sorted(v2f.items()))}')
        print(f'ST:  {dict(sorted(st.items()))}')
        with open(tmpf.name) as f: print(f.read()[:500])
        break
    if trial % 200 == 0: print(f'{trial} OK (n={n})', flush=True)
else:
    print('All 5000 passed!')
os.unlink(tmpf.name)
