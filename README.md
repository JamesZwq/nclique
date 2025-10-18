

# CBND — Clique‑Path Based $(r,s)$‑Nucleus Decomposition

CBND is a counting‑based implementation of $(r,s)$‑nucleus decomposition that builds and updates a **Clique Path Index (CPI)** instead of enumerating all $s$‑cliques. The approach keeps the index compact and supports path‑local updates used by the peeling routine and (optionally) the hierarchy construction.

> TL;DR: Build with CMake, then run
> ```bash
> ./build/bin/degeneracy_cliques <graph.edges> <s> <r>
> ```
> on an undirected edge list. No batch scripts are required.

---

## 1) Build

Requirements:
- Linux/macOS, CMake ≥ 3.16, a modern C++ compiler (GCC/Clang)

```bash
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The main binary will be created at:
```
build/bin/degeneracy_cliques
```

---

## 2) Input graph format

- **Undirected simple graph** as an edge list text file.
- One edge per line: `u v` (space or tab).
- Vertex IDs are integers; 0‑ or 1‑based are both fine as long as the file is consistent.
- Self‑loops and multi‑edges should be removed beforehand.

Example (`toy.edges`):
```
1 2
2 3
2 4
3 4
4 5
```

---

## 3) Quick start

Run $(r,s)=(3,4)$ on the toy graph:
```bash
./build/bin/degeneracy_cliques toy.edges 4 3
```

General form:
```bash
./build/bin/degeneracy_cliques <graph.edges> <s> <r>
```
- `<graph.edges>`: path to the edge‑list file
- `<s>`: target clique size (e.g., 4)
- `<r>`: base clique size for support (e.g., 3)

The program prints a short summary to `stdout` (and error messages to `stderr`).

> Tip: if you need wall‑clock and peak memory for a single run, wrap the command with `/usr/bin/time -v` (optional).

---

## 4) What gets computed

- $(r,s)$‑nucleus decomposition using CPI‑based counting.
- One‑shot support initialization from CPI, then path‑local updates as peeling proceeds.
- (If enabled in your build) a hierarchy can be maintained via a stack of disjoint‑set unions (DSUs). No special command‑line flags are required for basic runs.

---

## 5) Examples

```bash
# (2,3) — truss‑like setting
./build/bin/degeneracy_cliques graph.edges 3 2

# (3,4) — higher‑order setting
./build/bin/degeneracy_cliques graph.edges 4 3

# (1,3) — vertex‑centric higher‑order core
./build/bin/degeneracy_cliques graph.edges 3 1
```

---

## 6) Performance notes

- Build in `Release` mode for full optimizations.
- Larger graphs and higher parameters benefit from fast I/O and sufficient RAM.
- CPI avoids explicit $s$‑clique enumeration; in practice runtime is dominated by path‑local updates and support maintenance.

---

## 7) Troubleshooting

- **`std::bad_alloc` / out‑of‑memory**: ensure the graph is cleaned (no duplicates/self‑loops) and run in `Release` build.
- **Unexpected IDs**: verify the edge list contains only integer vertex IDs and consistent indexing.
- **Slow I/O**: place the graph on a local SSD and avoid network filesystems when benchmarking.

---

## 8) License

See `LICENSE` in this repository.