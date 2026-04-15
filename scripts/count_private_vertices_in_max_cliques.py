#!/usr/bin/env python3
"""Count vertices that appear in exactly one large maximal clique.

Default output is a compact three-line summary:
  vertices
  edges
  private_vertices_count
  private_works_as_cloud
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple


Vertex = int
Clique = Tuple[Vertex, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Enumerate maximal cliques, keep those with size >= min-size, "
            "and count how many vertices appear in exactly one such clique."
        )
    )
    parser.add_argument("graph", help="Path to graph file")
    parser.add_argument(
        "--min-size",
        type=int,
        default=4,
        help="Only count maximal cliques with size >= this value (default: 4)",
    )
    parser.add_argument(
        "--show-vertices",
        action="store_true",
        help="Print the vertices that appear in exactly one kept maximal clique",
    )
    parser.add_argument(
        "--per-clique",
        action="store_true",
        help="Print private-vertex counts for each kept maximal clique",
    )
    parser.add_argument(
        "--show-cliques",
        action="store_true",
        help="Print all kept maximal cliques",
    )
    return parser.parse_args()


def read_graph(path: Path) -> Dict[Vertex, Set[Vertex]]:
    lines: List[Tuple[int, int]] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"Invalid line: {raw.rstrip()}")
            lines.append((int(parts[0]), int(parts[1])))

    if not lines:
        return {}

    edges = lines

    # Support both plain edge list and "n m" header + edge list.
    first_u, first_v = lines[0]
    maybe_edges = lines[1:]
    if (
        first_u >= 0
        and first_v >= 0
        and len(maybe_edges) == first_v
        and maybe_edges
        and max(max(u, v) for u, v in maybe_edges) <= first_u
    ):
        edges = maybe_edges

    adj: Dict[Vertex, Set[Vertex]] = defaultdict(set)
    for u, v in edges:
        if u == v:
            continue
        if v in adj[u]:
            continue
        adj[u].add(v)
        adj[v].add(u)

    # Keep isolated vertices if the file had a header.
    if edges is maybe_edges and first_u > 0:
        for v in range(first_u):
            adj.setdefault(v, set())

    return dict(adj)


def bron_kerbosch(
    r: Set[Vertex],
    p: Set[Vertex],
    x: Set[Vertex],
    adj: Dict[Vertex, Set[Vertex]],
    min_size: int,
    out: List[Clique],
) -> None:
    if len(r) + len(p) < min_size:
        return

    if not p and not x:
        if len(r) >= min_size:
            out.append(tuple(sorted(r)))
        return

    pivot_candidates = p | x
    pivot = max(pivot_candidates, key=lambda v: len(p & adj[v])) if pivot_candidates else None
    blocked = adj[pivot] if pivot is not None else set()

    for v in list(p - blocked):
        bron_kerbosch(r | {v}, p & adj[v], x & adj[v], adj, min_size, out)
        p.remove(v)
        x.add(v)


def enumerate_maximal_cliques(adj: Dict[Vertex, Set[Vertex]], min_size: int) -> List[Clique]:
    cliques: List[Clique] = []
    bron_kerbosch(set(), set(adj), set(), adj, min_size, cliques)
    cliques.sort(key=lambda c: (-len(c), c))
    return cliques


def clique_edges(clique: Sequence[Vertex]) -> Set[Tuple[Vertex, Vertex]]:
    edges: Set[Tuple[Vertex, Vertex]] = set()
    for i, u in enumerate(clique):
        for v in clique[i + 1 :]:
            a, b = (u, v) if u < v else (v, u)
            edges.add((a, b))
    return edges


def main() -> None:
    args = parse_args()
    graph_path = Path(args.graph)
    adj = read_graph(graph_path)
    num_vertices = len(adj)
    num_edges = sum(len(nbrs) for nbrs in adj.values()) // 2

    cliques = enumerate_maximal_cliques(adj, args.min_size)
    membership = Counter()
    covered_edge_set: Set[Tuple[Vertex, Vertex]] = set()
    for clique in cliques:
        membership.update(clique)
        covered_edge_set.update(clique_edges(clique))

    covered_vertices = sorted(membership)
    covered_vertex_count = len(covered_vertices)
    covered_edge_count = len(covered_edge_set)
    private_vertices = sorted(v for v, cnt in membership.items() if cnt == 1)
    private_ratio = (
        len(private_vertices) / covered_vertex_count if covered_vertex_count else 0.0
    )
    hist = Counter(membership.values())

    if not (args.show_vertices or args.per_clique or args.show_cliques):
        print(f"vertices: {covered_vertex_count}")
        print(f"edges: {covered_edge_count}")
        print(f"private_vertices_count: {len(private_vertices)}")
        print(f"private_works_as_cloud: {private_ratio:.6%}")
        return

    print(f"graph: {graph_path}")
    print(f"original_vertices: {num_vertices}")
    print(f"original_edges: {num_edges}")
    print(f"vertices: {covered_vertex_count}")
    print(f"edges: {covered_edge_count}")
    print(f"min_size: {args.min_size}")
    print(f"maximal_cliques_kept: {len(cliques)}")
    print(f"covered_vertices: {covered_vertex_count}")
    print(f"private_vertices_count: {len(private_vertices)}")
    print(f"private_works_as_cloud: {private_ratio:.6%}")

    print("membership_histogram:")
    for freq in sorted(hist):
        print(f"  appears_in_{freq}_cliques: {hist[freq]}")

    if args.show_vertices:
        print("private_vertices:")
        print("  " + " ".join(map(str, private_vertices)) if private_vertices else "  (none)")

    if args.show_cliques:
        print("maximal_cliques:")
        for idx, clique in enumerate(cliques):
            print(f"  [{idx}] size={len(clique)} clique={' '.join(map(str, clique))}")

    if args.per_clique:
        print("per_clique_private_vertices:")
        for idx, clique in enumerate(cliques):
            private = [v for v in clique if membership[v] == 1]
            private_str = " ".join(map(str, private)) if private else "(none)"
            print(
                f"  [{idx}] size={len(clique)} private_count={len(private)} private={private_str}"
            )


if __name__ == "__main__":
    main()
