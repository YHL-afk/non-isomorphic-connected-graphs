#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from math import comb, ceil, sqrt
import csv
import sys
import networkx as nx
from networkx.generators.atlas import graph_atlas_g
import matplotlib.pyplot as plt


def connected_labeled_count(n: int) -> int:
    """
    c_n = 2^{C(n,2)} - sum_{k=1..n-1} C(n-1,k-1) * c_k * 2^{C(n-k,2)},  c_1=1
    """
    if n < 1:
        return 0
    c = [0] * (n + 1)
    c[1] = 1
    for m in range(2, n + 1):
        total = 2 ** comb(m, 2)
        s = 0
        for k in range(1, m):
            s += comb(m - 1, k - 1) * c[k] * (2 ** comb(m - k, 2))
        c[m] = total - s
    return c[n]


def get_unlabeled_connected(n: int):

    if n > 7:
        return []
    atlas = graph_atlas_g()  
    Gs = [G for G in atlas if G.number_of_nodes() == n and nx.is_connected(G)]
    return Gs

def bucket_by_degree_sequence(graphs):
    buckets = {}
    for G in graphs:
        degseq = tuple(sorted([d for _, d in G.degree()], reverse=True))
        buckets.setdefault(degseq, []).append(G)
    return buckets

def save_bucket_csv(buckets, path):
    rows = []
    for degseq, group in buckets.items():
        rows.append({
            "degseq": " ".join(map(str, degseq)),
            "count": len(group)
        })
    rows.sort(key=lambda r: (tuple(-int(x) for x in r["degseq"].split()), -r["count"]))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["degseq", "count"])
        w.writeheader()
        w.writerows(rows)


def save_graph6(graphs, path):
    with open(path, "w", encoding="utf-8") as f:
        for G in graphs:
            g6 = nx.to_graph6_bytes(G, header=False).decode().strip()
            f.write(g6 + "\n")


def draw_montage(graphs, out_path, cols=None, layout="spring",
                 node_size=10, edge_width=0.6, dpi=300):
    n = len(graphs)
    if n == 0:
        print("[Warn] No graphs to draw.", file=sys.stderr)
        return
    if cols is None:
        cols = int(ceil(sqrt(n)))  
    rows = int(ceil(n / cols))

    fig = plt.figure(figsize=(cols * 1.0, rows * 1.0), dpi=dpi)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0.05, hspace=0.05)

    for idx, G in enumerate(graphs, start=1):
        ax = fig.add_subplot(rows, cols, idx)
        ax.axis("off")
        if layout == "spring":
            pos = nx.spring_layout(G, seed=0, iterations=50)
        elif layout == "kamada":
            pos = nx.kamada_kawai_layout(G)
        elif layout == "circular":
            pos = nx.circular_layout(G)
        else:
            pos = nx.spectral_layout(G)
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_size, linewidths=0)
        nx.draw_networkx_edges(G, pos, ax=ax, width=edge_width)

    fig.savefig(out_path, bbox_inches="tight", pad_inches=0)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Count/draw connected simple graphs on n labeled/unlabeled vertices."
    )
    parser.add_argument("-n", type=int, default=7, help="number of vertices (default: 7)")
    parser.add_argument("--layout", type=str, default="spring",
                        choices=["spring", "kamada", "circular", "spectral"],
                        help="layout for montage (n<=7)")
    parser.add_argument("--cols", type=int, default=None, help="columns in montage grid")
    parser.add_argument("--g6", type=str, default=None,
                        help="export graph6 of non-isomorphic connected graphs (n<=7)")
    parser.add_argument("--bucket_csv", type=str, default=None,
                        help="export CSV (degree-sequence buckets) for n<=7")
    parser.add_argument("--montage", type=str, default=None,
                        help="output PNG for montage (n<=7)")
    args = parser.parse_args()

    n = args.n

    labeled = connected_labeled_count(n)
    print(f"[Connected] labeled count c_{n} = {labeled:,}")


    graphs = get_unlabeled_connected(n)
    if n <= 7:
        print(f"[Connected] unlabeled (non-isomorphic) count on n={n}: {len(graphs)}")

        if args.bucket_csv:
            buckets = bucket_by_degree_sequence(graphs)
            save_bucket_csv(buckets, args.bucket_csv)
            print(f"[Output] degree-sequence buckets -> {args.bucket_csv} "
                  f"(unique degseqs: {len(buckets)})")

        if args.g6:
            save_graph6(graphs, args.g6)
            print(f"[Output] graph6 list -> {args.g6}")

        if args.montage:
            draw_montage(graphs, args.montage, cols=args.cols, layout=args.layout)
            print(f"[Output] montage image -> {args.montage}")
    else:
        print("[Note] Unlabeled/atlas-based outputs are available only for n ≤ 7.")

if __name__ == "__main__":
    main()

