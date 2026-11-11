#!/usr/bin/env python3
"""Cluster plasmids by Mash distance using union-find."""
import argparse, csv, sys
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib.clustering import UnionFind
from epiquence_lib.mash import read_matrix, pairwise_dict


def parse_mash_matrix(path):
    """Parse Mash distance matrix TSV."""
    row_names, col_names, rows = read_matrix(path, basename=False, strip_fasta=False)
    dist = pairwise_dict(row_names, col_names, rows)
    return row_names, dist


def cluster_union_find(samples, dist, threshold):
    """Union-find clustering by distance threshold."""
    uf = UnionFind(samples)

    for i, a in enumerate(samples):
        for b in samples[i+1:]:
            d = dist.get((a, b), 1.0)
            if d <= threshold:
                uf.union(a, b)

    return uf.groups()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mash-matrix', required=True)
    parser.add_argument('--threshold', type=float, default=0.005)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    samples, dist = parse_mash_matrix(args.mash_matrix)
    if not samples:
        with open(args.output, 'w') as f:
            f.write('cluster_id\tplasmid\tcluster_size\n')
        return

    clusters = cluster_union_find(samples, dist, args.threshold)

    # Filter to clusters with 2+ members and assign IDs
    multi = {k: v for k, v in clusters.items() if len(v) >= 2}
    sorted_clusters = sorted(multi.items(), key=lambda x: -len(x[1]))

    with open(args.output, 'w') as f:
        f.write('cluster_id\tplasmid\tcluster_size\n')
        for idx, (root, members) in enumerate(sorted_clusters, 1):
            cid = f"hgt_cluster_{idx}"
            for m in sorted(members):
                f.write(f"{cid}\t{m}\t{len(members)}\n")

    print(f"Found {len(sorted_clusters)} plasmid clusters (2+ members) from {len(samples)} plasmids")
    for idx, (root, members) in enumerate(sorted_clusters, 1):
        print(f"  hgt_cluster_{idx}: {len(members)} plasmids")
        for m in sorted(members):
            print(f"    - {m}")


if __name__ == '__main__':
    main()