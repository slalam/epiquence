#!/usr/bin/env python3
"""Outbreak clustering: group samples by species, compute Mash distances, perform hierarchical
clustering, and prepare per-cluster directories."""
import argparse, csv, os, shutil, json
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib.clustering import UnionFind
from epiquence_lib.mash import read_matrix


def parse_species_map(gtdbtk_dir):
    """Parse GTDB-Tk parsed TSVs to build sample -> species map."""
    import re
    species_map = {}
    if not os.path.isdir(gtdbtk_dir):
        return species_map
    for f in os.listdir(gtdbtk_dir):
        if f.endswith('_gtdbtk_parsed.tsv'):
            with open(os.path.join(gtdbtk_dir, f)) as fh:
                for row in csv.DictReader(fh, delimiter='\t'):
                    sample = row.get('Sample', '').strip()
                    species = row.get('species', '').strip()
                    genus = row.get('genus', '').strip()
                    # Remove GTDB polyphyletic suffixes (e.g. _A, _B, _C, _D)
                    if species:
                        species = re.sub(r'_[A-Z]$', '', species)
                        species = re.sub(r'_[A-Z] ', ' ', species)
                    if genus:
                        genus = re.sub(r'_[A-Z]$', '', genus)
                    if sample and species:
                        species_map[sample] = species
                    elif sample and genus:
                        species_map[sample] = f"{genus} sp."
    return species_map


def parse_metadata(metadata_file):
    """Parse metadata CSV/TSV into dict keyed by sample_id."""
    metadata = {}
    if not metadata_file or not os.path.exists(metadata_file):
        return metadata
    sep = '\t' if metadata_file.lower().endswith('.tsv') else ','
    with open(metadata_file) as f:
        reader = csv.DictReader(f, delimiter=sep)
        for row in reader:
            sid = row.get('sample_id', row.get('Sample', row.get('sample', '')))
            if sid:
                metadata[sid.strip()] = {k: v for k, v in row.items()}
    return metadata


def parse_mash_matrix(matrix_file):
    """Parse Mash distance matrix TSV."""
    row_names, col_names, rows = read_matrix(matrix_file, basename=True, strip_fasta=True)
    df_data = dict(zip(row_names, rows))
    return col_names, df_data


def cluster_samples(sample_names, dist_matrix, sample_order, threshold=0.015):
    """Single-linkage clustering using union-find (no scipy required)."""
    n = len(sample_names)
    # Build index map from sample_order
    order_idx = {s: i for i, s in enumerate(sample_order)}

    uf = UnionFind(list(range(n)))

    # Merge samples within threshold
    for i in range(n):
        for j in range(i + 1, n):
            si, sj = sample_names[i], sample_names[j]
            oi, oj = order_idx.get(si), order_idx.get(sj)
            if oi is not None and oj is not None and si in dist_matrix:
                d = dist_matrix[si][oj]
                if d < threshold:
                    uf.union(i, j)

    groups = {}
    for i, smp in enumerate(sample_names):
        root = uf.find(i)
        groups.setdefault(root, []).append(smp)

    # Renumber clusters starting from 1
    result = {}
    for idx, (_, members) in enumerate(sorted(groups.items()), 1):
        result[idx] = members
    return result


def main():
    parser = argparse.ArgumentParser(description='Outbreak clustering by species')
    parser.add_argument('--gtdbtk-dir', required=True, help='Directory with GTDB-Tk parsed TSVs')
    parser.add_argument('--assemblies-dir', required=True, help='Directory with chromosome.fasta files')
    parser.add_argument('--mash-matrix', required=True, help='Mash distance matrix TSV')
    parser.add_argument('--metadata', default=None, help='Metadata CSV/TSV file')
    parser.add_argument('--threshold', type=float, default=0.015, help='Mash distance clustering threshold')
    parser.add_argument('--outdir', default='.', help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    species_map = parse_species_map(args.gtdbtk_dir)
    metadata = parse_metadata(args.metadata) if args.metadata else {}

    print(f"Species map: {species_map}")

    sample_order, dist_data = parse_mash_matrix(args.mash_matrix)
    print(f"Mash samples: {sample_order}")

    # Group samples by species (only species with >= 2 samples)
    species_groups = {}
    for sample, species in species_map.items():
        if species:
            species_groups.setdefault(species, []).append(sample)
    species_groups = {sp: samps for sp, samps in species_groups.items() if len(samps) >= 2}

    print(f"Species groups (2+ samples): {species_groups}")

    if not species_groups:
        print("No species with 2+ samples found. Skipping outbreak analysis.")
        with open(os.path.join(args.outdir, 'outbreak_summary.json'), 'w') as f:
            json.dump({'species_groups': {}, 'clusters': {}}, f, indent=2)
        with open(os.path.join(args.outdir, 'cluster_list.csv'), 'w') as f:
            f.write('cluster_name,cluster_dir,num_samples,species,run_parsnp,run_pairwise\n')
        return

    sample_order, dist_data = parse_mash_matrix(args.mash_matrix)

    # For each species group, cluster and create directories
    all_clusters = {}
    cluster_counter = 0

    for species, samples in sorted(species_groups.items()):
        # Filter to samples present in Mash matrix
        valid_samples = [s for s in samples if s in dist_data]
        if len(valid_samples) < 2:
            continue

        species_slug = species.replace(' ', '_').replace('.', '')
        groups = cluster_samples(valid_samples, dist_data, sample_order, args.threshold)

        species_cluster_num = 0
        for cluster_id, members in sorted(groups.items()):
            if len(members) < 2:
                continue  # Skip singleton clusters
            species_cluster_num += 1
            cluster_counter += 1
            cluster_name = f"{species_slug}_cluster_{species_cluster_num}"
            cluster_dir = os.path.join(args.outdir, 'clusters', cluster_name)
            os.makedirs(cluster_dir, exist_ok=True)

            for sample in members:
                src = os.path.join(args.assemblies_dir, f"{sample}.fasta")
                if os.path.exists(src):
                    shutil.copy(src, cluster_dir)

            all_clusters[cluster_name] = {
                'species': species,
                'samples': members,
                'size': len(members),
                'run_ska': True
            }

    summary = {
        'species_groups': {sp: samps for sp, samps in species_groups.items()},
        'clusters': all_clusters,
        'metadata': metadata,
        'threshold': args.threshold
    }
    with open(os.path.join(args.outdir, 'outbreak_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    # Write cluster list for Nextflow to iterate
    with open(os.path.join(args.outdir, 'cluster_list.csv'), 'w') as f:
        f.write('cluster_name,cluster_dir,num_samples,species,run_ska,run_cga\n')
        for name, info in sorted(all_clusters.items()):
            cdir = os.path.join(args.outdir, 'clusters', name)
            run_cga = info['size'] >= 5
            f.write(f"{name},{cdir},{info['size']},{info['species']},{info['run_ska']},{run_cga}\n")

    print(f"\n=== Outbreak Clustering Summary ===")
    print(f"Species with 2+ samples: {len(species_groups)}")
    print(f"Total clusters: {len(all_clusters)}")
    for name, info in sorted(all_clusters.items()):
        print(f"  {name}: {info['size']} samples ({info['species']})")
        for s in info['samples']:
            print(f"    - {s}")


if __name__ == '__main__':
    main()