#!/usr/bin/env python3
"""Parse VirulenceFinder results into a MultiQC-compatible TSV."""
import argparse
import csv
import os
import sys


def main():
    parser = argparse.ArgumentParser(description='Parse VirulenceFinder results')
    parser.add_argument('--sample', required=True, help='Sample name (e.g., SAMPLE-01_L)')
    parser.add_argument('--virfinder-dir', required=True, help='VirulenceFinder output directory')
    parser.add_argument('--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    # VirulenceFinder outputs results_tab.tsv with columns: Database, Virulence factor, Identity,
    # Query / Template length, Contig, Position in contig, Protein function, Accession number
    results_file = os.path.join(args.virfinder_dir, 'results_tab.tsv')

    virulence_factors = []
    if os.path.exists(results_file):
        with open(results_file) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                vf = row.get('Virulence factor', '').strip()
                if vf:
                    virulence_factors.append(vf)

    unique_factors = sorted(set(virulence_factors))
    vf_count = len(unique_factors)
    vf_list = '; '.join(unique_factors) if unique_factors else '-'

    with open(args.output, 'w') as f:
        f.write("Sample\tVirulence_Factor_Count\tVirulence_Factor_List\n")
        f.write(f"{args.sample}\t{vf_count}\t{vf_list}\n")


if __name__ == '__main__':
    main()