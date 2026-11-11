#!/usr/bin/env python3
"""Parse AMRFinderPlus output into a single-row MultiQC-compatible TSV."""
import argparse
import csv
import os

SUBTYPE_AMR = 'AMR'
SUBTYPE_POINT = 'POINT'
SCOPE_PLUS = 'plus'


def read_hits(path):
    """Yield AMRFinderPlus rows as dicts."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return
    with open(path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            yield row


def main():
    parser = argparse.ArgumentParser(description='Parse AMRFinderPlus results')
    parser.add_argument('--sample', required=True, help='Sample name (e.g. SAMPLE-01_S)')
    parser.add_argument('--input', required=True, help='AMRFinderPlus TSV')
    parser.add_argument('--output', required=True, help='Output single-row TSV')
    args = parser.parse_args()

    core_genes = set()
    plus_genes = set()
    mutations = set()
    gene_classes = set()

    for row in read_hits(args.input):
        subtype = (row.get('Subtype') or '').strip().upper()
        symbol = (row.get('Element symbol') or '').strip()
        if not symbol:
            continue
        if subtype == SUBTYPE_AMR:
            # Scope separates the curated core set from the --plus additions (efflux pumps, stress
            # response).
            scope = (row.get('Scope') or '').strip().lower()
            if scope == SCOPE_PLUS:
                plus_genes.add(symbol)
            else:
                core_genes.add(symbol)
            cls = (row.get('Class') or '').strip()
            if cls:
                gene_classes.add(cls)
        elif subtype == SUBTYPE_POINT:
            mutations.add(symbol)

    def joined(values):
        return '; '.join(sorted(values)) if values else '-'

    all_genes = core_genes | plus_genes

    with open(args.output, 'w') as out:
        out.write('Sample\tNum_AMR_Genes\tAMR_Genes\tNum_Core_AMR_Genes\tCore_AMR_Genes\t'
                  'Num_Plus_AMR_Genes\tPlus_AMR_Genes\tNum_Point_Mutations\t'
                  'Point_Mutations\tAMR_Classes\n')
        out.write('\t'.join([
            args.sample,
            str(len(all_genes)), joined(all_genes),
            str(len(core_genes)), joined(core_genes),
            str(len(plus_genes)), joined(plus_genes),
            str(len(mutations)), joined(mutations),
            joined(gene_classes),
        ]) + '\n')


if __name__ == '__main__':
    main()