#!/usr/bin/env python3
"""Maintain a cumulative per-species isolate metadata TSV for surveillance trees."""

import argparse
import csv
import os
import sys
from datetime import datetime

def now_in(tz_name):
    """Current time in the configured IANA zone, e.g. America/Chicago."""
    if tz_name:
        try:
            from zoneinfo import ZoneInfo
            return datetime.now(ZoneInfo(tz_name))
        except Exception as exc:
            print(f'WARN: timezone {tz_name!r} unavailable ({exc}); '
                  f'falling back to system local time', file=sys.stderr)
    return datetime.now()

FIELDS = [
    'sample',
    'isolate_id',
    'assembly_type',
    'species_name',
    'mlst_scheme',
    'mlst_st',
    'chromosomal_amr_genes',
    'chromosomal_amr_gene_count',
    'workflow_id',
    'group_id',
    'group_description',
    'last_seen_date',
    'first_seen_workflow',
    'first_seen_date',
]


def read_tsv(path):
    """Read a TSV into a list of dicts."""
    if not path or path in ('NO_PREV_ISOLATES', 'NO_AMR') or not os.path.exists(path):
        return []
    with open(path, newline='', encoding='utf-8') as f:
        return [r for r in csv.DictReader(f, delimiter='\t')]


def load_chromosomal_amr(path):
    """chromosomal_amr.tsv -> {sample: 'geneA;geneB'}."""
    amr = {}
    for row in read_tsv(path):
        sample = (row.get('Sample') or '').strip()
        genes = (row.get('Chromosomal_AMR_Genes') or '').strip()
        if sample and genes and genes != '-':
            amr[sample] = genes
    return amr


def gene_count(genes):
    if not genes:
        return '0'
    return str(len([g for g in genes.split(';') if g.strip()]))


def main():
    p = argparse.ArgumentParser(
        description='Merge this run into the cumulative per-species isolate TSV')
    p.add_argument('--current', required=True,
                   help='This run\'s isolates TSV from PREPARE_ISOLATES_TSV')
    p.add_argument('--prev', default=None,
                   help='Previous cumulative TSV (omit or NO_PREV_ISOLATES on first run)')
    p.add_argument('--chromosomal-amr', default=None,
                   help='chromosomal_amr.tsv from GENERATE_MULTIQC_CONFIG')
    p.add_argument('--species', required=True)
    p.add_argument('--workflow-id', default='')
    p.add_argument('--group-id', default='')
    p.add_argument('--group-description', default='')
    p.add_argument('--output', required=True)
    p.add_argument('--timezone', default=None,
                   help='IANA timezone for first_seen/last_seen stamps')
    args = p.parse_args()

    now = now_in(args.timezone).strftime('%Y-%m-%d %H:%M:%S %Z').strip()
    wf = args.workflow_id or 'unknown'

    # Existing cumulative rows, keyed by sample.
    merged = {}
    for row in read_tsv(args.prev):
        sample = (row.get('sample') or '').strip()
        if sample:
            merged[sample] = {k: row.get(k, '') for k in FIELDS}

    n_before = len(merged)

    amr_map = load_chromosomal_amr(args.chromosomal_amr)
    have_amr = bool(args.chromosomal_amr and args.chromosomal_amr != 'NO_AMR'
                    and os.path.exists(args.chromosomal_amr))

    n_new = 0
    n_updated = 0
    for row in read_tsv(args.current):
        sample = (row.get('sample') or '').strip()
        if not sample:
            continue

        prev_row = merged.get(sample, {})
        if have_amr:
            # AMR file present: absence from it means the isolate genuinely carries no chromosomal
            # resistance genes, so refresh unconditionally.
            genes = amr_map.get(sample, '')
        else:
            # No AMR file this run.
            genes = prev_row.get('chromosomal_amr_genes', '')
        rec = {
            'sample': sample,
            'isolate_id': (row.get('isolate_id') or '').strip(),
            'assembly_type': (row.get('assembly_type') or '').strip(),
            'species_name': args.species,
            'mlst_scheme': (row.get('mlst_scheme') or '').strip(),
            'mlst_st': (row.get('mlst_st') or '').strip(),
            'chromosomal_amr_genes': genes,
            'chromosomal_amr_gene_count': gene_count(genes),
            'workflow_id': wf,
            'group_id': args.group_id,
            'group_description': args.group_description,
            'last_seen_date': now,
        }

        prev = prev_row
        if prev and prev.get('first_seen_workflow'):
            # Sample seen before: preserve original detection, refresh the rest.
            rec['first_seen_workflow'] = prev['first_seen_workflow']
            rec['first_seen_date'] = prev.get('first_seen_date', '')
            n_updated += 1
        else:
            rec['first_seen_workflow'] = wf
            rec['first_seen_date'] = now
            n_new += 1

        merged[sample] = rec

    with open(args.output, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter='\t',
                           extrasaction='ignore', restval='')
        w.writeheader()
        for sample in sorted(merged):
            w.writerow(merged[sample])

    print(f'[{args.species}] cumulative isolates: {n_before} -> {len(merged)} '
          f'({n_new} new, {n_updated} refreshed)', file=sys.stderr)
    if amr_map:
        n_with_amr = sum(1 for r in merged.values() if r.get('chromosomal_amr_genes'))
        print(f'[{args.species}] chromosomal AMR genes on {n_with_amr}/{len(merged)} isolates',
              file=sys.stderr)


if __name__ == '__main__':
    main()