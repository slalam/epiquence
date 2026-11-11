#!/usr/bin/env python3
"""Turn `sylph query` output for one assembly into EPIquence's per-sample taxonomy file."""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from epiquence_lib.gtdb import RANKS, accession_of, load_lineages, parse_lineage

HEADER = (['Sample'] + RANKS + ['classification', 'closest_genome_reference',
          'closest_genome_ani', 'closest_genome_af', 'method', 'call_status',
          'second_species', 'second_ani', 'containment'])


def read_hits(path):
    """sylph query TSV -> list of {accession, ani, containment}, best first."""
    hits = []
    with open(path, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            raw = (row.get('Naive_ANI') or '').strip()
            if raw in ('', 'NA'):
                raw = (row.get('Adjusted_ANI') or '').strip()
            try:
                ani = float(raw)
            except ValueError:
                continue
            genome = row.get('Genome_file') or ''
            hits.append({
                'accession': accession_of(genome) or os.path.basename(genome),
                'ani': ani,
                'containment': (row.get('Containment_ind') or '').strip(),
            })
    hits.sort(key=lambda h: h['ani'], reverse=True)
    return hits


def classify(hits, lineages, species_ani, borderline_ani, min_margin):
    """Return the output row (minus Sample) for a list of hits."""
    empty = {r: '' for r in RANKS}
    if not hits:
        return empty, 'Unclassified', '', '', 'unclassified', '', '', ''

    top = hits[0]
    taxonomy, deepest = parse_lineage(lineages.get(top['accession'], ''))
    if not deepest:
        print(f"WARN: no GTDB lineage for {top['accession']}; check --taxonomy "
              f"matches the sylph database release", file=sys.stderr)

    if top['ani'] >= species_ani:
        status = 'species'
    elif top['ani'] >= borderline_ani:
        status = 'borderline'
    else:
        status = 'genus'
        taxonomy['species'] = ''
        deepest = taxonomy['genus'] or deepest

    # Best hit of a different species, for ambiguity and for the report.
    second_species, second_ani = '', ''
    top_species = parse_lineage(lineages.get(top['accession'], ''))[0]['species']
    for h in hits[1:]:
        sp = parse_lineage(lineages.get(h['accession'], ''))[0]['species']
        if sp and sp != top_species:
            second_species, second_ani = sp, f"{h['ani']:.2f}"
            if (status in ('species', 'borderline') and h['ani'] >= borderline_ani
                    and top['ani'] - h['ani'] < min_margin):
                status = 'ambiguous'
            break

    return (taxonomy, deepest, top['accession'], f"{top['ani']:.2f}", status,
            second_species, second_ani, top['containment'])


def main():
    p = argparse.ArgumentParser(description='Parse sylph query output for one assembly')
    p.add_argument('--sample', required=True, help='Sample name, e.g. SAMPLE-01_HL')
    p.add_argument('--query', required=True, help='sylph query TSV')
    p.add_argument('--taxonomy', nargs='+', required=True,
                   help='GTDB taxonomy TSV(s), e.g. bac120_taxonomy_r226.tsv.gz ar53_taxonomy_r226.tsv.gz')
    p.add_argument('--species-ani', type=float, default=95.0)
    p.add_argument('--borderline-ani', type=float, default=94.0)
    p.add_argument('--min-margin', type=float, default=0.5)
    p.add_argument('--output', required=True)
    args = p.parse_args()

    hits = read_hits(args.query)
    # Only the accessions that can influence the call need a lineage lookup.
    lineages = load_lineages(args.taxonomy, [h['accession'] for h in hits[:20]])
    taxonomy, deepest, ref, ani, status, sp2, ani2, cont = classify(
        hits, lineages, args.species_ani, args.borderline_ani, args.min_margin)

    with open(args.output, 'w', encoding='utf-8') as out:
        out.write('\t'.join(HEADER) + '\n')
        out.write('\t'.join([args.sample] + [taxonomy[r] for r in RANKS] +
                            [deepest, ref, ani, '-', 'sylph', status, sp2, ani2, cont]) + '\n')

    print(f'{args.sample}: {deepest or "unclassified"} (ANI {ani or "-"}, {status}, '
          f'{len(hits)} hit(s))', file=sys.stderr)


if __name__ == '__main__':
    main()
