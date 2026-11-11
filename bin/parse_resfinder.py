#!/usr/bin/env python3
"""Parse a ResFinder ``pheno_table.txt`` into a single-row MultiQC TSV."""
import argparse
import os

COLUMNS = ('Antimicrobial', 'Class', 'WGS-predicted phenotype', 'Match',
           'Genetic background')
EXCLUDED_CLASSES = {'under_development'}


def read_pheno_table(path):
    """Yield dict rows from a ResFinder pheno_table.txt."""
    if not path or not os.path.exists(path):
        return
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split('\t')
            row = {c: (parts[i].strip() if i < len(parts) else '')
                   for i, c in enumerate(COLUMNS)}
            yield row


def is_reportable(row):
    """False for unnameable drugs and not-yet-validated classes."""
    drug = row.get('Antimicrobial', '').strip().lower()
    cls = row.get('Class', '').strip().lower()
    if not drug or drug.startswith('unknown'):
        return False
    if cls in EXCLUDED_CLASSES:
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description='Parse ResFinder phenotype table')
    parser.add_argument('--sample', required=True, help='Sample name (e.g. SAMPLE-01_S)')
    parser.add_argument('--input', required=True, help='ResFinder pheno_table.txt')
    parser.add_argument('--output', required=True, help='Output single-row TSV')
    args = parser.parse_args()

    resistant = []
    classes = set()
    n_tested = 0

    for row in read_pheno_table(args.input):
        if not is_reportable(row):
            continue
        n_tested += 1
        if row.get('WGS-predicted phenotype', '').strip().lower() == 'resistant':
            resistant.append(row['Antimicrobial'])
            if row.get('Class'):
                classes.add(row['Class'])

    with open(args.output, 'w') as out:
        out.write('Sample\tNum_Antimicrobials_Tested\tNum_Resistant\t'
                  'Resistant_Antimicrobials\tResistant_Classes\n')
        out.write('\t'.join([
            args.sample,
            str(n_tested),
            str(len(resistant)),
            '; '.join(sorted(resistant)) if resistant else '-',
            '; '.join(sorted(classes)) if classes else '-',
        ]) + '\n')


if __name__ == '__main__':
    main()