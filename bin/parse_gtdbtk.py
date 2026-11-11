#!/usr/bin/env python3
"""Parse GTDB-Tk output into MultiQC-compatible TSV."""
import argparse, csv, glob

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--gtdbtk-dir', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    summary_files = glob.glob(f"{args.gtdbtk_dir}/*.summary.tsv")
    rank_prefixes = {'d__': 'domain', 'p__': 'phylum', 'c__': 'class',
                     'o__': 'order', 'f__': 'family', 'g__': 'genus', 's__': 'species'}

    with open(args.output, 'w') as out:
        out.write("Sample\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tclassification\tclosest_genome_reference\tclosest_genome_ani\tclosest_genome_af\n")

        if summary_files:
            for sf in summary_files:
                with open(sf, 'r') as f:
                    for row in csv.DictReader(f, delimiter='\t'):
                        taxonomy = {v: '' for v in rank_prefixes.values()}
                        final_classification = ''

                        for part in row.get('classification', '').split(';'):
                            part = part.strip()
                            for prefix, rank in rank_prefixes.items():
                                if part.startswith(prefix):
                                    value = part[3:]
                                    # Remove GTDB polyphyletic suffixes (e.g. _A, _B, _C, _D)
                                    if value and rank in ('genus', 'species'):
                                        import re
                                        # Strip suffix from genus: "Campylobacter_D" -> "Campylobacter"
                                        value = re.sub(r'_[A-Z]$', '', value)
                                        # Strip suffix mid-string: "Campylobacter_D jejuni" -> "Campylobacter jejuni"
                                        value = re.sub(r'_[A-Z] ', ' ', value)
                                    taxonomy[rank] = value
                                    if value:
                                        final_classification = value
                                    break

                        cols = [args.sample] + [taxonomy[r] for r in ['domain','phylum','class','order','family','genus','species']]
                        cols += [final_classification, row.get('closest_genome_reference',''),
                                 row.get('closest_genome_ani',''), row.get('closest_genome_af','')]
                        out.write('\t'.join(cols) + '\n')
        else:
            out.write(f"{args.sample}" + "\t" * 11 + "\n")

if __name__ == '__main__':
    main()