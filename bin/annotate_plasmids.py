#!/usr/bin/env python3
"""Annotate assembly contigs with plasmid cluster IDs from MOB-recon."""
import argparse, csv, os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--mob-dir', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    contig_report_file = f"{args.mob_dir}/contig_report.txt"

    # Build contig-to-cluster mapping from MOB-recon contig_report
    contig_to_cluster = {}
    if os.path.exists(contig_report_file) and os.path.getsize(contig_report_file) > 0:
        with open(contig_report_file, 'r') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                contig_id_full = row.get('contig_id', '').strip()
                contig_id = contig_id_full.split()[0] if contig_id_full else ''
                cluster_id = row.get('primary_cluster_id', '-').strip()
                molecule_type = row.get('molecule_type', '').strip()
                if contig_id and molecule_type == 'plasmid' and cluster_id and cluster_id != '-':
                    contig_to_cluster[contig_id] = cluster_id

    # Read assembly and append plasmid= tag only to plasmid contigs
    annotated = 0
    with open(args.input, 'r') as fin, open(args.output, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                header = line[1:].strip()
                contig_name = header.split(None, 1)[0]
                cluster = contig_to_cluster.get(contig_name, '')
                if cluster:
                    fout.write(f">{header} plasmid={cluster}\n")
                    annotated += 1
                else:
                    fout.write(line)
            else:
                fout.write(line)

    print(f"Annotated {annotated} plasmid contigs out of {len(contig_to_cluster)} mapped")

if __name__ == '__main__':
    main()