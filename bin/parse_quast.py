#!/usr/bin/env python3
"""Parse QUAST report into MultiQC-compatible TSV."""
import argparse, os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--quast-dir', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    report_file = f"{args.quast_dir}/report.tsv"
    metrics = {'L50': 0, 'Largest contig': 0, '# contigs (>= 0 bp)': 0,
               '# contigs': 0, 'Total length': 0, 'N50': 0, 'GC (%)': 0.0}

    if os.path.exists(report_file):
        with open(report_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2 and parts[0] in metrics:
                    try:
                        metrics[parts[0]] = float(parts[1]) if '.' in parts[1] else int(parts[1])
                    except ValueError:
                        pass

    num_contigs = metrics['# contigs (>= 0 bp)'] or metrics['# contigs']
    largest_mbp = metrics['Largest contig'] / 1e6
    total_mbp = metrics['Total length'] / 1e6
    n50_kbp = metrics['N50'] / 1e3

    with open(args.output, 'w') as out:
        out.write("Sample\tL50\tLargest_Contig_Mbp\tNum_Contigs\tTotal_Length_Mbp\tN50_kbp\tGC_Percent\n")
        out.write(f"{args.sample}\t{metrics['L50']}\t{largest_mbp:.3f}\t{num_contigs}\t{total_mbp:.3f}\t{n50_kbp:.1f}\t{metrics['GC (%)']:.2f}\n")

if __name__ == '__main__':
    main()