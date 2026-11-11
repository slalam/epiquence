#!/usr/bin/env python3
"""Parse SeqKit stats output into MultiQC-compatible TSV."""
import argparse, csv

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    with open(args.input, 'r') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))

    total_reads = sum(int(r.get('num_seqs', 0)) for r in rows)
    total_bases = sum(int(r.get('sum_len', 0)) for r in rows)
    min_len = min((int(r.get('min_len', 0)) for r in rows), default=0)
    max_len = max((int(r.get('max_len', 0)) for r in rows), default=0)
    avg_len = total_bases / total_reads if total_reads > 0 else 0
    gc_sum = sum(float(r.get('GC(%)', r.get('GC', 0))) * int(r.get('sum_len', 0)) for r in rows)
    avg_gc = gc_sum / total_bases if total_bases > 0 else 0

    with open(args.output, 'w') as f:
        f.write("Sample\tTotal_Reads\tTotal_Bases\tMin_Len\tMax_Len\tAvg_Len\tGC_pct\n")
        f.write(f"{args.sample}\t{total_reads}\t{total_bases}\t{min_len}\t{max_len}\t{avg_len:.1f}\t{avg_gc:.2f}\n")

if __name__ == '__main__':
    main()