#!/usr/bin/env python3
"""Parse fastplong JSON output into MultiQC-compatible TSV."""
import argparse, json, os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    reads_before = 0
    reads_after = 0
    pct_reads_passed = 0.0
    bases_before = 0
    bases_after = 0
    pct_bases_passed = 0.0

    if os.path.exists(args.input):
        with open(args.input, 'r') as f:
            data = json.load(f)

        summary = data.get('summary', {})
        before = summary.get('before_filtering', {})
        after = summary.get('after_filtering', {})

        reads_before = int(before.get('total_reads', 0))
        reads_after = int(after.get('total_reads', 0))
        pct_reads_passed = (reads_after / reads_before * 100) if reads_before > 0 else 0.0

        bases_before = int(before.get('total_bases', 0))
        bases_after = int(after.get('total_bases', 0))
        pct_bases_passed = (bases_after / bases_before * 100) if bases_before > 0 else 0.0

    with open(args.output, 'w') as out:
        out.write("Sample\tReads_Before_Filter\tReads_After_Filter\tPct_Reads_Passed\tBases_Before_Filter\tBases_After_Filter\tPct_Bases_Passed\n")
        out.write(f"{args.sample}\t{reads_before}\t{reads_after}\t{pct_reads_passed:.1f}\t{bases_before}\t{bases_after}\t{pct_bases_passed:.1f}\n")

if __name__ == '__main__':
    main()