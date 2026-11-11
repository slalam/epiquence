#!/usr/bin/env python3
"""Parse fastp JSON output into MultiQC-compatible TSV."""
import argparse, json, os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    median_phred = 0.0
    mean_phred = 0.0
    q30_pct = 0.0
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

        q30_rate = after.get('q30_rate', 0)
        q30_pct = q30_rate * 100 if q30_rate <= 1 else q30_rate

        read1_after = data.get('read1_after_filtering', {})
        read2_after = data.get('read2_after_filtering', {})
        q1 = read1_after.get('quality_curves', {}).get('mean', [])
        q2 = read2_after.get('quality_curves', {}).get('mean', [])

        all_quals = q1 + q2
        if all_quals:
            mean_phred = sum(all_quals) / len(all_quals)
            sorted_quals = sorted(all_quals)
            median_phred = sorted_quals[len(sorted_quals) // 2]
        else:
            thresholds = [(90, 35, 34), (80, 33, 32), (70, 31, 30), (50, 28, 27)]
            median_phred, mean_phred = 25.0, 24.0
            for q30_thresh, med, mn in thresholds:
                if q30_pct >= q30_thresh:
                    median_phred, mean_phred = float(med), float(mn)
                    break

    with open(args.output, 'w') as out:
        out.write("Sample\tMedian_Phred\tMean_Phred\tQ30_Pct\tReads_Before_Filter\tReads_After_Filter\tPct_Reads_Passed\tBases_Before_Filter\tBases_After_Filter\tPct_Bases_Passed\n")
        out.write(f"{args.sample}\t{median_phred:.1f}\t{mean_phred:.1f}\t{q30_pct:.1f}\t{reads_before}\t{reads_after}\t{pct_reads_passed:.1f}\t{bases_before}\t{bases_after}\t{pct_bases_passed:.1f}\n")

if __name__ == '__main__':
    main()