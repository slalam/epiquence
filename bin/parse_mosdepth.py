#!/usr/bin/env python3
"""Parse mosdepth output into MultiQC-compatible TSV."""
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--global-dist', required=True)
    parser.add_argument('--summary', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    coverage_dist = {}
    with open(args.global_dist, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] == 'total':
                coverage_dist[int(fields[1])] = float(fields[2])

    mean_cov = 0.0
    with open(args.summary, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] == 'total':
                try: mean_cov = float(fields[3])
                except (ValueError, IndexError): pass

    median_cov = 0
    for cov in sorted(coverage_dist.keys()):
        if coverage_dist[cov] >= 0.5:
            median_cov = cov
        else:
            break

    thresholds = [20, 30, 40, 50, 100]
    threshold_pcts = {t: coverage_dist.get(t, 0) * 100 for t in thresholds}

    with open(args.output, 'w') as out:
        out.write("Sample\tMean_Coverage\tMedian_Coverage\tPct_ge_20x\tPct_ge_30x\tPct_ge_40x\tPct_ge_50x\tPct_ge_100x\n")
        out.write(f"{args.sample}\t{round(mean_cov)}\t{median_cov}\t{threshold_pcts[20]:.2f}\t{threshold_pcts[30]:.2f}\t{threshold_pcts[40]:.2f}\t{threshold_pcts[50]:.2f}\t{threshold_pcts[100]:.2f}\n")

if __name__ == '__main__':
    main()