#!/usr/bin/env python3
"""Parse NanoPlot stats into MultiQC-compatible TSV."""
import argparse, os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    median_phred = 0.0
    mean_phred = 0.0

    if os.path.exists(args.input):
        with open(args.input, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if '\t' in line:
                    parts = line.split('\t')
                    metric = parts[0].strip().lower()
                    value = parts[-1].strip()
                elif ':' in line:
                    parts = line.split(':', 1)
                    metric = parts[0].strip().lower()
                    value = parts[1].strip() if len(parts) > 1 else ''
                else:
                    continue

                try:
                    val_clean = value.replace(',', '').strip()
                    if metric in ('median_qual', 'median read quality', 'median quality'):
                        median_phred = float(val_clean)
                    elif metric in ('mean_qual', 'mean read quality', 'mean quality', 'average read quality'):
                        mean_phred = float(val_clean)
                except (ValueError, IndexError):
                    pass

    with open(args.output, 'w') as out:
        out.write("Sample\tMedian_Phred\tMean_Phred\n")
        out.write(f"{args.sample}\t{median_phred:.1f}\t{mean_phred:.1f}\n")

if __name__ == '__main__':
    main()