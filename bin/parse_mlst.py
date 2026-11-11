#!/usr/bin/env python3
"""Parse MLST output into MultiQC-compatible TSV."""
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    with open(args.input, 'r') as f:
        parts = f.readline().strip().split('\t')

    scheme = parts[1] if len(parts) > 1 else '-'
    st = parts[2] if len(parts) > 2 else '-'
    alleles = [parts[i] if i < len(parts) and parts[i] else '-' for i in range(3, 10)]
    alleles.extend(['-'] * (7 - len(alleles)))

    with open(args.output, 'w') as out:
        out.write("Sample\tMLST_Scheme\tMLST_ST\tLocus_1\tLocus_2\tLocus_3\tLocus_4\tLocus_5\tLocus_6\tLocus_7\n")
        out.write(f"{args.sample}\t{scheme}\t{st}\t{chr(9).join(alleles)}\n")

if __name__ == '__main__':
    main()