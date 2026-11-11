#!/usr/bin/env python3
"""Rename assembly contigs to CONTIG_N format."""
import argparse
import re


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    # First pass: read all contigs to get sequences
    contigs = []  # list of (original_header, sequence_lines)
    current_header = None
    current_lines = []
    with open(args.input, 'r') as fin:
        for line in fin:
            if line.startswith('>'):
                if current_header is not None:
                    contigs.append((current_header, current_lines))
                current_header = line[1:].strip()
                current_lines = []
            else:
                current_lines.append(line)
        if current_header is not None:
            contigs.append((current_header, current_lines))

    # Second pass: write with new names, lengths, and circularity
    with open(args.output, 'w') as fout:
        for i, (header, seq_lines) in enumerate(contigs, 1):
            seq_len = sum(len(line.strip()) for line in seq_lines)

            # Parse existing metadata from header
            parts = header.split(None, 1)
            extra_info = parts[1] if len(parts) > 1 else ''

            # Detect circularity from header metadata
            is_circular = bool(re.search(
                r'rotated[= _][Tt]rue|circular[= _]true|topology[= _]circular'
                r'|\[circularised\]|circular_contig|\bcircular\b',
                extra_info, re.IGNORECASE
            ))

            new_header = f">CONTIG_{i} length={seq_len}"
            if is_circular:
                new_header += " circular=true"

            # Preserve other metadata, stripping rotated/circular/length flags we already handled
            cleaned_extra = re.sub(
                r'\s*rotated=[Tt]rue\s*|\s*circular=[Tt]rue\s*'
                r'|\s*topology=circular\s*|\s*\[circularised\]\s*'
                r'|\s*length=\d+\s*',
                ' ', extra_info
            ).strip()
            if cleaned_extra:
                new_header += f" {cleaned_extra}"

            fout.write(new_header + '\n')
            for seq_line in seq_lines:
                fout.write(seq_line)

    print(f"Renamed {len(contigs)} contigs")


if __name__ == '__main__':
    main()