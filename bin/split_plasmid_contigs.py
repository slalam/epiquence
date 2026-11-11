#!/usr/bin/env python3
"""Split a MOB-recon plasmid bin into one FASTA per contig."""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from epiquence_lib.names import contig_token  # noqa: E402


def read_fasta(path):
    """Yield (header_without_gt, [sequence lines])."""
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, seq
                header, seq = line[1:], []
            elif header is not None:
                seq.append(line)
    if header is not None:
        yield header, seq


def main():
    ap = argparse.ArgumentParser(description='Split a MOB plasmid bin into per-contig FASTAs')
    ap.add_argument('--input', required=True, help='MOB-recon plasmid_*.fasta')
    ap.add_argument('--prefix', required=True,
                    help='Output stem, e.g. SAMPLE-03_HL_plasmid_AA102')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--min-length', type=int, default=0,
                    help='Skip contigs shorter than this (0 = keep everything)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    written, skipped = 0, []
    seen = set()

    for header, seq in read_fasta(args.input):
        contig_id = header.split()[0] if header.split() else ''
        token = contig_token(contig_id)

        # Two contigs collapsing to one token would silently overwrite each other, so disambiguate
        # rather than lose a sequence.
        base_token = token
        n = 2
        while token in seen:
            token = f"{base_token}x{n}"
            n += 1
        seen.add(token)

        length = sum(len(s.strip()) for s in seq)
        if args.min_length and length < args.min_length:
            skipped.append(f"{contig_id} ({length} bp)")
            continue

        out = os.path.join(args.outdir, f"{args.prefix}_{token}.fasta")
        with open(out, 'w') as f:
            f.write('>' + header + '\n')
            for line in seq:
                f.write(line + '\n')
        written += 1

    print(f"  {os.path.basename(args.input)} -> {written} contig FASTA(s)", file=sys.stderr)
    if skipped:
        print(f"    below --min-length, skipped: {', '.join(skipped)}", file=sys.stderr)
    if written == 0:
        print(f"    WARN: nothing written from {args.input}", file=sys.stderr)


if __name__ == '__main__':
    main()