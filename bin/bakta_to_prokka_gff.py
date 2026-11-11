#!/usr/bin/env python3
"""Convert Bakta GFF3 + FNA files to Prokka-compatible combined GFF3 format for Panaroo."""
import argparse
import os
import sys


def convert_bakta_gff(input_gff, fna_path, output_gff):
    """Convert a single Bakta GFF3 + FNA pair to Prokka-compatible format."""
    gff_lines = []

    # Read GFF3: keep only ## directives and feature lines
    with open(input_gff) as f:
        for line in f:
            if line.startswith('##FASTA'):
                break
            if line.startswith('##'):
                gff_lines.append(line)
            elif not line.startswith('#') and line.strip():
                gff_lines.append(line)

    # Read FASTA from the .fna file
    fasta_lines = []
    if fna_path and os.path.exists(fna_path):
        with open(fna_path) as f:
            for line in f:
                fasta_lines.append(line)

    if not fasta_lines:
        print(f"  WARNING: No FASTA found for {os.path.basename(input_gff)}", file=sys.stderr)

    # Write combined GFF + FASTA (Prokka style)
    with open(output_gff, 'w') as out:
        for line in gff_lines:
            out.write(line)
        out.write('##FASTA\n')
        for line in fasta_lines:
            out.write(line)
        if fasta_lines and not fasta_lines[-1].endswith('\n'):
            out.write('\n')


def main():
    parser = argparse.ArgumentParser(description='Convert Bakta GFF3 + FNA to Prokka-compatible format')
    parser.add_argument('input_gffs', nargs='+', help='Input Bakta GFF3 files')
    parser.add_argument('--fna-dir', default='.', help='Directory containing matching .fna files')
    parser.add_argument('--outdir', default='converted_gffs', help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for gff in args.input_gffs:
        basename = os.path.basename(gff)
        prefix = basename.replace('.gff3', '').replace('.gff', '')
        fna_candidates = [
            os.path.join(args.fna_dir, f'{prefix}.fna'),
            os.path.join(args.fna_dir, f'{prefix}.fasta'),
            os.path.join(args.fna_dir, f'{prefix}.fa'),
        ]
        gff_dir = os.path.dirname(gff)
        if gff_dir:
            fna_candidates.extend([
                os.path.join(gff_dir, f'{prefix}.fna'),
                os.path.join(gff_dir, f'{prefix}.fasta'),
            ])

        fna_path = None
        for candidate in fna_candidates:
            if os.path.exists(candidate):
                fna_path = candidate
                break

        outpath = os.path.join(args.outdir, basename)
        convert_bakta_gff(gff, fna_path, outpath)

        n_features = 0
        n_seqs = 0
        with open(outpath) as f:
            in_fasta = False
            for line in f:
                if line.startswith('##FASTA'):
                    in_fasta = True
                elif in_fasta and line.startswith('>'):
                    n_seqs += 1
                elif not in_fasta and not line.startswith('#') and '\t' in line:
                    n_features += 1
        print(f"  {basename}: {n_features} features, {n_seqs} sequences, fna={'found' if fna_path else 'MISSING'}")

    print(f"Converted {len(args.input_gffs)} GFF3 files to {args.outdir}/")


if __name__ == '__main__':
    main()