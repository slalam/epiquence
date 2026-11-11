#!/usr/bin/env python3
"""Select AMR-carrying plasmids from a per-assembly-type plasmid directory."""
import argparse
import os
import shutil
import sys


def parse_fasta_contigs(path):
    """Return (contig_names, header_amr_genes)."""
    names, genes = [], set()
    try:
        with open(path) as fh:
            for line in fh:
                if not line.startswith('>'):
                    continue
                header = line[1:].strip()
                names.append(header.split()[0] if header else '')
                for token in header.split():
                    if token.startswith('resistance_genes='):
                        val = token.split('=', 1)[1].strip()
                        for g in val.replace(';', ',').split(','):
                            if g and g not in ('-', 'None', 'NA'):
                                genes.add(g.strip())
    except OSError as exc:
        print(f"  WARN: cannot read {path}: {exc}", file=sys.stderr)
    return names, genes


def find_amrfinder(plasmid_stem, search_dir):
    """Resolve the AMRFinderPlus TSV for a plasmid file stem."""
    candidates = [plasmid_stem]
    if '_plasmid_' in plasmid_stem:
        candidates.append(plasmid_stem.split('_plasmid_')[0])
    for cand in candidates:
        path = os.path.join(search_dir, f"{cand}_amrfinder.tsv")
        if os.path.exists(path):
            return path
    # Last resort: any AMRFinder file whose prefix this stem starts with
    for fn in sorted(os.listdir(search_dir)):
        if fn.endswith('_amrfinder.tsv'):
            prefix = fn[:-len('_amrfinder.tsv')]
            if plasmid_stem.startswith(prefix):
                return os.path.join(search_dir, fn)
    return None


def load_amrfinder(path):
    """Return ``{contig: [gene symbols]}`` from an AMRFinderPlus TSV."""
    by_contig = {}
    try:
        with open(path) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            try:
                si = header.index('Element symbol')
                ci = header.index('Contig id')
                bi = header.index('Subtype')
            except ValueError:
                return by_contig
            try:
                pi = header.index('Scope')
            except ValueError:
                pi = -1
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if max(si, ci, bi) >= len(parts):
                    continue
                if parts[bi].strip().upper() != 'AMR':
                    continue
                # Scope 'plus' (efflux, stress response) is excluded everywhere the report is concerned, so a
                # plasmid carrying only efflux genes is not treated as AMR-carrying.
                if 0 <= pi < len(parts) and parts[pi].strip().lower() == 'plus':
                    continue
                symbol = parts[si].strip()
                contig = parts[ci].strip().split()[0] if parts[ci].strip() else ''
                if symbol:
                    by_contig.setdefault(contig, []).append(symbol)
    except OSError as exc:
        print(f"  WARN: cannot read {path}: {exc}", file=sys.stderr)
    return by_contig


def newick_escape(label):
    """Newick-safe tip label."""
    return label.replace('(', '_').replace(')', '_').replace(',', '_').replace(':', '_').replace(';', '_')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--plasmid-dir', required=True, help='Directory of plasmid FASTAs + AMRFinderPlus TSVs')
    ap.add_argument('--outdir', required=True, help='Where to copy AMR-carrying plasmid FASTAs')
    ap.add_argument('--summary', default=None, help='Optional TSV listing selected plasmids and genes')
    ap.add_argument('--min-genomes', type=int, default=3,
                    help='Minimum plasmids for mashtree to build a tree (default 3)')
    ap.add_argument('--fallback-newick', default=None,
                    help='Write a placeholder Newick here when too few plasmids qualify')
    args = ap.parse_args()

    if not os.path.isdir(args.plasmid_dir):
        print(f"ERROR: no such directory: {args.plasmid_dir}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    fastas = sorted(f for f in os.listdir(args.plasmid_dir) if f.endswith('.fasta'))
    print(f"Scanning {len(fastas)} plasmid FASTA(s) in {args.plasmid_dir}", file=sys.stderr)

    selected = []
    for fn in fastas:
        stem = fn[:-len('.fasta')]
        path = os.path.join(args.plasmid_dir, fn)
        contigs, header_genes = parse_fasta_contigs(path)

        genes = set(header_genes)
        af = find_amrfinder(stem, args.plasmid_dir)
        if af:
            by_contig = load_amrfinder(af)
            contig_set = set(contigs)
            for contig, glist in by_contig.items():
                # Blank contig field means the hit was not localised; only trust it when this plasmid actually
                # contains the named contig.
                if contig and contig in contig_set:
                    genes.update(g for g in glist if g)

        if genes:
            shutil.copy(path, os.path.join(args.outdir, fn))
            selected.append((stem, sorted(genes)))
            print(f"  + {stem}: {', '.join(sorted(genes))}", file=sys.stderr)
        else:
            print(f"  - {stem}: no AMR genes", file=sys.stderr)

    if args.summary:
        with open(args.summary, 'w') as out:
            out.write("plasmid\tn_amr_genes\tamr_genes\n")
            for stem, genes in selected:
                out.write(f"{stem}\t{len(genes)}\t{';'.join(genes)}\n")

    n = len(selected)
    print(f"Selected {n} AMR-carrying plasmid(s) of {len(fastas)}", file=sys.stderr)

    if n < args.min_genomes and args.fallback_newick:
        # mashtree needs several genomes; emit a valid placeholder so the assembly type still appears
        # downstream instead of vanishing.
        tips = [newick_escape(stem) for stem, _ in selected]
        if not tips:
            tips = ['no_AMR_plasmids']
        if len(tips) == 1:
            newick = f"({tips[0]});"
        else:
            newick = "(" + ",".join(tips) + ");"
        with open(args.fallback_newick, 'w') as out:
            out.write(newick + "\n")
        print(f"Too few plasmids for a tree ({n} < {args.min_genomes}); "
              f"wrote placeholder Newick with {len(tips)} tip(s)", file=sys.stderr)

    # Exit code communicates whether mashtree should run
    sys.exit(0 if n >= args.min_genomes else 3)


if __name__ == '__main__':
    main()