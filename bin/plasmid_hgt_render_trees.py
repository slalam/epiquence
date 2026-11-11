#!/usr/bin/env python3
"""Render per-cluster mashtree Newick files to SVG using toytree."""
import argparse, csv, os, glob, re, hashlib
import toytree
import toyplot
import toyplot.svg
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib.tree_svg import (BADGE_PALETTES, BADGE_NEUTRAL,
    _text_color, _build_color_map, _add_badges_to_svg)
from epiquence_lib.reporting import display_name, split_name



def load_metadata(filepath):
    """Load plasmid_hgt_metadata.tsv: plasmid_name -> {size_bp, resistance_genes}"""
    meta = {}
    if not filepath or not os.path.exists(filepath):
        return meta
    with open(filepath) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            name = row['plasmid_name']
            meta[name] = {
                'size_bp': int(row.get('size_bp', 0) or 0),
                'resistance_genes': row.get('resistance_genes', ''),
            }
    return meta


def load_species_map(gtdbtk_dir):
    """Load species from GTDB-Tk parsed TSVs: sample -> species"""
    species_map = {}
    if not gtdbtk_dir or not os.path.isdir(gtdbtk_dir):
        return species_map
    for f in os.listdir(gtdbtk_dir):
        if f.endswith('_gtdbtk_parsed.tsv'):
            with open(os.path.join(gtdbtk_dir, f)) as fh:
                for row in csv.DictReader(fh, delimiter='\t'):
                    sample = row.get('Sample', '').strip()
                    species = (row.get('species', '') or
                               row.get('Species', '') or
                               row.get('GTDB_Species', '')).strip()
                    if sample and species and species != '-':
                        species = re.sub(r'^s__', '', species)
                        species = re.sub(r'_[A-Z]$', '', species)
                        species = re.sub(r'_[A-Z] ', ' ', species)
                        species_map[sample] = species
    return species_map


def render_tree(nwk_file, metadata, species_map):
    """Render a Newick tree file to SVG with colored badge annotations."""
    with open(nwk_file) as f:
        nwk = f.read().strip()
    if not nwk:
        return None

    tree = toytree.tree(nwk)
    tree = tree.ladderize()
    raw_tips = tree.get_tip_labels()
    ntips = tree.ntips

    tip_data = []
    plain_labels = []

    for tip in raw_tips:
        clean = tip.replace('.fasta', '')
        if '/' in clean:
            clean = clean.split('/')[-1]

        sample, asm, mob_id = split_name(clean)
        disp = display_name(clean)

        entry = {'sample': disp, 'badges': []}

        sample_asm = f"{sample}_{asm}" if asm else sample
        sp = species_map.get(sample_asm, '') or species_map.get(sample, '')
        if sp:
            entry['badges'].append(('species', sp))

        meta = metadata.get(clean, {})
        genes = meta.get('resistance_genes', '')
        if genes:
            for gene in genes.replace(';', ',').split(','):
                gene = gene.strip()
                if gene:
                    entry['badges'].append(('amr', gene))

        tip_data.append(entry)
        plain_labels.append(disp)

    color_map = _build_color_map(tip_data)

    max_badge_width = 0
    for td in tip_data:
        bw = sum(len(val) * 6.2 + 14 for _, val in td['badges'])
        if bw > max_badge_width:
            max_badge_width = bw
    max_label_len = max((len(s) for s in plain_labels), default=10)
    height = max(200, ntips * 24)
    width = max(600, 350 + max_label_len * 8 + max_badge_width + 60)

    canvas, axes, mark = tree.draw(
        width=width, height=height,
        tip_labels=plain_labels,
        tip_labels_style={"font-size": "11px", "font-weight": "bold"},
        node_sizes=0,
        fixed_order=raw_tips,
    )

    svg_out = toyplot.svg.render(canvas)
    if hasattr(svg_out, 'tag'):
        import xml.etree.ElementTree as ET
        svg_str = ET.tostring(svg_out, encoding='unicode')
    elif isinstance(svg_out, bytes):
        svg_str = svg_out.decode('utf-8')
    else:
        svg_str = str(svg_out)

    return _add_badges_to_svg(svg_str, tip_data, color_map, width)


def main():
    parser = argparse.ArgumentParser(description='Render HGT cluster trees with toytree')
    parser.add_argument('--tree-files', nargs='+', required=True, help='Newick tree files from mashtree')
    parser.add_argument('--metadata', default=None, help='Plasmid HGT metadata TSV')
    parser.add_argument('--gtdbtk-dir', default=None, help='Dir with GTDB-Tk parsed TSVs')
    parser.add_argument('--outdir', default='.', help='Output directory for SVGs')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    metadata = load_metadata(args.metadata)
    species_map = load_species_map(args.gtdbtk_dir)
    print(f"Loaded metadata for {len(metadata)} plasmids")
    print(f"Loaded species for {len(species_map)} samples")
    if species_map:
        for k, v in list(species_map.items())[:5]:
            print(f"  {k} -> {v}")

    for nwk_file in args.tree_files:
        if not os.path.exists(nwk_file):
            continue
        cluster_id = os.path.basename(nwk_file).replace('.nwk', '')

        svg = render_tree(nwk_file, metadata, species_map)
        if svg:
            out_path = os.path.join(args.outdir, f'{cluster_id}_tree.svg')
            with open(out_path, 'w') as f:
                f.write(svg)
            print(f"Rendered tree: {cluster_id}")


if __name__ == '__main__':
    main()