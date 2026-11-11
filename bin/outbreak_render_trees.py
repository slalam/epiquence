#!/usr/bin/env python3
"""Render Newick trees to SVG using toytree."""
import argparse, csv, os, glob, re, hashlib
import toytree
import toyplot
import toyplot.svg
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib.tree_svg import (BADGE_PALETTES, BADGE_NEUTRAL,
    _text_color, _build_color_map, _add_badges_to_svg)



def render_tree_svg(nwk_file, metadata=None, species_map=None, mlst_map=None, amr_map=None):
    """Render Newick tree to SVG string with colored metadata badges on leaf labels."""
    with open(nwk_file) as f:
        nwk = f.read().strip()
    if not nwk:
        return "<p><i>Empty tree file</i></p>"

    tree = toytree.tree(nwk)
    tree = tree.ladderize()
    ntips = tree.ntips
    raw_tips = tree.get_tip_labels()

    tip_data = []
    plain_labels = []

    for tip in raw_tips:
        sample = tip.replace('.fasta', '')
        if '/' in sample:
            sample = sample.split('/')[-1]
        base_id = re.sub(r'_(L|S|HL|HS)$', '', sample)

        entry = {'sample': sample, 'badges': []}

        if species_map:
            sp = species_map.get(sample, '') or species_map.get(base_id, '')
            if sp:
                entry['badges'].append(('species', sp))

        if mlst_map:
            st = mlst_map.get(sample, '') or mlst_map.get(base_id, '')
            if st and st != '-':
                entry['badges'].append(('mlst', f'ST{st}'))

        if amr_map:
            genes = amr_map.get(sample, '') or amr_map.get(base_id, '')
            if genes and genes != '-':
                for gene in genes.replace(';', ',').split(','):
                    gene = gene.strip()
                    if gene:
                        entry['badges'].append(('amr', gene))

        if metadata:
            meta = metadata.get(sample, {}) or metadata.get(base_id, {})
            for k, v in meta.items():
                if k.lower() not in ('sample_id', 'sample') and v:
                    entry['badges'].append(('meta', str(v)))

        tip_data.append(entry)
        plain_labels.append(sample)

    color_map = _build_color_map(tip_data)

    max_badge_width = 0
    for td in tip_data:
        bw = sum(len(val) * 6.2 + 14 for _, val in td['badges'])
        if bw > max_badge_width:
            max_badge_width = bw
    max_label_len = max((len(s) for s in plain_labels), default=10)
    width = max(600, 350 + max_label_len * 8 + max_badge_width + 60)
    height = max(300, ntips * 24)

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


def parse_metadata(metadata_file):
    """Parse metadata CSV/TSV."""
    metadata = {}
    if not metadata_file or not os.path.exists(metadata_file):
        return metadata
    sep = '\t' if metadata_file.lower().endswith('.tsv') else ','
    with open(metadata_file) as f:
        reader = csv.DictReader(f, delimiter=sep)
        for row in reader:
            sid = row.get('sample_id', row.get('Sample', row.get('sample', '')))
            if sid:
                metadata[sid.strip()] = {k: v for k, v in row.items()}
    return metadata


def main():
    parser = argparse.ArgumentParser(description='Render outbreak trees with toytree')
    parser.add_argument('--outdir', required=True, help='Output directory for rendered tree SVGs')
    parser.add_argument('--cluster-trees', nargs='*', default=[], help='Newick tree files')
    parser.add_argument('--overview-tree', default=None, help='Overview Newick tree')
    parser.add_argument('--metadata', default=None, help='Metadata CSV/TSV')
    parser.add_argument('--outbreak-summary', default=None, help='Outbreak summary JSON')
    parser.add_argument('--mlst-dir', default=None, help='MLST results directory')
    parser.add_argument('--chromosomal-amr', default=None, help='Chromosomal AMR summary TSV')
    parser.add_argument('--gtdbtk-dir', default=None, help='Dir with GTDB-Tk parsed TSVs')
    parser.add_argument('--prefix', default='', help='Prefix for output SVG filenames')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    metadata = parse_metadata(args.metadata)

    mlst_map = {}
    if args.mlst_dir and os.path.isdir(args.mlst_dir):
        for f in os.listdir(args.mlst_dir):
            if f.endswith('_mlst.tsv') or f.endswith('_mlst_parsed.tsv'):
                with open(os.path.join(args.mlst_dir, f)) as fh:
                    for row in csv.DictReader(fh, delimiter='\t'):
                        sample = row.get('Sample', '').strip()
                        st = (row.get('MLST_ST', '') or row.get('ST', '')).strip()
                        if sample and st and st != '-':
                            mlst_map[sample] = st

    amr_map = {}
    if args.chromosomal_amr and os.path.exists(args.chromosomal_amr):
        with open(args.chromosomal_amr) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                sample = row.get('Sample', '').strip()
                genes = row.get('Chromosomal_AMR_Genes', '').strip()
                if sample and genes:
                    amr_map[sample] = genes

    species_map = {}
    if args.gtdbtk_dir and os.path.isdir(args.gtdbtk_dir):
        for f in os.listdir(args.gtdbtk_dir):
            if f.endswith('_gtdbtk_parsed.tsv'):
                with open(os.path.join(args.gtdbtk_dir, f)) as fh:
                    for row in csv.DictReader(fh, delimiter='\t'):
                        sample = row.get('Sample', '').strip()
                        species = row.get('species', '').strip()
                        genus = row.get('genus', '').strip()
                        # Remove GTDB polyphyletic suffixes
                        if species:
                            species = re.sub(r'_[A-Z]$', '', species)
                            species = re.sub(r'_[A-Z] ', ' ', species)
                        if genus:
                            genus = re.sub(r'_[A-Z]$', '', genus)
                        if sample and species:
                            species_map[sample] = species
                        elif sample and genus:
                            species_map[sample] = f"{genus} sp."

    # Overlay with outbreak summary
    if args.outbreak_summary and os.path.exists(args.outbreak_summary):
        import json
        with open(args.outbreak_summary) as f:
            summary = json.load(f)
        for cluster_name, info in summary.get('clusters', {}).items():
            sp = info.get('species', '')
            if sp:
                sp = re.sub(r'_[A-Z]$', '', sp)
                sp = re.sub(r'_[A-Z] ', ' ', sp)
                for member in info.get('samples', []):
                    if member not in species_map:
                        species_map[member] = sp

    print(f"Species map: {len(species_map)} samples, MLST map: {len(mlst_map)}, AMR map: {len(amr_map)}")

    for nwk_file in (args.cluster_trees or []):
        if not os.path.exists(nwk_file):
            continue
        cluster_name = os.path.basename(nwk_file).replace('.treefile', '').replace('_iqtree', '')
        # Strip trailing prefix identifier to avoid duplication (e.g., _cga when prefix is cga_)
        if args.prefix:
            suffix = '_' + args.prefix.rstrip('_')
            if cluster_name.endswith(suffix):
                cluster_name = cluster_name[:-len(suffix)]
        svg = render_tree_svg(nwk_file, metadata=metadata, species_map=species_map,
                              mlst_map=mlst_map, amr_map=amr_map)
        if svg:
            out_path = os.path.join(args.outdir, f'{args.prefix}{cluster_name}_tree.svg')
            with open(out_path, 'w') as f:
                f.write(svg)
            print(f"Rendered: {cluster_name}")


if __name__ == '__main__':
    main()