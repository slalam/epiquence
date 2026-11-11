#!/usr/bin/env python3
"""Generate UMAP embedding from plasmid Mash distance matrix and visualize clusters using Vega-
Altair."""

import argparse
import csv
import json
import os
import re
import sys

import numpy as np
import umap
import altair as alt
import pandas as pd
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib.mash import read_matrix


def parse_mash_matrix(path):
    """Parse Mash distance matrix TSV."""
    row_names, col_names, rows = read_matrix(path, basename=True, strip_fasta=True)
    return row_names, np.array(rows)


def parse_clusters(path):
    """Parse plasmid_clusters.tsv -> {plasmid_name: cluster_id}."""
    clusters = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            clusters[row['plasmid']] = row['cluster_id']
    return clusters


def parse_plasmid_name(name):
    """Extract sample ID and MOB ID from plasmid name."""
    # Pattern: {sample}_{asm}_{plasmid}_{mob_id} or: {sample}_{asm}_plasmid_{mob_id}
    m = re.match(r'^(.+?)_(HL|HS|L|S)_(?:plasmid_)?(.+)$', name)
    if m:
        return m.group(1), m.group(2), m.group(3)
    return name, '', ''


def main():
    parser = argparse.ArgumentParser(description='UMAP embedding of plasmid Mash distances')
    parser.add_argument('--mash-matrix', required=True, help='Mash distance matrix TSV')
    parser.add_argument('--clusters', required=True, help='Plasmid clusters TSV')
    parser.add_argument('--asm-type', required=True, help='Assembly type (HL, HS, L, S)')
    parser.add_argument('--metadata', default=None, help='Plasmid HGT metadata TSV (sizes + resistance genes)')
    parser.add_argument('--n-neighbors', type=int, default=8, help='UMAP n_neighbors; lower = more local structure (default: 8)')
    parser.add_argument('--min-dist', type=float, default=0.01, help='UMAP min_dist; lower = tighter clusters (default: 0.01)')
    parser.add_argument('--spread', type=float, default=3.0, help='UMAP spread; higher = more space between clusters (default: 3.0)')
    parser.add_argument('--repulsion-strength', type=float, default=3.0, help='UMAP repulsion_strength; higher = more separation (default: 3.0)')
    parser.add_argument('--exclude-unclustered', action='store_true', help='Exclude unclustered plasmids from UMAP embedding')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for TSV and HTML')
    args = parser.parse_args()

    asm_type = args.asm_type

    names, dist_matrix = parse_mash_matrix(args.mash_matrix)
    n = len(names)

    if n < 4:
        print(f"Too few plasmids ({n}) for UMAP embedding, skipping")
        with open(f'{args.output_prefix}_umap_embedding.tsv', 'w') as f:
            f.write('plasmid\tUMAP1\tUMAP2\tcluster_id\tsample\tmob_id\n')
        with open(f'{args.output_prefix}_umap_plot.html', 'w') as f:
            f.write('<p style="color:#999;font-style:italic;">Too few plasmids for UMAP embedding.</p>')
        return

    cluster_map = parse_clusters(args.clusters)

    # Parse metadata for sizes and resistance genes
    meta_map = {}
    if args.metadata and os.path.exists(args.metadata):
        with open(args.metadata) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                pname = row.get('plasmid_name', row.get('plasmid', ''))
                if not pname:
                    continue
                try:
                    size_bp = int(row.get('size_bp', 0) or 0)
                except (ValueError, TypeError):
                    size_bp = 0
                meta_map[pname] = {
                    'size_kb': round(size_bp / 1000, 1) if size_bp else 0,
                    'resistance_genes': row.get('resistance_genes', ''),
                    'species': row.get('species', ''),
                }
        print(f"Loaded metadata for {len(meta_map)} plasmids")
        if meta_map:
            first_key = next(iter(meta_map))
            print(f"  Example key: '{first_key}' -> {meta_map[first_key]}")
        # Show first few Mash matrix names for comparison
        if names:
            print(f"  First Mash names: {names[:3]}")
            matched = sum(1 for n in names if n in meta_map)
            print(f"  Matched: {matched}/{len(names)}")

    dist_matrix = (dist_matrix + dist_matrix.T) / 2
    np.fill_diagonal(dist_matrix, 0)

    # Optionally filter out unclustered plasmids before embedding
    clustered_set = set(cluster_map.keys())
    if args.exclude_unclustered:
        keep_idx = [i for i, name in enumerate(names) if name in clustered_set]
        if len(keep_idx) < 4:
            print(f"Too few clustered plasmids ({len(keep_idx)}) for UMAP, including all")
        else:
            names = [names[i] for i in keep_idx]
            dist_matrix = dist_matrix[np.ix_(keep_idx, keep_idx)]
            n = len(names)
            print(f"Filtered to {n} clustered plasmids (excluded {len(keep_idx)} unclustered)")

    # Semi-supervised UMAP: cluster labels pull members together (target_weight balances labels vs distance).
    n_neighbors = min(args.n_neighbors, n - 1)

    # Build numeric cluster labels for supervised component.
    label_map = {}
    label_counter = 0
    target_labels = np.full(n, -1, dtype=int)
    for i, name in enumerate(names):
        craw = cluster_map.get(name, 'unclustered')
        if craw != 'unclustered':
            if craw not in label_map:
                label_map[craw] = label_counter
                label_counter += 1
            target_labels[i] = label_map[craw]

    # Use semi-supervised UMAP if we have cluster labels
    has_labels = (target_labels >= 0).sum() > 0
    if has_labels:
        print(f"Using semi-supervised UMAP with {label_counter} cluster labels (target_weight=0.4)")
        reducer = umap.UMAP(
            n_components=2,
            metric='precomputed',
            n_neighbors=n_neighbors,
            min_dist=args.min_dist,
            spread=args.spread,
            repulsion_strength=args.repulsion_strength,
            target_weight=0.4,
            random_state=42,
        )
        embedding = reducer.fit_transform(dist_matrix, y=target_labels)
    else:
        reducer = umap.UMAP(
            n_components=2,
            metric='precomputed',
            n_neighbors=n_neighbors,
            min_dist=args.min_dist,
            spread=args.spread,
            repulsion_strength=args.repulsion_strength,
            random_state=42,
        )
        embedding = reducer.fit_transform(dist_matrix)

    records = []
    for i, name in enumerate(names):
        sample, asm, mob_id = parse_plasmid_name(name)
        cluster_raw = cluster_map.get(name, 'unclustered')
        if cluster_raw != 'unclustered':
            # Extract number from hgt_cluster_1 format
            cnum = cluster_raw.replace('hgt_cluster_', '')
            cluster_label = f'{asm_type}_C{cnum}'
        else:
            cluster_label = 'unclustered'
        meta = meta_map.get(name, {})
        records.append({
            'plasmid': name,
            'UMAP1': round(embedding[i, 0], 4),
            'UMAP2': round(embedding[i, 1], 4),
            'cluster_id': cluster_label,
            'sample': sample,
            'mob_id': mob_id,
            'species': meta.get('species', ''),
            'size_kb': meta.get('size_kb', 0),
            'resistance_genes': meta.get('resistance_genes', ''),
        })

    df = pd.DataFrame(records)

    # Build legend labels with MOB IDs: "HL_C1 (pAA102, pAA364)"
    cluster_mob_map = {}  # cluster_id -> set of mob_ids
    for rec in records:
        cid = rec['cluster_id']
        if cid != 'unclustered' and rec['mob_id']:
            cluster_mob_map.setdefault(cid, set()).add(rec['mob_id'])

    def make_legend_label(row):
        cid = row['cluster_id']
        if cid == 'unclustered':
            return 'unclustered'
        mobs = cluster_mob_map.get(cid, set())
        if mobs:
            mob_str = ', '.join(f'p{m}' for m in sorted(mobs))
            return f'{cid} ({mob_str})'
        return cid

    df['cluster_legend'] = df.apply(make_legend_label, axis=1)

    tsv_path = f'{args.output_prefix}_umap_embedding.tsv'
    df.to_csv(tsv_path, sep='\t', index=False)
    print(f"Saved embedding: {tsv_path} ({len(df)} plasmids)")

    # ── Generate Altair plot ──
    # Separate clustered and unclustered
    df_clustered = df[df['cluster_id'] != 'unclustered'].copy()
    df_unclustered = df[df['cluster_id'] == 'unclustered'].copy()
    n_clustered = len(df_clustered)
    n_unclustered = len(df_unclustered)

    # Cluster color palette: use cluster_legend for display
    legend_labels = sorted(df_clustered['cluster_legend'].unique(),
                           key=lambda x: (x.split('_C')[0], int(x.split('_C')[-1].split(' ')[0]) if x.split('_C')[-1].split(' ')[0].isdigit() else 0))
    n_clusters = len(legend_labels)

    palette = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#bcbd22', '#17becf', '#aec7e8',
        '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94',
        '#f7b6d2', '#dbdb8d', '#9edae5', '#393b79', '#637939',
        '#8c6d31', '#843c39', '#7b4173', '#5254a3', '#6b6ecf',
    ]
    color_domain = legend_labels
    color_range = [palette[i % len(palette)] for i in range(len(legend_labels))]

    has_sizes = df['size_kb'].sum() > 0

    selection = alt.selection_point(fields=['cluster_legend'], bind='legend')

    tooltips = [
        alt.Tooltip('plasmid:N', title='Plasmid'),
        alt.Tooltip('sample:N', title='Sample'),
        alt.Tooltip('species:N', title='Species'),
        alt.Tooltip('cluster_id:N', title='Cluster'),
        alt.Tooltip('mob_id:N', title='MOB ID'),
        alt.Tooltip('size_kb:Q', title='Size (Kb)', format='.1f'),
        alt.Tooltip('resistance_genes:N', title='AMR Genes'),
    ]

    # Layer 1: Unclustered plasmids (background: small, grey, low opacity).
    all_umap = df[['UMAP1', 'UMAP2']]
    x_min, x_max = all_umap['UMAP1'].min(), all_umap['UMAP1'].max()
    y_min, y_max = all_umap['UMAP2'].min(), all_umap['UMAP2'].max()
    x_pad = (x_max - x_min) * 0.06
    y_pad = (y_max - y_min) * 0.06
    x_domain = [x_min - x_pad, x_max + x_pad]
    y_domain = [y_min - y_pad, y_max + y_pad]

    x_enc = alt.X('UMAP1:Q', title='UMAP 1',
                   scale=alt.Scale(domain=x_domain),
                   axis=alt.Axis(grid=True, gridOpacity=0.15))
    y_enc = alt.Y('UMAP2:Q', title='UMAP 2',
                   scale=alt.Scale(domain=y_domain),
                   axis=alt.Axis(grid=True, gridOpacity=0.15))

    layers = []
    if not df_unclustered.empty:
        unclustered_layer = alt.Chart(df_unclustered).mark_circle(
            size=25,
            color='#d0d0d0',
            stroke='#bbb',
            strokeWidth=0.3,
            opacity=0.35,
        ).encode(
            x=x_enc,
            y=y_enc,
            tooltip=tooltips,
        )
        layers.append(unclustered_layer)

    # Layer 2: Clustered plasmids (foreground - colored, prominent)
    if not df_clustered.empty:
        size_enc = (alt.Size('size_kb:Q', title='Size (Kb)',
                            scale=alt.Scale(range=[50, 400]),
                            legend=alt.Legend(symbolSize=100))
                    if has_sizes else alt.value(90))

        clustered_layer = alt.Chart(df_clustered).mark_circle(
            stroke='#333',
            strokeWidth=0.5,
        ).encode(
            x=x_enc,
            y=y_enc,
            color=alt.Color('cluster_legend:N', title='Cluster',
                            scale=alt.Scale(domain=color_domain, range=color_range),
                            legend=alt.Legend(columns=2, symbolSize=100)),
            size=size_enc,
            tooltip=tooltips,
            opacity=alt.condition(selection, alt.value(0.9), alt.value(0.12)),
        ).add_params(selection)
        layers.append(clustered_layer)

    _asm_labels = {'S': 'Short Read', 'L': 'Long Read', 'HL': 'Hybrid (Long)', 'HS': 'Hybrid (Short)'}
    _asm_label = _asm_labels.get(asm_type, asm_type)

    chart = alt.layer(*layers).properties(
        width=650,
        height=480,
        title=f'{_asm_label} - Plasmid UMAP ({n_clustered} clustered, {n_unclustered} unclustered, {n_clusters} clusters)'
    ).configure_view(
        strokeWidth=0,
    ).configure(
        font='Inter, Helvetica Neue, sans-serif',
        background='#fafafa',
    )

    html_path = f'{args.output_prefix}_umap_plot.html'
    chart.save(html_path)

    # Also save a minimal embeddable version with PNG download button
    spec = chart.to_json(indent=None)
    embed_path = f'{args.output_prefix}_umap_embed.html'
    with open(embed_path, 'w') as f:
        f.write('<!DOCTYPE html><html><head>')
        f.write('<script src="https://cdn.jsdelivr.net/npm/vega@5"></script>')
        f.write('<script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>')
        f.write('<script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>')
        f.write('<style>body{margin:0;background:#fafafa;font-family:sans-serif;}'
                '.dl-btn{position:absolute;top:8px;right:8px;display:inline-flex;align-items:center;'
                'justify-content:center;width:28px;height:28px;border-radius:4px;background:#f0f0f0;'
                'color:#555;text-decoration:none;border:1px solid #ddd;cursor:pointer;z-index:10;}'
                '.dl-btn:hover{background:#e0e0e0;}'
                '.dl-btn svg{width:16px;height:16px;}</style>')
        f.write('</head><body>')
        f.write('<div style="position:relative;">')
        f.write('<div id="vis"></div>')
        f.write('<button class="dl-btn" id="savePng" title="Download as PNG">'
                '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" '
                'stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">'
                '<path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/>'
                '<polyline points="7 10 12 15 17 10"/>'
                '<line x1="12" y1="15" x2="12" y2="3"/>'
                '</svg></button>')
        f.write('</div>')
        f.write(f'<script>'
                f'vegaEmbed("#vis",{spec},{{actions:false}}).then(function(result){{'
                f'document.getElementById("savePng").addEventListener("click",function(){{'
                f'result.view.toImageURL("png",2).then(function(url){{'
                f'var a=document.createElement("a");a.href=url;'
                f'a.download="{args.output_prefix}_umap.png";a.click();'
                f'}});}});}});'
                f'</script>')
        f.write('</body></html>')

    print(f"Saved plot: {html_path}")
    print(f"Saved embed: {embed_path}")
    print(f"  Clusters: {n_clusters} + {'unclustered' if n_unclustered > 0 else 'none unclustered'}")


if __name__ == '__main__':
    main()