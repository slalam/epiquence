#!/usr/bin/env python3
"""Build outbreak report HTML for MultiQC custom content."""
import argparse, csv, json, os, glob, base64
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from scipy.sparse.csgraph import minimum_spanning_tree
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib.reporting import make_download_button
from epiquence_lib.heatmap import render_heatmap_svg, snp_color


def render_snp_heatmap_html(snp_matrix_file):
    """Generate HTML heatmap from SNP distance matrix with hierarchical clustering."""
    samples = []
    matrix = []
    with open(snp_matrix_file) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        for row in reader:
            name = os.path.basename(row[0]).replace('.fasta', '')
            samples.append(name)
            matrix.append([int(x) for x in row[1:]])

    n = len(samples)
    if n == 0:
        return "<p><i>Empty SNP matrix</i></p>", [], []

    mat = np.array(matrix)
    max_val = mat.max() if mat.max() > 0 else 1

    if n > 2:
        np.fill_diagonal(mat, 0)
        sym = (mat + mat.T) / 2
        condensed = squareform(sym)
        Z = linkage(condensed, method='average')
        order = list(leaves_list(Z))
    else:
        order = list(range(n))

    ordered_samples = [samples[i] for i in order]
    ordered_matrix = mat[np.ix_(order, order)]

    cell_size = 40
    html = '<div style="overflow-x:auto;margin:10px 0;">'
    html += '<table style="border-collapse:collapse;font-family:monospace;font-size:11px;">'

    html += '<tr><td style="min-width:120px;"></td>'
    for s in ordered_samples:
        label = s[:12] + '..' if len(s) > 14 else s
        html += f'<td style="writing-mode:vertical-rl;text-align:left;padding:4px 2px;font-weight:600;max-width:{cell_size}px;">{label}</td>'
    html += '</tr>'

    for i, row_sample in enumerate(ordered_samples):
        html += f'<tr><td style="padding:2px 8px;font-weight:600;white-space:nowrap;">{row_sample}</td>'
        for j in range(len(ordered_samples)):
            val = int(ordered_matrix[i][j])
            if i == j:
                bg = '#e8e8e8'
                txt = '#999'
            else:
                # Dual-scale coloring for better resolution at low SNP distances
                if val <= 100:
                    # 0 SNPs = dark green (#1a7a2e), 100 SNPs = yellow (#e6c700)
                    t = val / 100.0
                    r = int(26 + t * (230 - 26))
                    g = int(122 + t * (199 - 122))
                    b = int(46 + t * (0 - 46))
                    bg = f'rgb({r},{g},{b})'
                    txt = '#fff' if t < 0.5 else '#333'
                else:
                    # 101 SNPs = yellow-orange, max SNPs = dark red
                    t = min(1.0, (val - 100) / max(max_val - 100, 1))
                    r = int(230 + t * (180 - 230))
                    g = int(180 - t * 150)
                    b = int(0 + t * 30)
                    bg = f'rgb({r},{g},{b})'
                    txt = '#fff' if t > 0.3 else '#333'
            html += f'<td style="background:{bg};color:{txt};text-align:center;padding:2px;min-width:{cell_size}px;height:{cell_size}px;font-size:10px;">{val}</td>'
        html += '</tr>'

    html += '</table></div>'
    return html, ordered_samples, ordered_matrix


def build_clustered_tsv(ordered_samples, ordered_matrix, species='', mlst_map=None, amr_map=None, metadata=None):
    """Build TSV string from clustered SNP matrix with annotation columns."""
    annot_cols = ['Isolate_ID', 'Species', 'MLST', 'Resistance_Genes']
    meta_cols = []
    if metadata:
        for s in ordered_samples:
            base_id = s.rsplit('_', 1)[0] if any(s.endswith(sfx) for sfx in ('_L', '_S', '_HL', '_HS')) else s
            meta = metadata.get(s) or metadata.get(base_id) or {}
            for k in meta:
                if k not in meta_cols:
                    meta_cols.append(k)

    header = annot_cols + meta_cols + ['Sample_ID'] + ordered_samples
    lines = ['\t'.join(header)]

    for i, sample in enumerate(ordered_samples):
        # Strip assembly suffix for isolate ID
        isolate_id = sample
        for sfx in ('_HL', '_HS', '_L', '_S'):
            if sample.endswith(sfx):
                isolate_id = sample[:-len(sfx)]
                break

        st = (mlst_map.get(sample, '') or mlst_map.get(isolate_id, '')) if mlst_map else ''
        amr = (amr_map.get(sample, '') or amr_map.get(isolate_id, '')) if amr_map else ''

        meta_vals = []
        if meta_cols and metadata:
            meta = metadata.get(sample) or metadata.get(isolate_id) or {}
            meta_vals = [meta.get(k, '') for k in meta_cols]

        row_vals = [str(int(ordered_matrix[i][j])) for j in range(len(ordered_samples))]
        line = [isolate_id, species, st, amr] + meta_vals + [sample] + row_vals
        lines.append('\t'.join(line))
    return '\n'.join(lines)


def build_heatmap_svg(ordered_samples, ordered_matrix):
    """Standalone SVG of the SNP distance heatmap with dual-scale colouring."""
    if len(ordered_samples) == 0:
        return ''
    return render_heatmap_svg(
        ordered_samples, ordered_matrix,
        color_fn=snp_color,
        cell_size=40,
        label_angle=-60,
    )

def render_mst_d3(snp_matrix_file, cluster_name, mlst_map=None, amr_map=None, metadata=None):
    """Build MST and render interactive MST with dynamic collapse threshold slider."""
    samples = []
    matrix = []
    with open(snp_matrix_file) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        for row in reader:
            name = os.path.basename(row[0]).replace('.fasta', '')
            samples.append(name)
            matrix.append([int(x) for x in row[1:]])

    n = len(samples)
    if n < 2:
        return ''

    mat = np.array(matrix, dtype=float)
    np.fill_diagonal(mat, 0)
    sym = (mat + mat.T) / 2

    eps = 0.001
    mst_input = sym.copy()
    zero_mask = (mst_input == 0)
    np.fill_diagonal(zero_mask, False)
    mst_input[zero_mask] = eps
    mst = minimum_spanning_tree(mst_input).toarray()

    # Build MST edges on individual samples
    mst_edges = []
    for i in range(n):
        for j in range(i + 1, n):
            w = max(mst[i][j], mst[j][i])
            if w > 0:
                mst_edges.append({'source': i, 'target': j, 'snps': int(round(sym[i][j]))})

    import re as _re
    meta_field_names = []
    if metadata:
        for s in samples:
            base_id = _re.sub(r'_(L|S|HL|HS)$', '', s)
            meta = metadata.get(s, {}) or metadata.get(base_id, {})
            for k in meta:
                if k.lower() not in ('sample_id', 'sample') and k not in meta_field_names:
                    meta_field_names.append(k)

    sample_data = []
    for i, s in enumerate(samples):
        base_id = _re.sub(r'_(L|S|HL|HS)$', '', s)
        info = {'id': i, 'name': s}
        info['st'] = (mlst_map.get(s, '') or mlst_map.get(base_id, '')) if mlst_map else ''
        info['amr'] = (amr_map.get(s, '') or amr_map.get(base_id, '')) if amr_map else ''
        info['metaFields'] = {}
        if metadata:
            meta = metadata.get(s, {}) or metadata.get(base_id, {})
            for k in meta_field_names:
                info['metaFields'][k] = meta.get(k, '')
        sample_data.append(info)

    # Symmetric distance matrix as flat array for JS
    dist_flat = sym.flatten().tolist()

    color_by_fields = ['MLST']
    if any(sd['amr'] for sd in sample_data):
        color_by_fields.append('Chromosomal AMR')
    color_by_fields.extend(meta_field_names)

    samples_json = json.dumps(sample_data)
    edges_json = json.dumps(mst_edges)
    dist_json = json.dumps(dist_flat)
    color_fields_json = json.dumps(color_by_fields)
    div_id = f'mst_{cluster_name}'.replace('-', '_').replace('.', '_')
    h = 500 if n <= 4 else 800

    _tmpl = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates', 'mst_d3.html')
    with open(_tmpl) as _tf:
        html = _tf.read()
    for _tok, _val in (("@@EPIQ_div_id@@", div_id), ("@@EPIQ_cluster_name@@", cluster_name),
                       ("@@EPIQ_h@@", h), ("@@EPIQ_samples_json@@", samples_json),
                       ("@@EPIQ_edges_json@@", edges_json), ("@@EPIQ_dist_json@@", dist_json),
                       ("@@EPIQ_color_fields_json@@", color_fields_json)):
        html = html.replace(_tok, str(_val))
    return html


def main():
    parser = argparse.ArgumentParser(description='Build outbreak report for MultiQC')
    parser.add_argument('--tree-svgs', required=True, help='Directory with rendered SVG tree files')
    parser.add_argument('--cluster-results', required=True, help='Directory with per-cluster SNP matrices and trees')
    parser.add_argument('--outbreak-summary', required=True, help='Outbreak summary JSON')
    parser.add_argument('--ref-dir', default=None, help='Directory with parsnp_reference.txt files')
    parser.add_argument('--mlst-dir', default=None, help='Directory with MLST parsed TSV files')
    parser.add_argument('--chromosomal-amr', default=None, help='Chromosomal AMR TSV file')
    parser.add_argument('--metadata', default=None, help='Metadata CSV/TSV file')
    parser.add_argument('--cga-tree-svgs', default=None, help='Directory with CGA rendered SVG tree files')
    parser.add_argument('--cga-cluster-results', default=None, help='Directory with CGA per-cluster SNP matrices and trees')
    parser.add_argument('--output', default='outbreak_section_mqc.yaml', help='Output YAML for MultiQC')
    parser.add_argument('--workflow-id', default='')
    parser.add_argument('--workflow-timestamp', default='')
    parser.add_argument('--group-id', default='')
    args = parser.parse_args()

    summary = {}
    if os.path.exists(args.outbreak_summary):
        with open(args.outbreak_summary) as f:
            summary = json.load(f)

    # Load MLST data: sample -> ST
    mlst_map = {}  # sample_name -> "scheme ST{n}"
    if args.mlst_dir and os.path.isdir(args.mlst_dir):
        for mf in os.listdir(args.mlst_dir):
            if not mf.endswith('.tsv'):
                continue
            try:
                with open(os.path.join(args.mlst_dir, mf)) as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        sample = row.get('Sample', row.get('sample', '')).strip()
                        st = row.get('MLST_ST', row.get('ST', row.get('st', ''))).strip()
                        scheme = row.get('MLST_Scheme', row.get('Scheme', row.get('scheme', ''))).strip()
                        if sample and st and st != '-':
                            mlst_map[sample] = f"ST{st}"
            except Exception:
                pass

    # Load chromosomal AMR data
    amr_map = {}
    if args.chromosomal_amr and os.path.exists(args.chromosomal_amr):
        with open(args.chromosomal_amr) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                sample = row.get('Sample', '').strip()
                genes = row.get('Chromosomal_AMR_Genes', '').strip()
                if sample and genes and genes != '-':
                    amr_map[sample] = genes

    meta_map = {}
    if args.metadata and os.path.exists(args.metadata):
        sep = '\t' if args.metadata.lower().endswith('.tsv') else ','
        with open(args.metadata) as f:
            reader = csv.DictReader(f, delimiter=sep)
            for row in reader:
                sid = row.get('sample_id', row.get('Sample', row.get('sample', '')))
                if sid:
                    meta_map[sid.strip()] = {k: v for k, v in row.items() if k.lower() not in ('sample_id', 'sample')}

    # Load reference sample IDs from parsnp_reference.txt files
    ref_map = {}  # cluster_name -> reference sample ID
    if args.ref_dir and os.path.isdir(args.ref_dir):
        for rf in glob.glob(os.path.join(args.ref_dir, '*_reference.txt')):
            # Filename: parsnp_reference.txt (all have same name, but content differs)
            try:
                with open(rf) as f:
                    ref_sample = f.read().strip()
                if ref_sample:
                    # Map all clusters to their reference (files staged flat)
                    ref_map[rf] = ref_sample
            except Exception:
                pass
    # If all reference files have the same name, try reading all in the dir
    if not ref_map and args.ref_dir and os.path.isdir(args.ref_dir):
        for rf in os.listdir(args.ref_dir):
            fp = os.path.join(args.ref_dir, rf)
            if os.path.isfile(fp):
                try:
                    with open(fp) as f:
                        ref_sample = f.read().strip()
                    if ref_sample:
                        ref_map[fp] = ref_sample
                except Exception:
                    pass
    # Build a simple lookup: collect all reference IDs
    all_refs = list(ref_map.values())

    html_sections = []
    # html_sections always receives the accordion CSS below, so its emptiness is not a usable
    # "were there clusters" signal.
    n_cluster_sections = 0

    html_sections.append('''<style>
    details summary::-webkit-details-marker { display: none; }
    details summary { list-style: none; }
    details[open] > summary span:first-child { transform: rotate(90deg); display: inline-block; }
    details summary span:first-child { transition: transform 0.2s; display: inline-block; }
    details[open] { border-color: #adb5bd !important; }
    </style>''')

    # Cluster summary table: write as native MultiQC table YAML
    clusters = summary.get('clusters', {})
    outbreak_summary_table = {}
    if clusters:
        for cluster_name in sorted(clusters.keys()):
            info = clusters[cluster_name]
            species = info.get('species', 'Unknown')
            members = info.get('samples', [])
            cluster_num = cluster_name.split('_')[-1]

            sts = set()
            for m in members:
                st = mlst_map.get(m, '')
                if st:
                    sts.add(st)
            st_str = '; '.join(sorted(sts)) if sts else '-'
            sample_str = '; '.join(sorted(members))

            row_key = f'{species.replace(" ", "_")}_cluster_{cluster_num}'
            outbreak_summary_table[row_key] = {
                'Species': species,
                'Cluster': cluster_num,
                'MLST': st_str,
                'Samples': len(members),
                'Sample_IDs': sample_str,
            }

    def yaml_escape(val):
        """Escape a value for safe YAML output."""
        s = str(val)
        if any(c in s for c in ":|,[]{}#&*!>?'\"\\@`%"):
            return '"' + s.replace('\\', '\\\\').replace('"', '\\"') + '"'
        return s

    # Write outbreak cluster summary as native MultiQC table.
    summary_yaml_path = os.path.join(os.path.dirname(args.output), 'outbreak_cluster_summary_mqc.yaml')
    import json as _json_yaml

    if not outbreak_summary_table:
        open(summary_yaml_path, 'w').close()
        print('No outbreak clusters: Species Cluster Summary section suppressed')
    else:
      yaml_data = {
        'id': 'outbreak_cluster_summary',
        'section_name': 'Species Cluster Summary',
        'description': '',
        'plot_type': 'table',
        'pconfig': {
            'id': 'outbreak_cluster_summary_table',
            'title': 'Species Cluster Summary',
            'sortRows': True,
            'col1_header': 'Cluster ID',
        },
        'headers': {
            'Species': {'title': 'Species', 'scale': False},
            'Cluster': {'title': 'Cluster', 'scale': False},
            'MLST': {'title': 'MLST', 'scale': False},
            'Samples': {'title': 'Samples', 'format': '{:,.0f}', 'scale': 'Blues'},
            'Sample_IDs': {'title': 'Sample IDs', 'scale': False},
        },
        'data': outbreak_summary_table,
    }
    
      # Use json to write valid YAML (JSON is valid YAML)
      with open(summary_yaml_path, 'w') as f:
          _json_yaml.dump(yaml_data, f, indent=2, default=str)
      print(f"Outbreak cluster summary table written to {summary_yaml_path}")

    # Per-cluster results: iterate ALL clusters from summary, not just those with result dirs
    for cluster_name in sorted(clusters.keys()):
        cluster_info = clusters[cluster_name]
        species = cluster_info.get('species', 'Unknown')
        n_samples = cluster_info.get('size', '?')
        cluster_num = cluster_name.split("_")[-1]

        # Find reference for this cluster
        cluster_members = cluster_info.get('samples', [])
        ref_id = ''
        for r in all_refs:
            if r in cluster_members:
                ref_id = r
                break
        if not ref_id and cluster_members:
            import re as _re
            for suffix in ['_HL', '_HS', '_L', '_S']:
                for m in sorted(cluster_members):
                    if m.endswith(suffix):
                        ref_id = m
                        break
                if ref_id:
                    break
            if not ref_id:
                ref_id = sorted(cluster_members)[0]

        # Get MLST for header
        sts = set()
        for m in cluster_members:
            st = mlst_map.get(m, '')
            if st:
                sts.add(st)
        st_label = f' <span style="font-size:11px;color:#999;font-weight:normal;">{"; ".join(sorted(sts))}</span>' if sts else ""

        section_html = (
            f'<details style="margin:8px 0;border:1px solid #dee2e6;border-radius:6px;overflow:hidden;">'
            f'<summary style="padding:10px 14px;background:#f8f9fa;cursor:pointer;font-size:14px;font-weight:600;list-style:none;">'
            f'<span style="margin-right:8px;">&#9654;</span>'
            f'{species} | Cluster {cluster_num} ({n_samples} samples){st_label}'
            f'</summary>'
            f'<div style="padding:12px 16px;">'
        )

        # SNP heatmap: check cluster_results directory
        snp_matrix = None
        if os.path.isdir(args.cluster_results):
            candidate = os.path.join(args.cluster_results, cluster_name, 'snp_matrix.tsv')
            if os.path.exists(candidate):
                snp_matrix = candidate
        if snp_matrix:
            heatmap_html, ordered_samples, ordered_matrix = render_snp_heatmap_html(snp_matrix)

            # Save clustered SNP matrix TSV with annotations
            clustered_tsv = build_clustered_tsv(
                ordered_samples, ordered_matrix,
                species=species, mlst_map=mlst_map, amr_map=amr_map, metadata=meta_map
            )
            clustered_tsv_path = os.path.join(args.cluster_results, cluster_name, 'snp_matrix_clustered.tsv')
            os.makedirs(os.path.dirname(clustered_tsv_path), exist_ok=True)
            with open(clustered_tsv_path, 'w') as cf:
                cf.write(clustered_tsv)

            # Generate SVG heatmap for download
            heatmap_svg = build_heatmap_svg(ordered_samples, ordered_matrix)

            dl_tsv = make_download_button(clustered_tsv, f'{cluster_name}_snp_matrix.tsv', 'Download SNP matrix (TSV)')
            dl_svg = make_download_button(heatmap_svg, f'{cluster_name}_snp_heatmap.svg', 'Download heatmap (SVG)') if heatmap_svg else ''
            section_html += f'<div style="margin:12px 0 4px;"><h4 style="display:inline-block;margin:0;">SNP Distance Matrix</h4>{dl_tsv}{dl_svg}</div>'
            section_html += heatmap_html
        else:
            section_html += '<p style="color:#95a5a6;font-style:italic;font-size:12px;">SNP distance matrix not available for this cluster (alignment may have failed).</p>'

        # Minimum Spanning Tree (D3.js interactive)
        if snp_matrix:
            mst_html = render_mst_d3(snp_matrix, cluster_name, mlst_map=mlst_map, amr_map=amr_map, metadata=meta_map)
            if mst_html:
                section_html += '<div style="margin:12px 0 4px;"><h4 style="display:inline-block;margin:0;">Likely Transmission Network</h4></div>'
                section_html += mst_html

        # IQ-TREE SVG (pre-rendered): only for 4+ sample clusters
        tree_svg = os.path.join(args.tree_svgs, f'{cluster_name}_tree.svg')
        if os.path.exists(tree_svg):
            with open(tree_svg) as f:
                svg_content = f.read()
            # Skip empty or placeholder tree files (< 3 samples can't produce ML tree)
            if svg_content.strip() and '<svg' in svg_content:
                dl_btn = make_download_button(svg_content, f'{cluster_name}_tree.svg', 'Download tree SVG')
                tree_id = f"obtree_{cluster_name}".replace('-', '_').replace(' ', '_')
                section_html += (
                    f'<div style="margin:12px 0 4px;display:flex;align-items:center;gap:8px;">'
                    f'<h4 style="display:inline-block;margin:0;">Maximum Likelihood Tree</h4>{dl_btn}'
                    f'<span style="margin-left:auto;display:flex;gap:4px;">'
                    f'<button onclick="treeZoom(\'{tree_id}\', 1.3)" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;width:24px;height:24px;font-size:14px;line-height:1;color:#555;" title="Zoom in">+</button>'
                    f'<button onclick="treeZoom(\'{tree_id}\', 0.77)" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;width:24px;height:24px;font-size:14px;line-height:1;color:#555;" title="Zoom out">-</button>'
                    f'<button onclick="treeReset(\'{tree_id}\')" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;height:24px;font-size:10px;padding:0 6px;color:#555;" title="Reset zoom">Reset</button>'
                    f'</span></div>'
                    f'<div id="{tree_id}" style="overflow:auto;max-height:500px;border:1px solid #eee;border-radius:6px;background:#fafafa;padding:8px;cursor:grab;">'
                    f'<div id="{tree_id}_inner" style="transform-origin:top left;transition:transform 0.2s;">'
                    f'{svg_content}'
                    f'</div></div>'
                    f'<script>'
                    f'(function(){{'
                    f'var scales={{}};'
                    f'window.treeZoom=window.treeZoom||function(id,f){{scales[id]=(scales[id]||1)*f;document.getElementById(id+"_inner").style.transform="scale("+scales[id]+")";'
                    f'}};'
                    f'window.treeReset=window.treeReset||function(id){{scales[id]=1;document.getElementById(id+"_inner").style.transform="scale(1)";'
                    f'}};'
                    f'var el=document.getElementById("{tree_id}");'
                    f'var dragging=false,startX,startY,scrollL,scrollT;'
                    f'el.addEventListener("mousedown",function(e){{dragging=true;startX=e.pageX-el.offsetLeft;startY=e.pageY-el.offsetTop;scrollL=el.scrollLeft;scrollT=el.scrollTop;el.style.cursor="grabbing";}});'
                    f'el.addEventListener("mouseleave",function(){{dragging=false;el.style.cursor="grab";}});'
                    f'el.addEventListener("mouseup",function(){{dragging=false;el.style.cursor="grab";}});'
                    f'el.addEventListener("mousemove",function(e){{if(!dragging)return;e.preventDefault();el.scrollLeft=scrollL-(e.pageX-el.offsetLeft-startX);el.scrollTop=scrollT-(e.pageY-el.offsetTop-startY);}});'
                    f'}})();'
                    f'</script>'
                )

        # ── Core Genome Analysis (CGA) section ──
        if args.cga_cluster_results and os.path.isdir(args.cga_cluster_results):
            cga_snp_file = os.path.join(args.cga_cluster_results, cluster_name, 'cga_snp_matrix.tsv')
            if os.path.exists(cga_snp_file):
                section_html += '<hr style="margin:16px 0;border-color:#d5dbdb;">'

                cga_heatmap_html, cga_ordered_samples, cga_ordered_matrix = render_snp_heatmap_html(cga_snp_file)
                if cga_ordered_samples:
                    cga_clustered_tsv = build_clustered_tsv(
                        cga_ordered_samples, cga_ordered_matrix,
                        species=species, mlst_map=mlst_map, amr_map=amr_map, metadata=meta_map
                    )
                    cga_heatmap_svg = build_heatmap_svg(cga_ordered_samples, cga_ordered_matrix)
                    dl_tsv = make_download_button(cga_clustered_tsv, f'{cluster_name}_cga_snp_matrix.tsv', 'Download CGA SNP matrix (TSV)')
                    dl_svg = make_download_button(cga_heatmap_svg, f'{cluster_name}_cga_snp_heatmap.svg', 'Download CGA heatmap (SVG)') if cga_heatmap_svg else ''
                    section_html += f'<div style="margin:12px 0 4px;"><h4 style="display:inline-block;margin:0;">Core Genome SNP Distance Matrix</h4>{dl_tsv}{dl_svg}</div>'
                    section_html += cga_heatmap_html

                cga_mst_html = render_mst_d3(cga_snp_file, f'{cluster_name}_cga', mlst_map=mlst_map, amr_map=amr_map, metadata=meta_map)
                if cga_mst_html:
                    section_html += '<div style="margin:12px 0 4px;"><h4 style="display:inline-block;margin:0;">Core Genome Likely Transmission Network</h4></div>'
                    section_html += cga_mst_html

                if args.cga_tree_svgs and os.path.isdir(args.cga_tree_svgs):
                    cga_tree_svg = os.path.join(args.cga_tree_svgs, f'cga_{cluster_name}_tree.svg')
                    if os.path.exists(cga_tree_svg):
                        with open(cga_tree_svg) as f:
                            cga_svg_content = f.read()
                        if cga_svg_content.strip() and '<svg' in cga_svg_content:
                            dl_btn = make_download_button(cga_svg_content, f'{cluster_name}_cga_tree.svg', 'Download CGA tree SVG')
                            cga_tree_id = f"cgatree_{cluster_name}".replace('-', '_').replace(' ', '_')
                            section_html += (
                                f'<div style="margin:12px 0 4px;display:flex;align-items:center;gap:8px;">'
                                f'<h4 style="display:inline-block;margin:0;">Core Genome Maximum Likelihood Tree</h4>{dl_btn}'
                                f'<span style="margin-left:auto;display:flex;gap:4px;">'
                                f'<button onclick="treeZoom(\'{cga_tree_id}\', 1.3)" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;width:24px;height:24px;font-size:14px;line-height:1;color:#555;" title="Zoom in">+</button>'
                                f'<button onclick="treeZoom(\'{cga_tree_id}\', 0.77)" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;width:24px;height:24px;font-size:14px;line-height:1;color:#555;" title="Zoom out">-</button>'
                                f'<button onclick="treeReset(\'{cga_tree_id}\')" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;height:24px;font-size:10px;padding:0 6px;color:#555;" title="Reset zoom">Reset</button>'
                                f'</span></div>'
                                f'<div id="{cga_tree_id}" style="overflow:auto;max-height:500px;border:1px solid #eee;border-radius:6px;background:#fafafa;padding:8px;cursor:grab;">'
                                f'<div id="{cga_tree_id}_inner" style="transform-origin:top left;transition:transform 0.2s;">'
                                f'{cga_svg_content}'
                                f'</div></div>'
                                f'<script>'
                                f'(function(){{'
                                f'window.treeZoom=window.treeZoom||function(id,f){{var s=window._ts=window._ts||{{}};s[id]=(s[id]||1)*f;document.getElementById(id+"_inner").style.transform="scale("+s[id]+")";'
                                f'}};'
                                f'window.treeReset=window.treeReset||function(id){{var s=window._ts=window._ts||{{}};s[id]=1;document.getElementById(id+"_inner").style.transform="scale(1)";'
                                f'}};'
                                f'var el=document.getElementById("{cga_tree_id}");'
                                f'var dragging=false,startX,startY,scrollL,scrollT;'
                                f'el.addEventListener("mousedown",function(e){{dragging=true;startX=e.pageX-el.offsetLeft;startY=e.pageY-el.offsetTop;scrollL=el.scrollLeft;scrollT=el.scrollTop;el.style.cursor="grabbing";}});'
                                f'el.addEventListener("mouseleave",function(){{dragging=false;el.style.cursor="grab";}});'
                                f'el.addEventListener("mouseup",function(){{dragging=false;el.style.cursor="grab";}});'
                                f'el.addEventListener("mousemove",function(e){{if(!dragging)return;e.preventDefault();el.scrollLeft=scrollL-(e.pageX-el.offsetLeft-startX);el.scrollTop=scrollT-(e.pageY-el.offsetTop-startY);}});'
                                f'}})();'
                                f'</script>'
                            )

        section_html += '</div></details>'
        html_sections.append(section_html)
        n_cluster_sections += 1

    # ── Generate combined pairwise distances TSV ──
    import re as _re_asm
    pairwise_rows = []
    for cluster_name in sorted(clusters.keys()):
        cluster_info = clusters[cluster_name]
        species = cluster_info.get('species', 'Unknown')
        cluster_num = cluster_name.split("_")[-1]
        cluster_members = cluster_info.get('samples', [])
        n_in_cluster = len(cluster_members)

        # Read SNP matrix for this cluster
        snp_file = None
        if os.path.isdir(args.cluster_results):
            candidate = os.path.join(args.cluster_results, cluster_name, 'snp_matrix.tsv')
            if os.path.exists(candidate):
                snp_file = candidate
        if not snp_file:
            continue

        # Parse the SNP matrix
        snp_samples = []
        snp_matrix = []
        with open(snp_file) as sf:
            reader = csv.reader(sf, delimiter='\t')
            header = next(reader)
            for row in reader:
                name = os.path.basename(row[0]).replace('.fasta', '')
                snp_samples.append(name)
                snp_matrix.append([int(x) for x in row[1:]])

        # Generate all pairwise combinations
        for i in range(len(snp_samples)):
            for j in range(i + 1, len(snp_samples)):
                s1 = snp_samples[i]
                s2 = snp_samples[j]
                dist = snp_matrix[i][j]
                base1 = _re_asm.sub(r'_(L|S|HL|HS)$', '', s1)
                base2 = _re_asm.sub(r'_(L|S|HL|HS)$', '', s2)
                st1 = mlst_map.get(s1, '') or mlst_map.get(base1, '')
                st2 = mlst_map.get(s2, '') or mlst_map.get(base2, '')
                # Extract ASM type from sample name
                asm_match1 = _re_asm.search(r'_(L|S|HL|HS)$', s1)
                asm_type = asm_match1.group(1) if asm_match1 else ''

                pairwise_rows.append({
                    'workflow_id': args.workflow_id,
                    'workflow_timestamp': args.workflow_timestamp,
                    'group_id': args.group_id,
                    'species': species,
                    'cluster': cluster_num,
                    'mlst_1': st1 if st1 else '-',
                    'mlst_2': st2 if st2 else '-',
                    'isolate_1': s1,
                    'isolate_2': s2,
                    'pairwise_distance': dist,
                    'sample_count': n_in_cluster,
                    'asm_type': asm_type,
                })

    # Write combined pairwise TSV
    pairwise_tsv = os.path.join(os.path.dirname(args.output), 'ska2_pairwise_distances.tsv')
    with open(pairwise_tsv, 'w', newline='') as pf:
        w = csv.DictWriter(pf, delimiter='\t',
            fieldnames=['workflow_id', 'workflow_timestamp', 'group_id', 'species', 'cluster', 'mlst_1', 'mlst_2', 'isolate_1', 'isolate_2', 'pairwise_distance', 'sample_count', 'asm_type'])
        w.writeheader()
        for row in pairwise_rows:
            w.writerow(row)
    print(f"Combined pairwise distances: {len(pairwise_rows)} pairs written to {pairwise_tsv}")

    # ── Generate combined CGA pairwise distances TSV ──
    cga_pairwise_rows = []
    if args.cga_cluster_results and os.path.isdir(args.cga_cluster_results):
        for cluster_name in sorted(clusters.keys()):
            cluster_info = clusters[cluster_name]
            species = cluster_info.get('species', 'Unknown')
            cluster_num = cluster_name.split("_")[-1]
            cluster_members = cluster_info.get('samples', [])
            n_in_cluster = len(cluster_members)

            cga_snp_file = os.path.join(args.cga_cluster_results, cluster_name, 'cga_snp_matrix.tsv')
            if not os.path.exists(cga_snp_file):
                continue

            snp_samples = []
            snp_matrix = []
            with open(cga_snp_file) as sf:
                reader = csv.reader(sf, delimiter='\t')
                header = next(reader)
                for row in reader:
                    name = os.path.basename(row[0]).replace('.fasta', '')
                    snp_samples.append(name)
                    snp_matrix.append([int(x) for x in row[1:]])

            for i in range(len(snp_samples)):
                for j in range(i + 1, len(snp_samples)):
                    s1 = snp_samples[i]
                    s2 = snp_samples[j]
                    dist = snp_matrix[i][j]
                    base1 = _re_asm.sub(r'_(L|S|HL|HS)$', '', s1)
                    base2 = _re_asm.sub(r'_(L|S|HL|HS)$', '', s2)
                    st1 = mlst_map.get(s1, '') or mlst_map.get(base1, '')
                    st2 = mlst_map.get(s2, '') or mlst_map.get(base2, '')
                    asm_match1 = _re_asm.search(r'_(L|S|HL|HS)$', s1)
                    asm_type = asm_match1.group(1) if asm_match1 else ''

                    cga_pairwise_rows.append({
                        'workflow_id': args.workflow_id,
                        'workflow_timestamp': args.workflow_timestamp,
                        'group_id': args.group_id,
                        'species': species,
                        'cluster': cluster_num,
                        'mlst_1': st1 if st1 else '-',
                        'mlst_2': st2 if st2 else '-',
                        'isolate_1': s1,
                        'isolate_2': s2,
                        'pairwise_distance': dist,
                        'sample_count': n_in_cluster,
                        'asm_type': asm_type,
                    })

    cga_pairwise_tsv = os.path.join(os.path.dirname(args.output), 'cga_pairwise_distances.tsv')
    with open(cga_pairwise_tsv, 'w', newline='') as pf:
        w = csv.DictWriter(pf, delimiter='\t',
            fieldnames=['workflow_id', 'workflow_timestamp', 'group_id', 'species', 'cluster', 'mlst_1', 'mlst_2', 'isolate_1', 'isolate_2', 'pairwise_distance', 'sample_count', 'asm_type'])
        w.writeheader()
        for row in cga_pairwise_rows:
            w.writerow(row)
    print(f"Combined CGA pairwise distances: {len(cga_pairwise_rows)} pairs written to {cga_pairwise_tsv}")

    # Write as MultiQC custom content YAML for proper section naming.
    if n_cluster_sections == 0:
        open(args.output, 'w').close()
        print('No outbreak clusters: Outbreak Detection section suppressed')
    else:
        full_html = '\n'.join(html_sections)
        # Manual YAML output to avoid PyYAML dependency
        with open(args.output, 'w') as f:
            f.write("id: 'outbreak_detection'\n")
            f.write("section_name: 'Outbreak Detection'\n")
            f.write("description: '<em style=\"color:#c0392b;\">The data presented here should be interpreted with caution and supplemented by further epidemiological investigation before taking any decisions</em>'\n")
            f.write("plot_type: 'html'\n")
            f.write("data: |\n")
            for line in full_html.split('\n'):
                f.write(f"  {line}\n")

        print(f"Outbreak report written to {args.output}")


if __name__ == '__main__':
    main()