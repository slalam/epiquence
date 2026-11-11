#!/usr/bin/env python3
"""Generate Plasmid HGT Detection report for MultiQC."""
import argparse, csv, json, os, glob, re
import base64
try:
    import yaml
except ImportError:
    yaml = None


def make_download_button(content, filename, icon_title="Download"):
    """Create a small download icon button using base64 data URI."""
    b64 = base64.b64encode(content.encode('utf-8')).decode('utf-8')
    mime = 'image/svg+xml' if filename.endswith('.svg') else 'text/tab-separated-values'
    icon_svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" '
        'fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">'
        '<path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/>'
        '<polyline points="7 10 12 15 17 10"/>'
        '<line x1="12" y1="15" x2="12" y2="3"/>'
        '</svg>'
    )
    return (
        f'<a href="data:{mime};base64,{b64}" download="{filename}" '
        f'title="{icon_title}" '
        f'style="display:inline-flex;align-items:center;justify-content:center;'
        f'width:26px;height:26px;border-radius:4px;background:#f0f0f0;'
        f'color:#555;text-decoration:none;margin-left:6px;vertical-align:middle;'
        f'border:1px solid #ddd;" '
        f'onmouseover="this.style.background=\'#e0e0e0\'" '
        f'onmouseout="this.style.background=\'#f0f0f0\'">'
        f'{icon_svg}</a>'
    )


def display_name(raw):
    """SAMPLE-03_HL_plasmid_AA102 -> SAMPLE-03_HL_pAA102"""
    return re.sub(r'_plasmid_([A-Z]+)', r'_p\1', raw)


def split_name(raw):
    """Split plasmid name into (sample_id, asm_type, mob_cluster)."""
    # Find _plasmid_ separator first
    m = re.match(r'^(.+)_plasmid_(.+)$', raw)
    if m:
        prefix = m.group(1)  # e.g., SAMPLE_A-06_HL_A or SAMPLE-03_HL
        mob = m.group(2)     # e.g., AA276
        # Extract asm_type from the end of prefix
        asm_m = re.match(r'^(.+)_(L|S|HL|HS|A)$', prefix)
        if asm_m:
            sample = asm_m.group(1)  # e.g., SAMPLE_A-06_HL or SAMPLE-03
            asm = asm_m.group(2)     # e.g., A or HL
            return f"{sample}_{asm}", asm, mob
        return prefix, '', mob
    return raw, '', ''


def load_metadata(filepath):
    """Load plasmid metadata TSV: plasmid_name -> {size_bp, resistance_genes, species}"""
    meta = {}
    if not filepath or not os.path.exists(filepath):
        return meta
    with open(filepath) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            name = row['plasmid_name']
            meta[name] = {
                'size_bp': int(row.get('size_bp', 0) or 0),
                'resistance_genes': row.get('resistance_genes', ''),
                'species': row.get('species', ''),
            }
    return meta


def build_cluster_label(cluster_idx, ani_data, asm_type='', metadata=None):
    """Build enriched cluster label with styled HTML: Bold: cluster ID | MOB IDs Grey: sample
    IDs Red: distinct AMR genes (if present)"""
    mob_ids = set()
    sample_ids = set()
    seen = set()
    plasmid_names = set()
    for d in ani_data:
        pair = tuple(sorted([d['query'], d['reference']]))
        if pair in seen:
            continue
        seen.add(pair)
        for name in (d['query'], d['reference']):
            s, a, m = split_name(name)
            if m:
                mob_ids.add(m)
            if s:
                sample_ids.add(s)
            plasmid_names.add(name)
    mob_str = '; '.join(f'p{m}' for m in sorted(mob_ids))
    sample_str = '; '.join(sorted(sample_ids))
    prefix = f"{asm_type}_C{cluster_idx}" if asm_type else f"C{cluster_idx}"
    n_samples = len(sample_ids)
    sample_label = f' <span style="font-size:11px;color:#999;font-weight:normal;">{sample_str}</span>' if sample_str else ''

    # Collect distinct AMR genes across all plasmids in this cluster
    amr_label = ''
    if metadata:
        all_genes = set()
        for pname in plasmid_names:
            genes_str = metadata.get(pname, {}).get('resistance_genes', '')
            if genes_str and genes_str != '-':
                for g in genes_str.replace(';', ',').split(','):
                    g = g.strip()
                    if g:
                        all_genes.add(g)
        if all_genes:
            genes_sorted = ', '.join(sorted(all_genes))
            amr_label = (
                f'<br><span style="font-size:10px;color:#c0392b;font-weight:500;'
                f'letter-spacing:0.2px;margin-left:24px;">'
                f'{genes_sorted}</span>'
            )

    return f"{prefix} | {mob_str} ({n_samples} samples){sample_label}{amr_label}"


def parse_ani_af(filepath):
    data = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ani_val = float(row.get('ANI', 0))
            snps_val = int(row.get('TotalSNPs', 0))
            af_q = float(row.get('AF_query', 0))

            # Fix dnadiff quirk: identical sequences report ANI=0.00
            if ani_val == 0.0 and snps_val == 0 and af_q > 0:
                ani_val = 100.0

            data.append({
                'query': row['query'],
                'reference': row['reference'],
                'ANI': ani_val,
                'AF_query': float(row.get('AF_query', 0)),
                'AF_ref': float(row.get('AF_ref', 0)),
                'TotalSNPs': snps_val,
                'TotalIndels': int(row.get('TotalIndels', 0)),
                'AlignedBases': row.get('AlignedBases', '0'),
                'query_size': int(row.get('query_size', 0) or 0),
                'ref_size': int(row.get('ref_size', 0) or 0),
            })
    return data


def build_table_yaml(all_rows, asm_type, asm_label, metadata):
    """Build native MultiQC table YAML with conditional formatting."""

    # Collect unique values for conditional formatting
    cluster_vals = sorted(set(r['cluster_label'] for r in all_rows))
    q_mob_vals = sorted(set(r['query_mob'] for r in all_rows if r['query_mob']))
    r_mob_vals = sorted(set(r['ref_mob'] for r in all_rows if r['ref_mob']))
    all_mob_vals = sorted(set(q_mob_vals) | set(r_mob_vals))

    # Color palette (same as Plasmid Resistance Markers)
    palette = [
        '#1565c0', '#c62828', '#2e7d32', '#6a1b9a', '#e65100',
        '#00838f', '#ad1457', '#4527a0', '#00695c', '#bf360c',
        '#283593', '#558b2f', '#d84315', '#0277bd', '#6d4c41',
        '#4e342e', '#37474f', '#827717', '#880e4f', '#1b5e20',
    ]

    # Build cluster conditional formatting
    cluster_rules = {}
    cluster_colours = []
    for i, cv in enumerate(cluster_vals):
        rule = f'hgt_cl_{cv}'
        cluster_rules[rule] = [{'s_eq': cv}]
        cluster_colours.append({rule: palette[i % len(palette)]})

    # Build mob ID conditional formatting
    mob_rules = {}
    mob_colours = []
    for i, mv in enumerate(all_mob_vals):
        rule = f'hgt_mob_{mv}'
        mob_rules[rule] = [{'s_eq': mv}]
        mob_colours.append({rule: palette[i % len(palette)]})

    headers = {
        'Cluster': {
            'title': 'Cluster', 'description': 'HGT cluster ID', 'scale': False,
            'cond_formatting_rules': dict(cluster_rules),
            'cond_formatting_colours': list(cluster_colours),
        },
        'Query_Sample':  {'title': 'Query Sample', 'description': 'Query sample ID', 'scale': False},
        'Query_Mob_ID': {
            'title': 'Query Mob ID', 'description': 'Query MOB-recon plasmid ID', 'scale': False,
            'cond_formatting_rules': dict(mob_rules),
            'cond_formatting_colours': list(mob_colours),
        },
        'Query_Size_KB': {'title': 'Query Size (KB)', 'description': 'Query plasmid size in kilobases', 'min': 0, 'scale': 'PuBu', 'format': '{:.1f}', 'suffix': ' KB'},
        'Ref_Sample':    {'title': 'Ref Sample', 'description': 'Reference sample ID', 'scale': False},
        'Ref_Mob_ID': {
            'title': 'Ref Mob ID', 'description': 'Reference MOB-recon plasmid ID', 'scale': False,
            'cond_formatting_rules': dict(mob_rules),
            'cond_formatting_colours': list(mob_colours),
        },
        'Ref_Size_KB': {'title': 'Ref Size (KB)', 'description': 'Reference plasmid size in kilobases', 'min': 0, 'scale': 'PuBu', 'format': '{:.1f}', 'suffix': ' KB'},
        'ANI':        {'title': 'ANI (%)', 'description': 'Average Nucleotide Identity', 'min': 99, 'max': 100, 'scale': 'RdYlGn', 'format': '{:.2f}', 'suffix': '%'},
        'AF_QtoR':    {'title': 'AF Q\u2192R (%)', 'description': 'Alignment fraction query to reference', 'min': 0, 'max': 100, 'scale': 'Blues', 'format': '{:.2f}', 'suffix': '%'},
        'AF_RtoQ':    {'title': 'AF R\u2192Q (%)', 'description': 'Alignment fraction reference to query', 'min': 0, 'max': 100, 'scale': 'Blues', 'format': '{:.2f}', 'suffix': '%'},
        'SNPs':       {'title': 'SNPs', 'description': 'Total SNPs between pair', 'scale': 'OrRd', 'format': '{:,d}'},
        'Indels':     {'title': 'Indels', 'description': 'Total indels between pair', 'scale': 'OrRd', 'format': '{:,d}'},
        'Query_Res_Genes': {'title': 'Query Res Genes', 'description': 'Resistance genes on query plasmid', 'scale': False},
        'Ref_Res_Genes': {'title': 'Ref Res Genes', 'description': 'Resistance genes on reference plasmid', 'scale': False},
        'ASM_Type':   {'title': 'ASM Type', 'description': 'Assembly type', 'scale': False},
    }

    # Strip empty formatting if no values
    if not cluster_rules:
        del headers['Cluster']['cond_formatting_rules']
        del headers['Cluster']['cond_formatting_colours']
    if not mob_rules:
        del headers['Query_Mob_ID']['cond_formatting_rules']
        del headers['Query_Mob_ID']['cond_formatting_colours']
        del headers['Ref_Mob_ID']['cond_formatting_rules']
        del headers['Ref_Mob_ID']['cond_formatting_colours']

    data = {}
    for i, r in enumerate(all_rows):
        q_disp = display_name(r['query_raw'])
        r_disp = display_name(r['ref_raw'])
        row_key = f"{q_disp} vs {r_disp}"

        # Look up metadata for sizes and resistance genes
        q_meta = metadata.get(r['query_raw'], {})
        r_meta = metadata.get(r['ref_raw'], {})
        # Size: prefer metadata (from FASTA), fall back to dnadiff TotalBases
        q_size_bp = q_meta.get('size_bp', 0) or r.get('query_size', 0)
        r_size_bp = r_meta.get('size_bp', 0) or r.get('ref_size', 0)
        q_size_kb = round(q_size_bp / 1000, 1) if q_size_bp else 0.0
        r_size_kb = round(r_size_bp / 1000, 1) if r_size_bp else 0.0
        q_genes = q_meta.get('resistance_genes', '')
        r_genes = r_meta.get('resistance_genes', '')
        # Convert semicolons to commas for display
        q_genes_disp = q_genes.replace(';', ', ') if q_genes else '-'
        r_genes_disp = r_genes.replace(';', ', ') if r_genes else '-'

        data[row_key] = {
            'Cluster': r['cluster_label'],
            'ASM_Type': r['asm_type'],
            'Query_Sample': r['query_sample'],
            'Query_Mob_ID': r['query_mob'],
            'Query_Size_KB': q_size_kb,
            'Query_Res_Genes': q_genes_disp,
            'Ref_Sample': r['ref_sample'],
            'Ref_Mob_ID': r['ref_mob'],
            'Ref_Size_KB': r_size_kb,
            'Ref_Res_Genes': r_genes_disp,
            'ANI': r['ANI'],
            'AF_QtoR': r['AF_query'],
            'AF_RtoQ': r['AF_ref'],
            'SNPs': r['TotalSNPs'],
            'Indels': r['TotalIndels'],
        }

    return {
        'id': f'plasmid_hgt_1_table_{asm_type}',
        'parent_id': 'plasmid_hgt_detection',
        'parent_name': 'Plasmid HGT Detection',
        'parent_description': 'Plasmid-level horizontal gene transfer detection based on Mash distance clustering and phylogenetic analysis',
        'section_name': f'{asm_label} HGT | Pairwise Plasmid Comparison',
        'description': '',
        'plot_type': 'table',
        'pconfig': {
            'id': f'plasmid_hgt_1_table_{asm_type}',
            'title': f'{asm_label} HGT | Pairwise Plasmid Comparison',
            'sortRows': True,
            'col1_header': 'Comparison',
            'no_violin': True,
        },
        'headers': headers,
        'data': data,
    }


def generate_scatter(ani_data, div_id):
    """ANI vs AF scatter using Altair/Vega-Lite: hover shows ANI, AF, SNPs, Indels."""
    import altair as alt
    import pandas as pd
    import numpy as np

    if not ani_data:
        return ''

    seen = set()
    rows = []
    for d in ani_data:
        pair = tuple(sorted([d['query'], d['reference']]))
        if pair in seen:
            continue
        seen.add(pair)
        q_disp = display_name(d['query'])
        r_disp = display_name(d['reference'])
        rows.append({
            'AF_query': d['AF_query'],
            'ANI': d['ANI'],
            'AF_ref': d['AF_ref'],
            'SNPs': d['TotalSNPs'],
            'Indels': d['TotalIndels'],
            'Comparison': f"{q_disp} vs {r_disp}",
        })

    if not rows:
        return ''

    df = pd.DataFrame(rows)

    # Add jitter to separate overlapping points.
    np.random.seed(42)
    af_range = max(1, df['AF_query'].max() - df['AF_query'].min())
    ani_range = max(0.001, df['ANI'].max() - df['ANI'].min())
    jitter_af = af_range * 0.012   # ~1.2% of AF range
    jitter_ani = ani_range * 0.015  # ~1.5% of ANI range
    df['AF_plot'] = df['AF_query'] + np.random.uniform(-jitter_af, jitter_af, len(df))
    df['ANI_plot'] = df['ANI'] + np.random.uniform(-jitter_ani, jitter_ani, len(df))

    max_snps = max(1, df['SNPs'].max())

    chart = alt.Chart(df).mark_circle(
        size=60,
        stroke='#2c3e50',
        strokeWidth=0.5,
    ).encode(
        x=alt.X('AF_plot:Q',
                 title='Alignment Fraction Q\u2192R (%)',
                 scale=alt.Scale(domain=[
                     max(0, df['AF_query'].min() - 2),
                     min(102, df['AF_query'].max() + 2)
                 ])),
        y=alt.Y('ANI_plot:Q',
                 title='ANI (%)',
                 scale=alt.Scale(domain=[
                     max(90.0, df['ANI'].min() - 0.1),
                     min(100.05, df['ANI'].max() + 0.05)
                 ])),
        color=alt.Color('SNPs:Q',
                        title='SNPs',
                        scale=alt.Scale(
                            domain=[0, max_snps],
                            scheme='redyellowgreen',
                            reverse=True
                        ),
                        legend=alt.Legend(orient='right')),
        tooltip=[
            alt.Tooltip('Comparison:N', title='Pair'),
            alt.Tooltip('ANI:Q', title='ANI (%)', format='.2f'),
            alt.Tooltip('AF_query:Q', title='AF Q\u2192R (%)', format='.2f'),
            alt.Tooltip('AF_ref:Q', title='AF R\u2192Q (%)', format='.2f'),
            alt.Tooltip('SNPs:Q', title='SNPs'),
            alt.Tooltip('Indels:Q', title='Indels'),
        ]
    ).properties(
        width='container',
        height=220,
    ).configure_view(
        strokeWidth=0
    ).configure(
        font='Inter, -apple-system, Helvetica Neue, sans-serif',
        background='#fafafa',
    )

    spec_json = chart.to_json(indent=None)

    return f"""
    <div id="{div_id}" style="width:100%;"></div>
    <script>
    (function() {{
        var spec = {spec_json};
        vegaEmbed('#{div_id}', spec, {{actions: false, renderer: 'svg'}});
    }})();
    </script>
    """


def generate_amr_tree_section(asm_type, tree_svgs_dir):
    """Single Mash tree across all AMR-carrying plasmids for this assembly type."""
    if not tree_svgs_dir:
        return ''
    svg_path = os.path.join(tree_svgs_dir, f'{asm_type}_amr_plasmids_tree.svg')
    if not os.path.exists(svg_path):
        return ''
    with open(svg_path) as f:
        svg_content = f.read()

    dl_btn = make_download_button(svg_content, f'{asm_type}_amr_plasmids_tree.svg',
                                  'Download tree SVG')
    tree_id = f"tree_{asm_type}_amr_plasmids".replace('-', '_')

    return (
        f'<details style="margin:12px 0;border:1px solid #e0e0e0;border-radius:8px;">'
        f'<summary style="cursor:pointer;padding:10px 14px;font-weight:600;background:#f7f7f7;'
        f'border-radius:8px 8px 0 0;list-style:none;">'
        f'<span style="margin-right:8px;">&#9654;</span>AMR Plasmid Tree ({asm_type})</summary>'
        f'<div style="padding:12px 16px;">'
        f'<div style="margin:8px 0 4px;display:flex;align-items:center;gap:8px;">'
        f'{dl_btn}'
        f'<span style="margin-left:auto;display:flex;gap:4px;">'
        f'<button onclick="treeZoom(\'{tree_id}\', 1.3)" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;width:24px;height:24px;font-size:14px;line-height:1;color:#555;" title="Zoom in">+</button>'
        f'<button onclick="treeZoom(\'{tree_id}\', 0.77)" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;width:24px;height:24px;font-size:14px;line-height:1;color:#555;" title="Zoom out">-</button>'
        f'<button onclick="treeReset(\'{tree_id}\')" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;height:24px;font-size:10px;padding:0 6px;color:#555;" title="Reset zoom">Reset</button>'
        f'</span></div>'
        f'<div id="{tree_id}" style="overflow:auto;max-height:600px;border:1px solid #eee;'
        f'border-radius:6px;background:#fafafa;padding:8px;cursor:grab;">'
        f'<div id="{tree_id}_inner" style="transform-origin:top left;transition:transform 0.2s;">'
        f'{svg_content}'
        f'</div></div>'
        f'<script>'
        f'(function(){{'
        f'if(!window.treeZoom){{var scales={{}};'
        f'window.treeZoom=function(id,fac){{scales[id]=(scales[id]||1)*fac;'
        f'document.getElementById(id+"_inner").style.transform="scale("+scales[id]+")";}};'
        f'window.treeReset=function(id){{scales[id]=1;'
        f'document.getElementById(id+"_inner").style.transform="scale(1)";}};}}'
        f'}})();'
        f'</script>'
        f'</div></details>'
    )


def generate_amr_matrix_section(html_path):
    """AMR plasmid Mash distance matrix, pre-rendered by plasmid_amr_matrix.py."""
    if not html_path or not os.path.exists(html_path):
        return ''
    if os.path.getsize(html_path) < 50:
        return ''
    with open(html_path) as f:
        return f.read()


def generate_cluster_section(ani_data, asm_type, cluster_id, cluster_idx, tree_svgs_dir, comparison_maps_dir=None, metadata=None):
    """Generate full HTML section for a cluster: header + tree SVG + scatter + comparison map."""
    if not ani_data:
        return ''

    cluster_label = build_cluster_label(cluster_idx, ani_data, asm_type, metadata=metadata)
    div_id = f"scatter_{asm_type}_{cluster_id}"

    samples = set()
    for d in ani_data:
        s1, _, _ = split_name(d['query'])
        s2, _, _ = split_name(d['reference'])
        if s1: samples.add(s1)
        if s2: samples.add(s2)
    n_samples = len(samples)

    # Body is assembled first: a cluster whose only content would be the header (no tree SVG, and
    # no comparison map: either because --hgt_maps is off or because the cluster carries no AMR
    # and so got no map) is dropped entirely rather than rendered as an empty collapsible.
    parts = []

    # Load pre-rendered tree SVG from toytree
    if tree_svgs_dir:
        tree_svg_path = os.path.join(tree_svgs_dir, f'{cluster_id}_tree.svg')
        if os.path.exists(tree_svg_path):
            with open(tree_svg_path) as f:
                svg_content = f.read()
            dl_btn = make_download_button(svg_content, f'{cluster_id}_tree.svg', 'Download tree SVG')

            # Unique ID for this tree container
            tree_id = f"tree_{asm_type}_{cluster_id}".replace('-', '_')

            parts.append(f'<div style="margin:8px 0 4px;display:flex;align-items:center;gap:8px;">'
                         f'<span style="font-size:12px;font-weight:600;color:#555;">Mash Cluster Tree</span>'
                         f'{dl_btn}'
                         f'<span style="margin-left:auto;display:flex;gap:4px;">'
                         f'<button onclick="treeZoom(\'{tree_id}\', 1.3)" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;width:24px;height:24px;font-size:14px;line-height:1;color:#555;" title="Zoom in">+</button>'
                         f'<button onclick="treeZoom(\'{tree_id}\', 0.77)" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;width:24px;height:24px;font-size:14px;line-height:1;color:#555;" title="Zoom out">-</button>'
                         f'<button onclick="treeReset(\'{tree_id}\')" style="cursor:pointer;border:1px solid #ddd;background:#f8f8f8;border-radius:4px;height:24px;font-size:10px;padding:0 6px;color:#555;" title="Reset zoom">Reset</button>'
                         f'</span>'
                         f'</div>')

            # Scrollable + zoomable container
            parts.append(
                f'<div id="{tree_id}" style="overflow:auto;max-height:500px;border:1px solid #eee;'
                f'border-radius:6px;background:#fafafa;padding:8px;cursor:grab;">'
                f'<div id="{tree_id}_inner" style="transform-origin:top left;transition:transform 0.2s;">'
                f'{svg_content}'
                f'</div></div>'
                f'<script>'
                f'(function(){{'
                f'var scales={{}};'
                f'window.treeZoom=function(id,f){{scales[id]=(scales[id]||1)*f;'
                f'document.getElementById(id+"_inner").style.transform="scale("+scales[id]+")";'
                f'}};'
                f'window.treeReset=function(id){{scales[id]=1;'
                f'document.getElementById(id+"_inner").style.transform="scale(1)";'
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

    # Comparison map - embed as iframe via srcdoc ── Comparison maps (supports multiple groups per
    # cluster) ──
    if comparison_maps_dir:
        import glob as _glob
        # Find maps for this cluster only: exact match with optional group suffix Pattern:
        # {cluster_id}_map.html OR {cluster_id}_{a-z}_map.html
        all_htmls = sorted(_glob.glob(os.path.join(comparison_maps_dir, f'{cluster_id}_*.html')))
        map_files = []
        for mf in all_htmls:
            fname = os.path.basename(mf)
            # Match exact cluster: cluster_id_map.html or cluster_id_X_map.html (single letter group)
            if fname == f'{cluster_id}_map.html':
                map_files.append(mf)
            elif re.match(rf'^{re.escape(cluster_id)}_[a-z]_map\.html$', fname):
                map_files.append(mf)

        if len(map_files) == 1:
            # Single map: show directly with title
            with open(map_files[0]) as f:
                map_content = f.read()
            # Extract title for display
            title_match = re.search(r'<title>([^<]+)</title>', map_content)
            map_title = title_match.group(1) if title_match else ''
            detail_parts = []
            if 'circular' in map_title.lower():
                detail_parts.append('circular')
            amr_match = re.search(r'(\d+)\s+AMR\s+gene', map_title)
            if amr_match:
                detail_parts.append(f"{amr_match.group(1)} AMR genes")
            len_match = re.search(r'([\d,]+)\s*bp', map_title)
            if len_match:
                detail_parts.append(f"{len_match.group(1)} bp")
            detail_str = ' / '.join(detail_parts)
            header_html = f'<div style="font-size:12px;color:#666;margin-bottom:6px;">{detail_str}</div>' if detail_str else ''

            escaped = map_content.replace('&', '&amp;').replace('"', '&quot;')
            parts.append(
                f'<div style="margin:16px 0;">'
                f'{header_html}'
                f'<iframe srcdoc="{escaped}" '
                f'style="width:100%;height:800px;border:none;border-radius:8px;'
                f'box-shadow:0 2px 8px rgba(0,0,0,0.06);" '
                f'sandbox="allow-scripts allow-downloads"></iframe>'
                f'</div>'
            )
        elif len(map_files) > 1:
            # Multiple maps - tab selector
            tab_cid = f"maptabs_{asm_type}_{cluster_id}".replace('-', '_')
            tab_btns = []
            tab_pnls = []
            for idx, mf in enumerate(map_files):
                fname = os.path.basename(mf)
                suffix_match = re.search(rf'{re.escape(cluster_id)}_([a-z])_map\.html', fname)
                glbl = suffix_match.group(1).upper() if suffix_match else str(idx + 1)
                with open(mf) as f:
                    mc = f.read()
                tmatch = re.search(r'<title>([^<]+)</title>', mc)
                mtitle = tmatch.group(1) if tmatch else ''
                dp = []
                if 'circular' in mtitle.lower():
                    dp.append('circular')
                else:
                    dp.append('linear')
                am = re.search(r'(\d+)\s+AMR\s+gene', mtitle)
                if am:
                    dp.append(f"{am.group(1)} AMR")
                lm = re.search(r'([\d,]+)\s*bp', mtitle)
                if lm:
                    dp.append(f"{lm.group(1)} bp")
                ds = ' / '.join(dp)
                abg = '#fff' if idx == 0 else '#f5f5f5'
                afw = '600' if idx == 0 else '400'
                acl = '#333' if idx == 0 else '#888'
                abb = '2px solid #4a90d9' if idx == 0 else '1px solid #ddd'
                dsp = 'block' if idx == 0 else 'none'
                tid = f"{tab_cid}_{idx}"
                # Extract AMR gene names from the HTML content
                amr_names_match = re.search(r'var amrGenes\s*=\s*"([^"]*)"', mc)
                amr_names = amr_names_match.group(1) if amr_names_match and amr_names_match.group(1) else ''

                tab_btns.append(
                    f'<button class="mt" data-t="{tid}" data-c="{tab_cid}" onclick="smt(this)" '
                    f'style="padding:6px 14px;border:1px solid #ddd;background:{abg};border-radius:6px 6px 0 0;'
                    f'font-size:12px;font-weight:{afw};color:{acl};cursor:pointer;border-bottom:{abb};margin-right:4px;text-align:left;">'
                    f'Map {glbl} <span style="font-size:10px;color:#999;margin-left:6px;">{ds}</span>'
                    f'{"<br><span style=font-size:10px;color:#c0392b;font-weight:400;>" + amr_names + "</span>" if amr_names else ""}'
                    f'</button>')
                esc = mc.replace('&', '&amp;').replace('"', '&quot;')
                tab_pnls.append(
                    f'<div id="{tid}" class="mp" data-c="{tab_cid}" style="display:{dsp};">'
                    f'<iframe srcdoc="{esc}" style="width:100%;height:800px;border:none;border-radius:0 0 8px 8px;'
                    f'box-shadow:0 2px 8px rgba(0,0,0,0.06);" sandbox="allow-scripts allow-downloads"></iframe></div>')
            parts.append(
                f'<div style="margin:16px 0;">'
                f'<div style="display:flex;flex-wrap:wrap;border-bottom:1px solid #ddd;">{"".join(tab_btns)}</div>'
                f'{"".join(tab_pnls)}</div>'
                f'<script>window.smt=window.smt||function(b){{'
                f'var c=b.getAttribute("data-c"),t=b.getAttribute("data-t");'
                f'document.querySelectorAll(".mp[data-c=\'"+c+"\']").forEach(function(p){{p.style.display="none";}});'
                f'document.querySelectorAll(".mt[data-c=\'"+c+"\']").forEach(function(x){{'
                f'x.style.background="#f5f5f5";x.style.fontWeight="400";x.style.color="#888";x.style.borderBottom="1px solid #ddd";}});'
                f'document.getElementById(t).style.display="block";'
                f'b.style.background="#fff";b.style.fontWeight="600";b.style.color="#333";b.style.borderBottom="2px solid #4a90d9";}};'
                f'</script>'
            )

    if not parts:
        return ''

    header = (
        f'<details style="margin:8px 0;border:1px solid #dee2e6;border-radius:6px;overflow:hidden;">'
        f'<summary style="padding:10px 14px;background:#f8f9fa;cursor:pointer;font-size:14px;font-weight:600;list-style:none;">'
        f'<span style="margin-right:8px;">&#9654;</span>'
        f'{cluster_label}'
        f'</summary>'
        f'<div style="padding:12px 16px;">'
    )
    return '\n'.join([header] + parts + ['</div></details>'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--asm-type', required=True)
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--metadata', default=None, help='Plasmid metadata TSV (sizes + resistance genes)')
    parser.add_argument('--tree-svgs', default=None, help='Directory with pre-rendered tree SVGs')
    parser.add_argument('--comparison-maps', default=None, help='Directory with per-cluster CGView comparison map HTMLs')
    parser.add_argument('--umap-html', default=None, help='UMAP embedding HTML file')
    parser.add_argument('--amr-matrix-html', default=None,
                        help='AMR plasmid Mash distance matrix HTML fragment')
    parser.add_argument('--table-yaml', required=True, help='Output path for MultiQC table YAML')
    parser.add_argument('--plots-yaml', required=True, help='Output path for MultiQC plots YAML')
    parser.add_argument('--combined-tsv', default=None)
    parser.add_argument('--workflow-id', default='')
    parser.add_argument('--group-id', default='')
    args = parser.parse_args()

    asm_type = args.asm_type
    asm_label = {
        'S': 'Short Read (S)', 'L': 'Long Read (L)',
        'HL': 'Hybrid Long First (HL)', 'HS': 'Hybrid Short First (HS)'
    }.get(asm_type, asm_type)

    metadata = load_metadata(args.metadata)

    import re as _re_sort
    def _natural_sort_key(path):
        name = os.path.basename(path)
        parts = _re_sort.split(r'(\d+)', name)
        return [int(p) if p.isdigit() else p.lower() for p in parts]

    cluster_dirs = sorted(glob.glob(os.path.join(args.results_dir, 'hgt_cluster_*')), key=_natural_sort_key)

    all_table_rows = []
    scatter_sections = []
    tsv_rows = []

    for idx, cdir in enumerate(cluster_dirs, 1):
        cluster_id = os.path.basename(cdir)
        # Extract cluster number from hgt_cluster_N format
        import re as _re
        _cnum_match = _re.search(r'(\d+)$', cluster_id)
        cnum = _cnum_match.group(1) if _cnum_match else str(idx)
        cluster_label = f"{asm_type}_C{cnum}"
        ani_file = os.path.join(cdir, f'{cluster_id}_ani_af.tsv')

        if not os.path.exists(ani_file):
            continue

        ani_data = parse_ani_af(ani_file)
        if not ani_data:
            continue

        seen = set()
        for d in ani_data:
            pair = tuple(sorted([d['query'], d['reference']]))
            if pair in seen:
                continue
            seen.add(pair)

            q_sample, q_asm, q_mob = split_name(d['query'])
            r_sample, r_asm, r_mob = split_name(d['reference'])

            row = {
                'cluster_label': cluster_label,
                'cluster_id': cluster_id,
                'asm_type': asm_type,
                'query_raw': d['query'],
                'ref_raw': d['reference'],
                'query_sample': q_sample,
                'query_mob': q_mob,
                'ref_sample': r_sample,
                'ref_mob': r_mob,
                'ANI': d['ANI'],
                'AF_query': d['AF_query'],
                'AF_ref': d['AF_ref'],
                'TotalSNPs': d['TotalSNPs'],
                'TotalIndels': d['TotalIndels'],
                'AlignedBases': d['AlignedBases'],
                'query_size': d['query_size'],
                'ref_size': d['ref_size'],
            }
            all_table_rows.append(row)
            tsv_rows.append(row)

        scatter_sections.append(
            generate_cluster_section(ani_data, asm_type, cluster_id, cnum, None, args.comparison_maps, metadata=metadata)
        )

    # Subheading over the per-cluster comparison maps.
    if any(sec.strip() for sec in scatter_sections):
        scatter_sections.insert(0,
            f'<div style="margin:20px 0 4px;">'
            f'<h4 style="margin:0;font-size:15px;font-weight:600;color:#333;">'
            f'AMR Plasmid Comparison Maps ({asm_type})</h4></div>'
        )

    # Prepend the assembly-type-wide AMR plasmid views (replace per-cluster trees).
    amr_matrix_html = generate_amr_matrix_section(args.amr_matrix_html)
    if amr_matrix_html:
        scatter_sections.insert(0, amr_matrix_html)

    amr_tree_html = generate_amr_tree_section(asm_type, args.tree_svgs)
    if amr_tree_html:
        scatter_sections.insert(0, amr_tree_html)

    table_data = build_table_yaml(all_table_rows, asm_type, asm_label, metadata)
    if yaml:
        with open(args.table_yaml, 'w') as f:
            yaml.dump(table_data, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
    else:
        # Fallback: write JSON (MultiQC also accepts JSON with _mqc suffix)
        with open(args.table_yaml, 'w') as f:
            json.dump(table_data, f, indent=2, ensure_ascii=False)

    html_parts = [
        '<style>'
        'details summary::-webkit-details-marker { display: none; }'
        'details summary { list-style: none; }'
        'details[open] > summary span:first-child { transform: rotate(90deg); display: inline-block; }'
        'details summary span:first-child { transition: transform 0.2s; display: inline-block; }'
        'details[open] { border-color: #adb5bd !important; }'
        '</style>',
        '<script src="https://cdn.jsdelivr.net/npm/vega@5"></script>',
        '<script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>',
        '<script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>',
    ]

    # UMAP embedding plot at the top
    if args.umap_html and os.path.exists(args.umap_html) and os.path.getsize(args.umap_html) > 100:
        with open(args.umap_html) as f:
            umap_content = f.read()
        escaped_umap = umap_content.replace('&', '&amp;').replace('"', '&quot;')
        html_parts.append(
            '<details style="margin:12px 0;border:1px solid #e0e0e0;border-radius:8px;">'
            '<summary style="cursor:pointer;padding:10px 14px;font-weight:600;background:#f7f7f7;'
            'border-radius:8px 8px 0 0;list-style:none;">'
            '<span style="margin-right:8px;">&#9654;</span>Plasmid UMAP Embedding</summary>'
            '<div style="padding:12px 16px;">'
            f'<iframe srcdoc="{escaped_umap}" '
            'style="width:100%;height:550px;border:1px solid #e0e0e0;border-radius:4px;" '
            'sandbox="allow-scripts allow-downloads"></iframe>'
            '</div></details>'
        )

    # Build cluster summary as native MultiQC table YAML
    hgt_summary_table = {}
    if all_table_rows:
        cluster_summary = {}
        for r in all_table_rows:
            cl = r['cluster_label']
            if cl not in cluster_summary:
                cluster_summary[cl] = {'mob_ids': set(), 'sample_ids': set(), 'raw_names': set(), 'res_genes': set(), 'species': set()}
            for name_key in ['query_raw', 'ref_raw']:
                s, a, m = split_name(r[name_key])
                if m:
                    cluster_summary[cl]['mob_ids'].add(m)
                if s:
                    cluster_summary[cl]['sample_ids'].add(s)
                cluster_summary[cl]['raw_names'].add(r[name_key])
                if metadata:
                    meta = metadata.get(r[name_key], {})
                    genes_str = meta.get('resistance_genes', '')
                    if genes_str and genes_str != '-':
                        for g in genes_str.replace(';', ',').split(','):
                            g = g.strip()
                            if g:
                                cluster_summary[cl]['res_genes'].add(g)
                    sp = meta.get('species', '')
                    if sp and sp != '-':
                        cluster_summary[cl]['species'].add(sp)

        def cluster_sort_key(x):
            parts = x.rsplit('_C', 1)
            prefix = parts[0]
            num = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0
            return (prefix, num)

        for cl in sorted(cluster_summary.keys(), key=cluster_sort_key):
            info = cluster_summary[cl]
            species_str = '; '.join(sorted(info['species'])) if info['species'] else '-'
            mob_str = '; '.join(f'p{m}' for m in sorted(info['mob_ids']))
            sample_str = '; '.join(sorted(info['sample_ids']))
            plasmid_count = len(info['raw_names'])
            genes_str = '; '.join(sorted(info['res_genes'])) if info['res_genes'] else '-'

            hgt_summary_table[cl] = {
                'MOB_IDs': mob_str,
                'Plasmids': plasmid_count,
                'Species': species_str,
                'Samples': sample_str,
                'Resistance_Genes': genes_str,
            }

    # Write HGT cluster summary as native MultiQC table
    summary_yaml_path = os.path.join(os.path.dirname(args.plots_yaml), f'plasmid_hgt_{asm_type}_summary_mqc.yaml')
    import json as _json_yaml
    
    yaml_data = {
        'id': f'plasmid_hgt_summary_{asm_type}',
        'parent_id': 'plasmid_hgt_detection',
        'parent_name': 'Plasmid HGT Detection',
        'parent_description': 'Plasmid-level horizontal gene transfer detection based on Mash distance clustering and phylogenetic analysis',
        'section_name': f'{asm_label} HGT | Cluster Summary',
        'description': '',
        'plot_type': 'table',
        'pconfig': {
            'id': f'hgt_cluster_summary_table_{asm_type}',
            'title': f'{asm_label} Plasmid HGT Cluster Summary',
            'sortRows': True,
            'col1_header': 'Cluster ID',
        },
        'headers': {
            'MOB_IDs': {'title': 'MOB IDs', 'scale': False},
            'Plasmids': {'title': 'Plasmids', 'format': '{:,.0f}', 'scale': 'Blues'},
            'Species': {'title': 'Species', 'scale': False},
            'Samples': {'title': 'Samples', 'scale': False},
            'Resistance_Genes': {'title': 'Resistance Genes', 'scale': False},
        },
        'data': hgt_summary_table,
    }
    
    with open(summary_yaml_path, 'w') as f:
        _json_yaml.dump(yaml_data, f, indent=2, default=str)
    print(f"HGT cluster summary table written to {summary_yaml_path}")

    if scatter_sections:
        for s in scatter_sections:
            if s:
                html_parts.append(s)
    else:
        html_parts.append('<p style="color:#95a5a6;font-style:italic;">No plasmid HGT clusters detected.</p>')

    full_html = '\n'.join(html_parts)

    with open(args.plots_yaml, 'w') as f:
        f.write(f"id: 'plasmid_hgt_2_plots_{asm_type}'\n")
        f.write("parent_id: 'plasmid_hgt_detection'\n")
        f.write("parent_name: 'Plasmid HGT Detection'\n")
        f.write("parent_description: 'Plasmid-level horizontal gene transfer detection based on Mash distance clustering and phylogenetic analysis'\n")
        f.write(f"section_name: '{asm_label} HGT | Cluster Detection and Plots'\n")
        f.write("description: ''\n")
        f.write("plot_type: 'html'\n")
        f.write("data: |\n")
        for line in full_html.split('\n'):
            f.write(f"  {line}\n")

    if args.combined_tsv and tsv_rows:
        with open(args.combined_tsv, 'w', newline='') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['workflow_id', 'group_id', 'assembly_type', 'cluster',
                        'query_sample', 'query_mob_id', 'query_size_bp', 'query_res_genes',
                        'ref_sample', 'ref_mob_id', 'ref_size_bp', 'ref_res_genes',
                        'ANI', 'AF_query', 'AF_ref', 'TotalSNPs', 'TotalIndels', 'AlignedBases'])
            for r in tsv_rows:
                q_meta = metadata.get(r['query_raw'], {})
                r_meta = metadata.get(r['ref_raw'], {})
                q_size = q_meta.get('size_bp', 0) or r.get('query_size', 0)
                r_size = r_meta.get('size_bp', 0) or r.get('ref_size', 0)
                q_genes = q_meta.get('resistance_genes', '')
                r_genes = r_meta.get('resistance_genes', '')
                w.writerow([args.workflow_id, args.group_id, r['asm_type'], r['cluster_label'],
                            r['query_sample'], r['query_mob'], q_size, q_genes,
                            r['ref_sample'], r['ref_mob'], r_size, r_genes,
                            r['ANI'], r['AF_query'], r['AF_ref'],
                            r['TotalSNPs'], r['TotalIndels'], r['AlignedBases']])

    print(f"HGT report: {asm_type} | {len(cluster_dirs)} clusters, {len(all_table_rows)} pairs")


if __name__ == '__main__':
    main()