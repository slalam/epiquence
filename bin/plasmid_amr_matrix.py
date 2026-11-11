#!/usr/bin/env python3
"""AMR-carrying plasmid Mash distance matrix with hierarchical clustering."""

import argparse
import csv
import os
import sys

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib import heatmap as hm
from epiquence_lib.mash import read_matrix
from epiquence_lib.names import plasmid_display_name, split_plasmid_name
from epiquence_lib.reporting import make_download_button


def amr_plasmid_names(amr_dir):
    """FASTA stems in the SELECT_AMR_PLASMIDS output directory."""
    if not amr_dir or not os.path.isdir(amr_dir):
        return []
    return sorted(
        f[:-len('.fasta')] for f in os.listdir(amr_dir) if f.endswith('.fasta')
    )


def load_metadata(path):
    """plasmid_hgt_metadata.tsv -> {plasmid_name: {size_bp, resistance_genes, species}}."""
    meta = {}
    if not path or not os.path.exists(path) or os.path.basename(path).startswith('NO_'):
        return meta
    try:
        with open(path) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                name = (row.get('plasmid_name') or '').strip()
                if not name:
                    continue
                try:
                    size = int(row.get('size_bp', 0) or 0)
                except (ValueError, TypeError):
                    size = 0
                meta[name] = {
                    'size_bp': size,
                    'resistance_genes': (row.get('resistance_genes') or '').strip(),
                    'species': (row.get('species') or '').strip(),
                }
    except OSError as exc:
        print(f"  WARN: cannot read metadata {path}: {exc}", file=sys.stderr)
    return meta


def subset_matrix(matrix_path, keep):
    """Subset a Mash matrix to ``keep``, preserving matrix row order."""
    row_names, col_names, rows = read_matrix(matrix_path, basename=True, strip_fasta=True)
    keep_set = set(keep)

    col_index = {name: j for j, name in enumerate(col_names)}
    names, out = [], []
    for i, rname in enumerate(row_names):
        if rname not in keep_set:
            continue
        names.append(rname)

    for rname in names:
        i = row_names.index(rname)
        row = []
        for cname in names:
            j = col_index.get(cname)
            row.append(rows[i][j] if j is not None and j < len(rows[i]) else 1.0)
        out.append(row)
    return names, out


def annotate(name, meta, html=True):
    """Left-hand row label: display name + the AMR genes the plasmid carries."""
    info = meta.get(name, {})
    disp = plasmid_display_name(name)
    genes = [g.strip() for g in info.get('resistance_genes', '').replace(';', ',').split(',') if g.strip()]
    if not genes:
        return disp
    joined = '; '.join(genes)
    if not html:
        return f'{disp} ({joined})'
    return f'{disp} <span style="color:#888;font-weight:400;">({joined})</span>'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mash-matrix', required=True, help='Full plasmid mash_dist_matrix.tsv')
    ap.add_argument('--amr-dir', required=True, help='SELECT_AMR_PLASMIDS output directory')
    ap.add_argument('--asm-type', required=True)
    ap.add_argument('--metadata', default=None, help='plasmid_hgt_metadata.tsv')
    ap.add_argument('--threshold', type=float, default=0.005,
                    help='Mash clustering threshold; sets the colour-ramp break point')
    ap.add_argument('--max-plasmids', type=int, default=120,
                    help='Skip rendering the table above this many plasmids (TSV is still written)')
    ap.add_argument('--output-prefix', required=True)
    args = ap.parse_args()

    prefix = args.output_prefix
    html_path = f'{prefix}_amr_mash_matrix.html'
    tsv_path = f'{prefix}_amr_mash_matrix.tsv'
    svg_path = f'{prefix}_amr_mash_matrix.svg'
    order_path = f'{prefix}_amr_mash_order.txt'

    def write_placeholder(msg):
        with open(html_path, 'w') as f:
            f.write(f'<p style="color:#95a5a6;font-style:italic;">{msg}</p>')
        open(tsv_path, 'a').close()
        open(svg_path, 'a').close()
        open(order_path, 'a').close()
        print(msg, file=sys.stderr)

    amr_names = amr_plasmid_names(args.amr_dir)
    print(f"AMR-carrying plasmids in {args.amr_dir}: {len(amr_names)}", file=sys.stderr)

    if len(amr_names) < 2:
        write_placeholder(
            f'Fewer than two AMR-carrying plasmids in {args.asm_type} '
            f'({len(amr_names)}); no distance matrix to show.'
        )
        return

    if not os.path.exists(args.mash_matrix) or os.path.getsize(args.mash_matrix) < 20:
        write_placeholder('No Mash distance matrix available for this assembly type.')
        return

    names, rows = subset_matrix(args.mash_matrix, amr_names)
    missing = sorted(set(amr_names) - set(names))
    if missing:
        print(f"  WARN: {len(missing)} AMR plasmid(s) absent from the Mash matrix: "
              f"{', '.join(missing[:5])}{'...' if len(missing) > 5 else ''}", file=sys.stderr)

    if len(names) < 2:
        write_placeholder('Fewer than two AMR-carrying plasmids present in the Mash matrix.')
        return

    meta = load_metadata(args.metadata)

    # --- Hierarchical (UPGMA) leaf order ---
    order = hm.cluster_order(rows, method='average')
    ordered_names, ordered_matrix = hm.reorder(names, rows, order)
    n = len(ordered_names)

    with open(order_path, 'w') as f:
        for nm in ordered_names:
            f.write(nm + '\n')

    # --- Clustered TSV with annotation columns ---
    with open(tsv_path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['plasmid', 'sample', 'mob_cluster', 'size_bp', 'species',
                    'resistance_genes'] + [plasmid_display_name(x) for x in ordered_names])
        for i, nm in enumerate(ordered_names):
            info = meta.get(nm, {})
            sample, asm, mob = split_plasmid_name(nm)
            w.writerow([
                plasmid_display_name(nm), sample, mob,
                info.get('size_bp', ''), info.get('species', ''),
                info.get('resistance_genes', ''),
            ] + [f'{ordered_matrix[i][j]:.5f}' for j in range(n)])

    # --- Report fragment ---
    max_val = max(
        (ordered_matrix[i][j] for i in range(n) for j in range(n) if i != j),
        default=0.0,
    )

    with open(tsv_path) as f:
        tsv_content = f.read()
    dl_tsv = make_download_button(tsv_content, f'{args.asm_type}_amr_mash_matrix.tsv',
                                  'Download clustered matrix TSV')

    # Vector version of the same heatmap.
    svg = hm.render_heatmap_svg(
        ordered_names, ordered_matrix,
        color_fn=lambda v, mx: hm.mash_color(v, mx, args.threshold),
        fmt=lambda v: f'{v:.4f}'.lstrip('0') if v else '0',
        cell_size=40,
        col_labels=[plasmid_display_name(nm) for nm in ordered_names],
        row_labels=[annotate(nm, meta, html=False) for nm in ordered_names],
    )
    with open(svg_path, 'w') as f:
        f.write(svg)
    dl_svg = make_download_button(svg, f'{args.asm_type}_amr_mash_matrix.svg',
                                  'Download matrix SVG')

    if n > args.max_plasmids:
        body = (
            f'<p style="color:#95a5a6;font-style:italic;">'
            f'{n} plasmids exceeds the inline rendering limit of {args.max_plasmids}; '
            f'the clustered matrix is available as a TSV download, and an SVG of the '
            f'full heatmap is in the results directory.</p>'
            f'<div style="margin:8px 0;">{dl_tsv}</div>'
        )
    else:
        def tooltip(row_label, col_label, val):
            a = plasmid_display_name(row_label)
            b = plasmid_display_name(col_label)
            return f'{a} vs {b}: {val:.5f}'

        table = hm.render_heatmap_html(
            ordered_names, ordered_matrix,
            color_fn=lambda v, mx: hm.mash_color(v, mx, args.threshold),
            fmt=lambda v: f'{v:.4f}'.lstrip('0') if v else '0',
            cell_size=40,
            col_labels=[plasmid_display_name(nm) for nm in ordered_names],
            row_labels=[annotate(nm, meta) for nm in ordered_names],
            tooltips=tooltip,
            max_label_chars=None,          # show column names in full
        )
        body = (f'<div style="margin:8px 0;display:flex;align-items:center;gap:6px;">'
                f'{dl_tsv}{dl_svg}</div>{table}')

    # The arrow rotation CSS lives in the report page, but this fragment is also published
    # standalone, so carry a scoped copy of the marker rules with it.
    marker_css = (
        '<style>'
        '.eq-amr-mx > summary::-webkit-details-marker{display:none;}'
        '.eq-amr-mx > summary{list-style:none;}'
        '.eq-amr-mx > summary span:first-child{display:inline-block;transition:transform .2s;}'
        '.eq-amr-mx[open] > summary span:first-child{transform:rotate(90deg);}'
        '</style>'
    )
    section = (
        f'{marker_css}'
        f'<details class="eq-amr-mx" style="margin:12px 0;border:1px solid #e0e0e0;border-radius:8px;">'
        f'<summary style="cursor:pointer;padding:10px 14px;font-weight:600;background:#f7f7f7;'
        f'border-radius:8px 8px 0 0;list-style:none;">'
        f'<span style="margin-right:8px;">&#9654;</span>'
        f'AMR Plasmid Mash Distance Matrix ({args.asm_type})</summary>'
        f'<div style="padding:12px 16px;">{body}</div>'
        f'</details>'
    )

    with open(html_path, 'w') as f:
        f.write(section)

    print(f"AMR plasmid Mash matrix: {args.asm_type} | {n} plasmids, "
          f"max distance {max_val:.5f}, break at {args.threshold:g}", file=sys.stderr)


if __name__ == '__main__':
    main()