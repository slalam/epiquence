#!/usr/bin/env python3
"""Render one cumulative chromosomal distance matrix per species, covering every isolate."""

import argparse
import base64
import csv
import os
import struct
import sys
import zlib
from datetime import datetime

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
# Private palette constants on purpose: restating the hex values here would let this heatmap
# and the outbreak one drift apart on a future palette change.
from epiquence_lib.heatmap import (cluster_order, reorder, render_heatmap_svg,
                                   dual_scale_color, _GREEN, _YELLOW,
                                   _YELLOW_SNP_HI, _RED)
from epiquence_lib.surveillance_scope import make_filter, scope_label

ASM_SUFFIXES = ('_HL', '_HS', '_L', '_S', '_A')
LABEL_MIN_PX = 7          # below this a label is unreadable, so labels are dropped
TARGET_PX = 1600          # nominal on-screen size of the cell body
EMBED_TSV_MAX = 2 * 1024 * 1024


def now_in(tz_name):
    if tz_name:
        try:
            from zoneinfo import ZoneInfo
            return datetime.now(ZoneInfo(tz_name))
        except Exception as exc:
            print(f'WARN: timezone {tz_name!r} unavailable ({exc}); using system local time',
                  file=sys.stderr)
    return datetime.now()


def resolve_columns(headers):
    """Same column resolution as load_surveillance_db and surveillance_tree."""
    hm = {h.lower().replace(' ', '_').replace('(', '').replace(')', ''): h for h in headers}
    s1 = hm.get('sample1') or hm.get('sample_1') or hm.get('isolate_1')
    s2 = hm.get('sample2') or hm.get('sample_2') or hm.get('isolate_2')
    d = (hm.get('distance') or hm.get('snp_distance') or hm.get('snps')
         or hm.get('snpdistance'))
    return s1, s2, d


def stream_pairs(path, keep=None):
    """Yield (a, b, distance); ``keep`` optionally restricts both ends to a scope."""
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            return
        s1, s2, dc = resolve_columns(reader.fieldnames)
        if not (s1 and s2 and dc):
            sys.exit(f'ERROR: cannot resolve sample/distance columns in {path}; '
                     f'found {reader.fieldnames}')
        for row in reader:
            a, b = (row.get(s1) or '').strip(), (row.get(s2) or '').strip()
            raw = (row.get(dc) or '').strip()
            if not a or not b or a == b or not raw:
                continue
            if keep is not None and not (keep(a) and keep(b)):
                continue
            try:
                yield a, b, int(float(raw))
            except ValueError:
                continue



def assert_no_duplicate_pairs(path, limit=5, keep=None):
    """Fail loudly if the distances file names the same unordered pair twice."""
    seen = {}
    dupes = []
    for a, b, d in stream_pairs(path, keep):
        key = (a, b) if a <= b else (b, a)
        if key in seen and seen[key] != d:
            dupes.append((key[0], key[1], seen[key], d))
        seen[key] = d
    if dupes:
        sys.stderr.write(
            f"ERROR: {len(dupes)} pair(s) appear more than once with different "
            f"distances in {path}.\n"
            f"The merged sketch almost certainly holds a stale entry for a "
            f"re-sequenced isolate under the same sample name.\n"
            f"Rebuild it from the per-isolate sketches:\n"
            f"    rm <sketch_dir>/<Species>/merged/<Species>_merged.skf\n"
            f"and re-run; SKA_MERGE rebuilds from sketches/, which holds exactly "
            f"one sketch per sample.\n")
        for a, b, d1, d2 in dupes[:limit]:
            sys.stderr.write(f"    {a} vs {b}: {d1} and {d2}\n")
        if len(dupes) > limit:
            sys.stderr.write(f"    ... and {len(dupes) - limit} more\n")
        sys.exit(1)


def load_matrix(path, max_isolates, keep=None):
    """Two passes: index the samples, then fill an N x N int32 array."""
    index = {}
    for a, b, _ in stream_pairs(path, keep):
        if a not in index:
            index[a] = len(index)
        if b not in index:
            index[b] = len(index)
    n = len(index)
    if n > max_isolates:
        sys.exit(f'ERROR: {n} isolates exceeds --max-isolates ({max_isolates}). '
                 f'A {n}x{n} matrix is {n * n:,} cells; raise the limit only if you '
                 f'have the memory for it.')

    mat = np.zeros((n, n), dtype=np.int32)
    for a, b, d in stream_pairs(path, keep):
        i, j = index[a], index[b]
        mat[i, j] = d
        mat[j, i] = d

    names = [''] * n
    for name, i in index.items():
        names[i] = name
    return names, mat


def build_palette():
    """256-entry LUT for the dual-scale ramp."""
    pal = np.zeros((256, 3), dtype=np.uint8)
    for i in range(128):
        t = i / 127.0
        for c in range(3):
            pal[i, c] = int(_GREEN[c] + t * (_YELLOW[c] - _GREEN[c]))
    for i in range(128):
        t = i / 127.0
        for c in range(3):
            pal[128 + i, c] = int(_YELLOW_SNP_HI[c] + t * (_RED[c] - _YELLOW_SNP_HI[c]))
    return pal


PALETTE = build_palette()


def colour_indices(mat, break_at):
    """Vectorised dual-scale ramp -> (N, N) uint8 palette indices."""
    m = mat.astype(np.float64)
    max_val = float(m.max()) if m.size else 0.0
    brk = float(break_at) if break_at > 0 else (max_val or 1.0)

    lo = np.clip(m / brk, 0.0, 1.0)
    span = max(max_val - brk, 1.0)          # min_span=1
    hi = np.clip((m - brk) / span, 0.0, 1.0)

    idx = np.where(m <= brk,
                   np.rint(lo * 127.0),
                   128.0 + np.rint(hi * 127.0))
    return idx.astype(np.uint8), max_val


def png_data_uri(idx, palette):
    """Minimal palette PNG encoder (colour type 3), so no Pillow is needed."""
    h, w = idx.shape
    raw = np.empty((h, w + 1), dtype=np.uint8)
    raw[:, 0] = 0                                    # filter type 0 per scanline
    raw[:, 1:] = idx

    def chunk(tag, data):
        return (struct.pack('>I', len(data)) + tag + data
                + struct.pack('>I', zlib.crc32(tag + data) & 0xffffffff))

    png = (b'\x89PNG\r\n\x1a\n'
           + chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 3, 0, 0, 0))
           + chunk(b'PLTE', palette.tobytes())
           + chunk(b'IDAT', zlib.compress(raw.tobytes(), 9))
           + chunk(b'IEND', b''))
    return 'data:image/png;base64,' + base64.b64encode(png).decode('ascii')


def _esc(t):
    return (str(t).replace('&', '&amp;').replace('<', '&lt;')
            .replace('>', '&gt;').replace('"', '&quot;'))


def snp_color_at(break_at):
    """snp_color with a configurable break point."""
    def fn(val, max_val):
        return dual_scale_color(float(val), float(max_val), float(break_at),
                                yellow_hi=_YELLOW_SNP_HI, min_span=1)
    return fn


def build_vector_svg(labels, mat, break_at):
    """The outbreak cluster heatmap, unchanged: one <rect> per cell with the distance printed
    in it."""
    svg = render_heatmap_svg(labels, mat, color_fn=snp_color_at(break_at),
                             cell_size=40)
    # id so the download button can serialise whichever variant is on the page
    return svg.replace('<svg xmlns=', '<svg id="epiq-heatmap" xmlns=', 1)


def build_svg(labels, mat, break_at, max_val):
    """SVG with a raster cell body and vector labels/legend."""
    n = len(labels)
    cell = max(1, min(24, TARGET_PX // n if n else 24))
    body = n * cell
    show_labels = cell >= LABEL_MIN_PX

    idx, _ = colour_indices(mat, break_at)
    if cell > 1:
        idx = np.repeat(np.repeat(idx, cell, axis=0), cell, axis=1)
    uri = png_data_uri(idx, PALETTE)

    if show_labels:
        gutter = min(240, max(120, int(max(len(l) for l in labels) * 6.6)))
        top = left = gutter
    else:
        top = left = 12
    legend_h = 58
    w = left + body + 30
    h = top + body + legend_h + 20

    parts = [
        f'<svg id="epiq-heatmap" xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{w}" height="{h}" viewBox="0 0 {w} {h}" font-family="Inter,Arial,sans-serif">',
        f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
        f'<image x="{left}" y="{top}" width="{body}" height="{body}" '
        f'image-rendering="pixelated" xlink:href="{uri}"/>',
        f'<rect x="{left}" y="{top}" width="{body}" height="{body}" fill="none" stroke="#ccc"/>',
    ]

    if show_labels:
        fs = min(11, cell - 1)
        for i, lab in enumerate(labels):
            cy = top + i * cell + cell / 2 + 3
            parts.append(f'<text x="{left - 6}" y="{cy:.1f}" font-size="{fs}" '
                         f'text-anchor="end" fill="#333">{_esc(lab)}</text>')
            cx = left + i * cell + cell / 2
            parts.append(f'<text x="{cx:.1f}" y="{top - 6}" font-size="{fs}" '
                         f'text-anchor="start" fill="#333" '
                         f'transform="rotate(-60 {cx:.1f} {top - 6})">{_esc(lab)}</text>')

    # Legend uses the same ramp, so the break point is visible rather than implied.
    ly = top + body + 26
    grad_w, grad_h = 260, 12
    stops = []
    for k in range(21):
        t = k / 20.0
        val = t * max_val
        if val <= break_at:
            u = val / break_at if break_at else 0.0
            col = tuple(int(_GREEN[c] + u * (_YELLOW[c] - _GREEN[c])) for c in range(3))
        else:
            u = min(1.0, (val - break_at) / max(max_val - break_at, 1.0))
            col = tuple(int(_YELLOW_SNP_HI[c] + u * (_RED[c] - _YELLOW_SNP_HI[c])) for c in range(3))
        stops.append(f'<stop offset="{t * 100:.0f}%" stop-color="rgb{col}"/>')
    parts.append('<defs><linearGradient id="snpramp" x1="0" x2="1">'
                 + ''.join(stops) + '</linearGradient></defs>')
    parts.append(f'<rect x="{left}" y="{ly}" width="{grad_w}" height="{grad_h}" '
                 f'fill="url(#snpramp)" stroke="#ccc"/>')
    parts.append(f'<text x="{left}" y="{ly + grad_h + 13}" font-size="10" fill="#666">0</text>')
    parts.append(f'<text x="{left + grad_w}" y="{ly + grad_h + 13}" font-size="10" '
                 f'text-anchor="end" fill="#666">{int(max_val)}</text>')
    if 0 < break_at < max_val:
        bx = left + grad_w * (break_at / max_val)
        parts.append(f'<line x1="{bx:.1f}" y1="{ly - 3}" x2="{bx:.1f}" y2="{ly + grad_h + 3}" '
                     f'stroke="#333" stroke-width="1.5"/>')
        parts.append(f'<text x="{bx:.1f}" y="{ly - 6}" font-size="10" text-anchor="middle" '
                     f'fill="#333">{int(break_at)}</text>')

    parts.append('</svg>')
    return '\n'.join(parts), cell, show_labels


def base_isolate_id(sample):
    for sfx in ASM_SUFFIXES:
        if sample.endswith(sfx):
            return sample[:-len(sfx)]
    return sample


def write_square_tsv(path, labels, mat, species, isolates):
    """Square matrix with annotation columns, in the clustered leaf order."""
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Isolate_ID', 'Species', 'MLST', 'Chromosomal_AMR_Genes',
                    'First_Seen', 'Sample_ID'] + labels)
        for i, s in enumerate(labels):
            iso = isolates.get(s, {})
            w.writerow([base_isolate_id(s), species, iso.get('mlst_st', ''),
                        iso.get('chromosomal_amr_genes', ''),
                        iso.get('first_seen_workflow', ''), s]
                       + [int(v) for v in mat[i]])


def load_isolates(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline='', encoding='utf-8') as f:
        return {(r.get('sample') or '').strip(): r
                for r in csv.DictReader(f, delimiter='\t')
                if (r.get('sample') or '').strip()}


def main():
    p = argparse.ArgumentParser(
        description='Render the cumulative chromosomal distance matrix for one species')
    p.add_argument('--distances', required=True)
    p.add_argument('--isolates', required=True)
    p.add_argument('--species', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--matrix-out', required=True,
                   help='Sidecar square-matrix TSV covering every isolate')
    p.add_argument('--color-break', type=int, default=100,
                   help='Distance at which the ramp turns from green/yellow to '
                        'yellow/red; 100 matches the outbreak cluster heatmap')
    p.add_argument('--vector-max', type=int, default=300,
                   help='Below this isolate count the outbreak-style vector matrix '
                        'is drawn, with the distance printed in every cell. At or '
                        'above it the cells are rasterised, since one rect plus one '
                        'text node per cell stops being viable.')
    p.add_argument('--max-isolates', type=int, default=4000)
    p.add_argument('--timezone', default=None)
    p.add_argument('--group', default=None,
                   help='Restrict to isolates whose cumulative group_id equals this '
                        '(an empty string selects isolates without a group). '
                        'Omit for the species-wide view.')
    args = p.parse_args()

    if not os.path.exists(args.distances):
        print(f'[{args.species}] no distances file; skipping', file=sys.stderr)
        return

    isolates = load_isolates(args.isolates)
    keep = make_filter(isolates, args.group)
    label = scope_label(isolates, args.group)
    tag = args.species + (f' / group {args.group or "(none)"}' if args.group is not None else '')

    assert_no_duplicate_pairs(args.distances, keep=keep)
    names, mat = load_matrix(args.distances, args.max_isolates, keep)
    n = len(names)
    if n < 2:
        print(f'[{tag}] {n} isolate(s) in scope; skipping', file=sys.stderr)
        return

    order = cluster_order(mat)
    labels, omat = reorder(names, mat, order)
    omat = np.rint(omat).astype(np.int32)
    max_val = float(omat.max())
    vector = n < args.vector_max

    write_square_tsv(args.matrix_out, labels, omat, args.species, isolates)
    tsv_size = os.path.getsize(args.matrix_out)

    if vector:
        svg = build_vector_svg(labels, omat, args.color_break)
        cell, labelled = 40, True
    else:
        svg, cell, labelled = build_svg(labels, omat, args.color_break, max_val)

    safe = args.species.replace(' ', '_')

    # Serialise the inline SVG in the browser rather than base64-embedding a second copy.
    _btn_css = ("display:inline-flex;align-items:center;gap:5px;font-size:11px;"
                "line-height:1.4;box-sizing:border-box;color:#4e79a7;"
                "text-decoration:none;border:1px solid #ccc;border-radius:4px;"
                "padding:5px 10px;background:#fff;cursor:pointer;"
                "font-family:inherit;font-weight:400;")
    dl_svg = (f'<button type="button" id="epiq-dl-svg" style="{_btn_css}">'
              f'Download heatmap (SVG)</button>')

    matrix_name = os.path.basename(args.matrix_out)
    if tsv_size <= EMBED_TSV_MAX:
        with open(args.matrix_out, encoding='utf-8') as f:
            _b64 = base64.b64encode(f.read().encode('utf-8')).decode('ascii')
        dl_tsv = (f'<a href="data:text/tab-separated-values;base64,{_b64}" '
                  f'download="{matrix_name}" style="{_btn_css}">'
                  f'Download distance matrix (TSV)</a>')
    else:
        # Too large to base64 into the page; it sits next to this file instead.
        dl_tsv = (f'<a href="{matrix_name}" download="{matrix_name}" style="{_btn_css}">'
                  f'Download distance matrix (TSV, {tsv_size / 1048576:.0f} MB)</a>')

    note = '' if vector else (
        f'<div style="margin-top:8px;font-size:11px;color:#a8590c;background:#fdf6ec;'
        f'border-left:3px solid #e8a33d;padding:6px 10px;">Axis labels are omitted at '
        f'{n} isolates; at or above {args.vector_max} the cells are rasterised '
        f'and each is {cell}px, so per-cell values and labels would be illegible. '
        f'The TSV is in the same clustered order as the heatmap, so row <i>i</i> there '
        f'is row <i>i</i> here.</div>')

    generated = now_in(args.timezone).strftime('%Y-%m-%d %H:%M %Z').strip()
    filename = f'{safe}_distance_matrix.svg'
    doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_esc(args.species)}{(' | ' + label) if label else ''} | chromosomal distance matrix</title>
<style>
  body {{ margin:0; padding:16px; background:#f4f4f4;
         font-family:Inter,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; color:#333; }}
  .epiq-wrap {{ max-width:1800px; margin:0 auto; }}
  .epiq-card {{ background:#fff; border:1px solid #e0e0e0; border-radius:4px;
                padding:14px 18px; margin-bottom:10px; }}
  .epiq-title {{ font-size:17px; font-weight:600; color:#2c3e50; font-style:italic; }}
  .epiq-meta {{ margin-top:6px; font-size:11px; color:#888; display:flex;
                flex-wrap:wrap; gap:16px; }}
  .epiq-meta b {{ color:#4e79a7; font-weight:600; }}
  .epiq-plot {{ overflow:auto; }}
  .epiq-scope {{ margin-top:2px; font-size:12px; color:#3d6b5e; font-weight:600; }}
</style>
</head>
<body>
<div class="epiq-wrap">
  <div class="epiq-card">
    <div class="epiq-title">{_esc(args.species)}</div>
    {('<div class="epiq-scope">' + label + '</div>') if label else ''}
    <div class="epiq-meta">
      <span><b>{n}</b> isolates</span>
      <span>max distance <b>{int(max_val)}</b></span>
      <span>generated <b>{generated}</b></span>
    </div>
    <div style="display:flex;gap:6px;margin-top:10px;">{dl_svg}{dl_tsv}</div>
    {note}
  </div>
  <div class="epiq-card epiq-plot">
{svg}
  </div>
</div>
<script>
(function () {{
  var btn = document.getElementById('epiq-dl-svg');
  var svgEl = document.getElementById('epiq-heatmap');
  if (!btn || !svgEl) return;

  btn.addEventListener('click', function () {{
    try {{
      var src = new XMLSerializer().serializeToString(svgEl);
      if (src.indexOf('xmlns=') === -1) {{
        src = src.replace('<svg', '<svg xmlns="http://www.w3.org/2000/svg"');
      }}
      src = '<?xml version="1.0" encoding="UTF-8"?>' + String.fromCharCode(10) + src;

      var a = document.createElement('a');
      a.download = {filename!r};
      a.style.display = 'none';
      // Firefox ignores click() on an anchor that is not in the document, which
      // is why a detached anchor appears to do nothing at all.
      document.body.appendChild(a);

      if (window.URL && window.URL.createObjectURL) {{
        var url = URL.createObjectURL(new Blob([src], {{type: 'image/svg+xml;charset=utf-8'}}));
        a.href = url;
        a.click();
        setTimeout(function () {{
          URL.revokeObjectURL(url);
          if (a.parentNode) a.parentNode.removeChild(a);
        }}, 1000);
      }} else {{
        // Very old browsers: data URI instead of a blob.
        a.href = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(src);
        a.click();
        if (a.parentNode) a.parentNode.removeChild(a);
      }}
    }} catch (err) {{
      // Silent failure is what made this hard to diagnose the first time.
      alert('SVG download failed: ' + (err && err.message ? err.message : err));
    }}
  }});
}})();
</script>
</body>
</html>
"""
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(doc)

    print(f'[{tag}] distance matrix: {n} isolates, '
          f'{"vector" if vector else "raster"}, {cell}px cells, '
          f'html {os.path.getsize(args.output) / 1048576:.1f} MB, '
          f'tsv {tsv_size / 1048576:.1f} MB', file=sys.stderr)


if __name__ == '__main__':
    main()