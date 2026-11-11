"""Distance-matrix heatmap rendering with hierarchical (UPGMA) leaf ordering."""

import math

import numpy as np

try:
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform
    _HAVE_SCIPY = True
except ImportError:  # pragma: no cover - scipy is present in the report containers
    _HAVE_SCIPY = False



def symmetrize(matrix):
    """Return a float ndarray with a zero diagonal and ``(M + M.T) / 2``."""
    mat = np.array(matrix, dtype=float)
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError(f"expected a square matrix, got shape {mat.shape}")
    np.fill_diagonal(mat, 0)
    return (mat + mat.T) / 2


def cluster_order(matrix, method='average'):
    """UPGMA leaf order for a square distance matrix."""
    n = len(matrix)
    if n < 3 or not _HAVE_SCIPY:
        return list(range(n))
    sym = symmetrize(matrix)
    condensed = squareform(sym, checks=False)
    return list(leaves_list(linkage(condensed, method=method)))


def reorder(labels, matrix, order):
    """Apply a leaf order to labels and matrix together."""
    mat = np.array(matrix, dtype=float)
    return [labels[i] for i in order], mat[np.ix_(order, order)]



def _ramp(t, lo, hi):
    """Linear interpolation between two RGB triples."""
    t = min(1.0, max(0.0, t))
    return tuple(int(lo[i] + t * (hi[i] - lo[i])) for i in range(3))


_GREEN = (26, 122, 46)          # #1a7a2e
_YELLOW = (230, 199, 0)         # #e6c700: top of the lower ramp
_YELLOW_SNP_HI = (230, 180, 0)  # #e6b400: bottom of the outbreak upper ramp
_RED = (180, 30, 30)            # #b41e1e


def dual_scale_color(val, max_val, break_at, yellow_hi=None, min_span=None):
    """``(css_background, css_text_color)`` for one cell."""
    if break_at <= 0:
        break_at = max_val if max_val > 0 else 1.0
    if yellow_hi is None:
        yellow_hi = _YELLOW

    if val <= break_at:
        t = val / break_at if break_at else 0.0
        r, g, b = _ramp(t, _GREEN, _YELLOW)
        txt = '#fff' if t < 0.5 else '#333'
    else:
        span = max_val - break_at
        span = max(span, min_span) if min_span is not None else max(span, 1e-12)
        t = min(1.0, (val - break_at) / span)
        r, g, b = _ramp(t, yellow_hi, _RED)
        txt = '#fff' if t > 0.3 else '#333'
    return f'rgb({r},{g},{b})', txt


def snp_color(val, max_val):
    """Outbreak SNP ramp: break at 100 SNPs, byte-identical to the original."""
    return dual_scale_color(float(val), float(max_val), 100.0,
                            yellow_hi=_YELLOW_SNP_HI, min_span=1)


def mash_color(val, max_val, threshold=0.005):
    """Plasmid Mash ramp: break at the HGT clustering threshold, continuous."""
    return dual_scale_color(float(val), float(max_val), float(threshold))



def _esc(text):
    """Escape a string for use inside a double-quoted HTML attribute."""
    return (str(text).replace('&', '&amp;').replace('<', '&lt;')
            .replace('>', '&gt;').replace('"', '&quot;'))


def render_heatmap_html(labels, matrix, color_fn, fmt=None, cell_size=40,
                        row_labels=None, col_labels=None, tooltips=None,
                        max_label_chars=14):
    """Render a clustered distance matrix as a self-contained HTML table."""
    n = len(labels)
    if n == 0:
        return '<p style="color:#95a5a6;font-style:italic;">Empty distance matrix</p>'

    mat = np.array(matrix, dtype=float)
    off_diagonal = mat[~np.eye(n, dtype=bool)] if n > 1 else np.array([0.0])
    max_val = float(off_diagonal.max()) if off_diagonal.size else 0.0
    if max_val <= 0:
        max_val = 1.0

    if fmt is None:
        def fmt(v):
            return f'{int(round(v))}'

    rlabels = row_labels or labels
    clabels = col_labels or labels

    html = ['<div style="overflow:auto;max-height:720px;margin:10px 0;">',
            '<table style="border-collapse:collapse;font-family:monospace;font-size:11px;">']

    # Header row: vertical column labels
    html.append('<tr><td style="min-width:220px;"></td>')
    for lab, disp in zip(labels, clabels):
        if max_label_chars is None or len(disp) <= max_label_chars:
            shown = disp
        else:
            shown = disp[:max_label_chars - 2] + '..'
        html.append(
            f'<td style="writing-mode:vertical-rl;text-align:left;padding:4px 2px;'
            f'font-weight:600;max-width:{cell_size}px;white-space:nowrap;" '
            f'title="{_esc(lab)}">{_esc(shown)}</td>'
        )
    html.append('</tr>')

    for i in range(n):
        html.append(
            f'<tr><td style="padding:2px 8px;font-weight:600;white-space:nowrap;" '
            f'title="{_esc(labels[i])}">{rlabels[i]}</td>'
        )
        for j in range(n):
            val = float(mat[i][j])
            if i == j:
                bg, txt = '#e8e8e8', '#999'
            else:
                bg, txt = color_fn(val, max_val)
            tip = ''
            if tooltips:
                tip = f' title="{_esc(tooltips(labels[i], labels[j], val))}"'
            html.append(
                f'<td style="background:{bg};color:{txt};text-align:center;padding:2px;'
                f'min-width:{cell_size}px;height:{cell_size}px;font-size:10px;"{tip}>'
                f'{fmt(val)}</td>'
            )
        html.append('</tr>')

    html.append('</table></div>')
    return ''.join(html)



def _xml_esc(text):
    """Escape a string for XML text content."""
    return (str(text).replace('&', '&amp;').replace('<', '&lt;')
            .replace('>', '&gt;').replace('"', '&quot;').replace("'", '&apos;'))


def render_heatmap_svg(labels, matrix, color_fn, fmt=None, cell_size=40,
                       label_width=None, label_angle=-60,
                       row_labels=None, col_labels=None):
    """Render the clustered matrix as a standalone SVG."""
    n = len(labels)
    if n == 0:
        return ''

    mat = np.array(matrix, dtype=float)
    max_val = float(mat.max()) if mat.max() > 0 else 1.0

    if fmt is None:
        def fmt(v):
            return f'{int(round(v))}'

    rlabels = [str(x) for x in (row_labels or labels)]
    clabels = [str(x) for x in (col_labels or labels)]

    cell = cell_size
    char_width = 6.5                      # approx char width at font-size 10

    rad = math.radians(abs(label_angle) % 180)
    col_px = max((len(x) for x in clabels), default=10) * char_width
    row_px = max((len(x) for x in rlabels), default=10) * char_width

    # Vertical reach of the rotated column labels sets the header height.
    header_h = max(int(col_px * math.sin(rad) + 30), 80)

    if label_width is not None:
        label_w = label_width
    else:
        label_w = max(160, int(row_px + 14))

    # Horizontal reach of the LAST column label sets the right edge.
    last_cx = label_w + (n - 1) * cell + cell / 2
    w = max(label_w + n * cell + 20,
            int(last_cx + col_px * math.cos(rad) + 12))
    h = header_h + n * cell + 20

    svg = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">\n'
    svg += '<style>text { font-family: monospace; font-size: 11px; }</style>\n'
    svg += f'<rect width="{w}" height="{h}" fill="white"/>\n'

    # Column headers, rotated so long names climb away from the grid
    for j, s_lab in enumerate(clabels):
        cx = label_w + j * cell + cell / 2
        cy = header_h
        svg += (f'<text x="{cx}" y="{cy - 6}" transform="rotate({label_angle},{cx},{cy - 6})" '
                f'text-anchor="start" font-weight="bold" font-size="10">{_xml_esc(s_lab)}</text>\n')

    for i in range(n):
        y = header_h + i * cell
        svg += (f'<text x="{label_w - 6}" y="{y + cell / 2 + 4}" text-anchor="end" '
                f'font-weight="bold" font-size="10">{_xml_esc(rlabels[i])}</text>\n')

        for j in range(n):
            val = float(mat[i][j])
            x = label_w + j * cell
            if i == j:
                bg, txt = '#e8e8e8', '#999'
            else:
                bg, txt = color_fn(val, max_val)
            svg += (f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{bg}" '
                    f'stroke="white" stroke-width="1"/>\n')
            svg += (f'<text x="{x + cell / 2}" y="{y + cell / 2 + 4}" text-anchor="middle" '
                    f'fill="{txt}" font-size="10">{_xml_esc(fmt(val))}</text>\n')

    svg += '</svg>'
    return svg