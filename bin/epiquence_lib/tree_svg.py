"""Phylogenetic-tree SVG badge rendering shared by the tree renderers."""

import re
from collections import Counter

# Soft pastel palettes per category - more distinct colors
BADGE_PALETTES = {
    'species': ['#a5d6a7', '#81c784', '#66bb6a', '#4caf50', '#aed581', '#9ccc65', '#7cb342', '#558b2f', '#33691e', '#2e7d32', '#43a047', '#689f38'],
    'mlst':    ['#90caf9', '#42a5f5', '#1e88e5', '#bbdefb', '#64b5f6', '#1976d2', '#0d47a1', '#82b1ff', '#448aff', '#2962ff', '#1565c0', '#0277bd'],
    'amr':     ['#ef9a9a', '#e57373', '#ef5350', '#f48fb1', '#f06292', '#ec407a', '#e53935', '#d32f2f', '#c62828', '#b71c1c', '#ad1457', '#880e4f'],
    'meta':    ['#ffcc80', '#ffb74d', '#ffa726', '#ff9800', '#ffab40', '#ff9100', '#fb8c00', '#f57c00', '#ef6c00', '#e65100', '#d84315', '#bf360c'],
}

BADGE_NEUTRAL = '#b0bec5'


def _text_color(bg_hex):
    """Return dark or light text color based on background luminance."""
    bg = bg_hex.lstrip('#')
    r, g, b = int(bg[:2], 16), int(bg[2:4], 16), int(bg[4:6], 16)
    lum = 0.299 * r + 0.587 * g + 0.114 * b
    return '#333' if lum > 160 else '#fff'


def _build_color_map(tip_data):
    """Build per-value color decisions."""
    cat_counts = {}
    for td in tip_data:
        for cat, val in td['badges']:
            if cat not in cat_counts:
                cat_counts[cat] = Counter()
            cat_counts[cat][val] += 1

    color_map = {}
    for cat, counts in cat_counts.items():
        palette = BADGE_PALETTES.get(cat, BADGE_PALETTES['meta'])
        # Get repeated values sorted for deterministic color assignment
        repeated = sorted([v for v, c in counts.items() if c >= 2])
        for i, val in enumerate(repeated):
            color_map[(cat, val)] = palette[i % len(palette)]
        # Unique values get neutral
        for val, count in counts.items():
            if count < 2:
                color_map[(cat, val)] = BADGE_NEUTRAL
    return color_map


def _add_badges_to_svg(svg_str, tip_data, color_map, initial_width):
    """Post-process SVG to add colored badges, remove clip paths, expand viewBox."""
    import xml.etree.ElementTree as ET
    ET.register_namespace('', 'http://www.w3.org/2000/svg')
    root = ET.fromstring(svg_str)

    # Remove ALL clipPath definitions and clip-path references to prevent trimming
    for clip in list(root.findall('.//{http://www.w3.org/2000/svg}clipPath')):
        for p in root.iter():
            if clip in list(p):
                p.remove(clip)
                break

    for el in root.iter():
        if 'clip-path' in el.attrib:
            del el.attrib['clip-path']
        style = el.get('style', '')
        if 'clip-path' in style:
            style = re.sub(r'clip-path\s*:\s*[^;]+;?', '', style)
            if style.strip():
                el.set('style', style)
            else:
                del el.attrib['style']

    max_x_extent = initial_width

    all_texts = root.findall('.//{http://www.w3.org/2000/svg}text')

    for text_el in all_texts:
        label = (text_el.text or '').strip()
        matching = [td for td in tip_data if td['sample'] == label]
        if not matching or not matching[0]['badges']:
            continue
        td = matching[0]

        x = float(text_el.get('x', 0))
        y = float(text_el.get('y', 0))

        # Account for parent <g> transform offset
        transform_x = 0
        parent = None
        for p in root.iter():
            if text_el in list(p):
                parent = p
                break
        if parent is None:
            continue

        transform = parent.get('transform', '')
        tm = re.search(r'translate\(\s*([\d.eE+-]+)', transform)
        if tm:
            transform_x = float(tm.group(1))

        char_w = 6.5
        badge_x = x + len(label) * char_w + 8

        for cat, val in td['badges']:
            bg = color_map.get((cat, val), BADGE_NEUTRAL)
            fg = _text_color(bg)
            text_w = len(val) * 6.2 + 10
            badge_h = 15
            badge_y = y - 11

            rect = ET.SubElement(parent, '{http://www.w3.org/2000/svg}rect')
            rect.set('x', str(badge_x))
            rect.set('y', str(badge_y))
            rect.set('width', str(text_w))
            rect.set('height', str(badge_h))
            rect.set('rx', '3')
            rect.set('fill', bg)

            badge_text = ET.SubElement(parent, '{http://www.w3.org/2000/svg}text')
            badge_text.set('x', str(badge_x + 5))
            badge_text.set('y', str(y - 1))
            badge_text.set('fill', fg)
            badge_text.set('font-size', '9px')
            badge_text.set('font-family', 'monospace')
            badge_text.text = val

            badge_x += text_w + 4

        if badge_x + transform_x > max_x_extent:
            max_x_extent = badge_x + transform_x

    # Expand SVG viewBox and dimensions to show all content
    max_x_extent += 30
    vb = root.get('viewBox', '')
    if vb:
        parts = vb.split()
        if len(parts) == 4:
            parts[2] = str(max_x_extent)
            root.set('viewBox', ' '.join(parts))
    root.set('width', f'{max_x_extent}px')
    # Remove preserveAspectRatio so SVG doesn't get scaled down to fit container
    if root.get('preserveAspectRatio'):
        del root.attrib['preserveAspectRatio']

    return ET.tostring(root, encoding='unicode')