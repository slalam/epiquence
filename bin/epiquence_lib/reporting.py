"""Shared HTML report helpers."""

import base64

from . import names as _names


def make_download_button(content, filename, icon_title="Download"):
    """Create a small download icon button using base64 data URI."""
    b64 = base64.b64encode(content.encode('utf-8')).decode('utf-8')
    mime = 'image/svg+xml' if filename.endswith('.svg') else 'text/tab-separated-values'
    # Small download arrow icon as inline SVG
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
    return _names.plasmid_display_name(raw)


def split_name(raw):
    """Split plasmid name into (sample_id, asm_type, mob_cluster)."""
    return _names.split_plasmid_name(raw)