"""Scope selection for the cumulative surveillance views."""
import html
import re

ALL_SCOPE = 'all'
UNGROUPED_DIR = 'ungrouped'
_SUPERSEDED = re.compile(r'_O\d+$')
_UNSAFE = re.compile(r'[^A-Za-z0-9._-]')


def scope_dir(group_id):
    """Folder name for a group's views; identical to scopeDir() in subworkflows/surveillance.nf."""
    g = (group_id or '').strip()
    if not g:
        return UNGROUPED_DIR
    safe = _UNSAFE.sub('_', g)
    return 'group_all' if safe == ALL_SCOPE else safe


def base_sample(name):
    """Strip a superseded-version suffix: ``SAMPLE-01_S_O2`` -> ``SAMPLE-01_S``."""
    return _SUPERSEDED.sub('', name)


def make_filter(isolates, group):
    """Return a predicate name -> bool for the requested scope, or None for ``all``."""
    if group is None:
        return None
    want = group.strip()
    members = {s for s, row in isolates.items()
               if (row.get('group_id') or '').strip() == want}
    return lambda name: base_sample(name) in members


def scope_label(isolates, group):
    """Human-readable scope for page headers and titles, HTML-escaped; '' for ``all``."""
    if group is None:
        return ''
    want = group.strip()
    if not want:
        return 'Isolates without a group'
    desc = ''
    for row in isolates.values():
        if (row.get('group_id') or '').strip() == want:
            desc = (row.get('group_description') or '').strip()
            if desc:
                break
    text = f'Group {want}' + (f' | {desc}' if desc else '')
    return html.escape(text)
