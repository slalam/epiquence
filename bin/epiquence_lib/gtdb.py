"""GTDB lineage parsing shared by the GTDB-Tk and sylph taxonomy parsers."""
import gzip
import re

RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
_PREFIX = {'d__': 'domain', 'p__': 'phylum', 'c__': 'class', 'o__': 'order',
           'f__': 'family', 'g__': 'genus', 's__': 'species'}
_ACCESSION = re.compile(r'(GC[AF]_\d{9}\.\d+)')


def strip_suffix(value):
    """Remove GTDB polyphyletic suffixes from a genus or species name."""
    value = re.sub(r'_[A-Z]$', '', value)
    return re.sub(r'_[A-Z] ', ' ', value)


def parse_lineage(lineage):
    """``d__...;s__...`` -> ({rank: name}, deepest_name)."""
    taxonomy = {r: '' for r in RANKS}
    deepest = ''
    for part in (lineage or '').split(';'):
        part = part.strip()
        rank = _PREFIX.get(part[:3])
        if not rank:
            continue
        value = part[3:]
        if value and rank in ('genus', 'species'):
            value = strip_suffix(value)
        taxonomy[rank] = value
        if value:
            deepest = value
    return taxonomy, deepest


def accession_of(text):
    """First NCBI assembly accession (GCA_/GCF_ with version) in a string, or ''."""
    m = _ACCESSION.search(text or '')
    return m.group(1) if m else ''


def load_lineages(paths, wanted):
    """Read GTDB taxonomy TSVs (``accession<TAB>lineage``; plain or .gz) and return {accession:
    lineage} for the accessions in ``wanted`` only."""
    wanted = set(wanted)
    twins = {w.replace('GCF_', 'GCA_') if w.startswith('GCF_') else w.replace('GCA_', 'GCF_'): w
             for w in wanted}
    found, via_twin = {}, {}
    for path in paths:
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt', encoding='utf-8') as fh:
            for line in fh:
                key, _, lineage = line.rstrip('\n').partition('\t')
                acc = accession_of(key)
                if acc in wanted and acc not in found:
                    found[acc] = lineage
                elif acc in twins and twins[acc] not in via_twin:
                    via_twin[twins[acc]] = lineage
                if len(found) == len(wanted):
                    return found
    for acc, lineage in via_twin.items():
        found.setdefault(acc, lineage)
    return found
