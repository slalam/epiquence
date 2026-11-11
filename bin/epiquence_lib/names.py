"""Sample-name, assembly-type and plasmid-name parsing helpers."""

import re

# Assembly-type suffixes understood by the pipeline.
ASM_TYPES = ("L", "S", "HL", "HS", "A")


def species_safe(species_name):
    """Filesystem-safe species name: collapse whitespace runs to underscores."""
    return re.sub(r"\s+", "_", species_name or "")


def split_sample_asm(sample_full):
    """Split ``{sample_id}_{asm_type}`` into ``(sample_id, asm_type)``."""
    parts = sample_full.split("_")
    if len(parts) < 2:
        return sample_full, ""
    return "_".join(parts[:-1]), parts[-1]


def parse_gtdbtk_filename(fname):
    """Extract ``(sample_id, asm_type)`` from a GTDB-Tk parsed filename."""
    stem = fname
    if stem.endswith("_gtdbtk_parsed.tsv"):
        stem = stem[: -len("_gtdbtk_parsed.tsv")]
    return split_sample_asm(stem)


# Contig-level plasmid names carry a _CG<n> suffix, e.g. SAMPLE-03_HL_plasmid_AA102_CG3.
CONTIG_SUFFIX_RE = re.compile(r"_(CG\d+)$")


def split_contig_suffix(raw):
    """``..._plasmid_AA102_CG3`` -> ``('..._plasmid_AA102', 'CG3')``."""
    m = CONTIG_SUFFIX_RE.search(raw or '')
    if m:
        return raw[:m.start()], m.group(1)
    return raw, ''


def contig_token(contig_id):
    """``CONTIG_3`` -> ``CG3``."""
    m = re.match(r'^CONTIG_(\d+)$', (contig_id or '').strip(), re.I)
    if m:
        return f"CG{m.group(1)}"
    return re.sub(r'\W+', '', (contig_id or '')) or 'CG0'


def plasmid_display_name(raw):
    """``SAMPLE-03_HL_plasmid_AA102`` -> ``SAMPLE-03_HL_pAA102``
    ``SAMPLE-03_HL_plasmid_AA102_CG3`` -> ``SAMPLE-03_HL_pAA102_CG3``"""
    base, contig = split_contig_suffix(raw)
    disp = re.sub(r"_plasmid_([A-Z]+)", r"_p\1", base)
    return f"{disp}_{contig}" if contig else disp


def split_plasmid_name(raw):
    """Split a plasmid FASTA stem into ``(sample_with_asm, asm_type, mob_cluster)``."""
    raw, _contig = split_contig_suffix(raw)
    m = re.match(r"^(.+)_plasmid_(.+)$", raw)
    if m:
        prefix = m.group(1)
        mob = m.group(2)
        asm_m = re.match(r"^(.+)_(L|S|HL|HS|A)$", prefix)
        if asm_m:
            sample = asm_m.group(1)
            asm = asm_m.group(2)
            return f"{sample}_{asm}", asm, mob
        return prefix, "", mob
    return raw, "", ""