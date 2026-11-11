#!/usr/bin/env python3
"""Generate MultiQC configuration files and custom data tables."""
import argparse, os, csv, re, sys, yaml, json
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from epiquence_lib.names import contig_token

def parse_tsv_dir(directory, suffix, callback):
    """Parse all TSV files matching suffix in a directory."""
    if not os.path.exists(directory):
        return
    for f in os.listdir(directory):
        if f.endswith(suffix):
            with open(f'{directory}/{f}', 'r') as fh:
                for row in csv.DictReader(fh, delimiter='\t'):
                    callback(row)

def extract_asm_type(sample_name):
    """Extract assembly type from sample name like SAMPLE-03_HL -> HL."""
    import re
    m = re.search(r'_(L|S|HL|HS)$', sample_name)
    return m.group(1) if m else '-'


def ncbi_genome_link(value):
    """Generate NCBI genome hyperlink for genome references."""
    if value and value != '-':
        accession = str(value).replace('RS_', '').replace('GB_', '')
        return f'<a href="https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/" target="_blank">{value}</a>'
    return value

def ncbi_nuccore_link(value):
    """Generate NCBI nuccore hyperlink for plasmid accessions."""
    if value and value != '-':
        return f'<a href="https://www.ncbi.nlm.nih.gov/nuccore/{value}" target="_blank">{value}</a>'
    return value

# Tokens that must stay upper-case after title-casing: single letters and the roman numerals
# used in streptogramin/pristinamycin naming ("pristinamycin ia" -> "Pristinamycin IA",
# "streptogramin a" -> "Streptogramin A").
_UPPER_TOKENS = {'i', 'ii', 'iii', 'iv', 'v', 'ia', 'ib', 'iia', 'iib'}


def title_case(text):
    """Title-case a ResFinder antimicrobial or class name for display."""
    if not text:
        return text
    parts = []
    for word in str(text).split():
        if word.lower() in _UPPER_TOKENS:
            parts.append(word.upper())
        else:
            # Title-case each hyphen/plus separated component so "amoxicillin+clavulanic acid" and "beta-
            # lactam" both read correctly.
            parts.append(re.sub(r'[A-Za-z]+', lambda m: m.group(0).capitalize(), word))
    return ' '.join(parts)


# Background colours for categorical columns (antimicrobial class, plasmid cluster id).
CATEGORICAL_PALETTE = [
    '#1565c0', '#bf360c', '#2e7d32', '#c7239e', '#4e342e',
    '#2323c7', '#00838f', '#160f57', '#946000', '#880e4f',
    '#915298', '#7b1fa2', '#193e6b', '#6e171b', '#4d4d19',
    '#004d40', '#a45c3d', '#98526b', '#455a64', '#53132d',
    '#bf23c7', '#641268', '#827717', '#b03b46', '#896d43',
    '#3f1c4a', '#554f92', '#283593', '#0277bd', '#4527a0',
    '#ad1457', '#c62828', '#1b5e20', '#00695c', '#6d4c41',
    '#a8320f', '#b34000', '#33691e', '#37474f', '#6a1b9a',
]


def write_yaml(filepath, data):
    """Write YAML with preserved key order."""
    with open(filepath, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

def write_versions(filepath, versions):
    """Write software versions TSV file."""
    with open(filepath, 'w') as f:
        f.write("Software\tVersion\n")
        for sw, ver in versions.items():
            f.write(f"{sw}\t{ver}\n")

AMRFINDER_SUBTYPE_AMR = 'AMR'
AMRFINDER_SUBTYPE_POINT = 'POINT'
# AMRFinderPlus splits hits by Scope: 'core' is NCBI's curated resistance set, 'plus' is what
# --plus adds (efflux pumps, stress response).
AMRFINDER_SCOPE_PLUS = 'plus'


def read_amrfinder_files(amrfinder_path):
    """Read every staged AMRFinderPlus TSV."""
    gene_hits = {}
    point_hits = {}
    all_rows = {}
    if not os.path.exists(amrfinder_path):
        return gene_hits, point_hits, all_rows

    for fn in sorted(os.listdir(amrfinder_path)):
        if not fn.endswith('_amrfinder.tsv'):
            continue
        sample_name = fn[:-len('_amrfinder.tsv')]
        path = os.path.join(amrfinder_path, fn)
        gene_hits[sample_name] = defaultdict(list)
        point_hits[sample_name] = defaultdict(list)
        all_rows[sample_name] = []
        with open(path, newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                # all_rows is the database export and stays unfiltered.
                all_rows[sample_name].append(row)
                if (row.get('Scope') or '').strip().lower() == AMRFINDER_SCOPE_PLUS:
                    continue
                symbol = (row.get('Element symbol') or '').strip()
                contig = (row.get('Contig id') or '').strip().split()[0] if row.get('Contig id') else ''
                if not symbol or not contig:
                    continue
                subtype = (row.get('Subtype') or '').strip().upper()
                if subtype == AMRFINDER_SUBTYPE_AMR:
                    # One row per hit, so the same symbol can repeat on a contig.
                    if symbol not in gene_hits[sample_name][contig]:
                        gene_hits[sample_name][contig].append(symbol)
                elif subtype == AMRFINDER_SUBTYPE_POINT:
                    if symbol not in point_hits[sample_name][contig]:
                        point_hits[sample_name][contig].append(symbol)
    return gene_hits, point_hits, all_rows


RESFINDER_COLUMNS = ('Antimicrobial', 'Class', 'WGS-predicted phenotype', 'Match',
                     'Genetic background')
RESFINDER_EXCLUDED_CLASSES = {'under_development'}


def read_resfinder_pheno(pheno_path):
    """Read every staged ResFinder pheno_table.txt."""
    results = {}
    if not os.path.exists(pheno_path):
        return results

    for fn in sorted(os.listdir(pheno_path)):
        if not fn.endswith('_pheno_table.txt'):
            continue
        sample_name = fn[:-len('_pheno_table.txt')]
        rows = []
        with open(os.path.join(pheno_path, fn)) as f:
            for line in f:
                line = line.rstrip('\n')
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.split('\t')
                rows.append({c: (parts[i].strip() if i < len(parts) else '')
                             for i, c in enumerate(RESFINDER_COLUMNS)})
        results[sample_name] = rows
    return results


def merge_preserving_order(primary, secondary):
    """Union two gene/mutation lists, primary first, de-duplicated case-insensitively."""
    merged, seen = [], set()
    for value in list(primary) + list(secondary):
        text = (value or '').strip()
        if not text or text == '-':
            continue
        key = text.lower()
        if key not in seen:
            seen.add(key)
            merged.append(text)
    return merged


def split_list_field(value):
    """'a; b, c' -> ['a', 'b', 'c']."""
    text = (value or '').strip()
    if not text or text == '-':
        return []
    return [p.strip() for p in re.split(r'[;,]', text) if p.strip()]


def read_resfinder_hits(hits_path):
    """Read every staged ResFinder_results_tab.txt."""
    by_contig, raw = {}, {}
    if not os.path.exists(hits_path):
        return by_contig, raw
    for fn in sorted(os.listdir(hits_path)):
        if not fn.endswith('_resfinder_hits.txt'):
            continue
        sample_name = fn[:-len('_resfinder_hits.txt')]
        by_contig[sample_name] = defaultdict(list)
        raw[sample_name] = []
        with open(os.path.join(hits_path, fn), newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                raw[sample_name].append(row)
                gene = (row.get('Resistance gene') or '').strip()
                contig_full = (row.get('Contig') or '').strip()
                contig = contig_full.split()[0] if contig_full else ''
                if gene and contig and gene not in by_contig[sample_name][contig]:
                    by_contig[sample_name][contig].append(gene)
    return by_contig, raw


def parse_genetic_background(value):
    """Pull the gene names out of a ResFinder "Genetic background" cell."""
    genes = []
    for entry in (value or '').split(','):
        entry = entry.strip()
        if not entry:
            continue
        gene = entry.rsplit(' (', 1)[0].strip() if ' (' in entry else entry
        if gene and gene not in genes:
            genes.append(gene)
    return '; '.join(genes)


POINTFINDER_COLUMNS = ('Mutation', 'Nucleotide change', 'Amino acid change',
                       'Resistance', 'PMID')


def format_point_mutation(mutation):
    """``gyrA p.S91F`` -> ``gyrA_S91F``."""
    m = re.match(r'^\s*(\S+)\s+(?:[a-z]\.)?(\S+)\s*$', str(mutation or ''))
    if m:
        return f"{m.group(1)}_{m.group(2)}"
    return str(mutation or '').strip()


def read_pointfinder(pointfinder_path):
    """Read every staged PointFinder_results.txt."""
    results = {}
    if not os.path.exists(pointfinder_path):
        return results
    for fn in sorted(os.listdir(pointfinder_path)):
        if not fn.endswith('_pointfinder.txt'):
            continue
        sample_name = fn[:-len('_pointfinder.txt')]
        rows = []
        with open(os.path.join(pointfinder_path, fn)) as f:
            header = f.readline().rstrip('\n').split('\t')
            for line in f:
                if not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                rows.append({c: (parts[i].strip() if i < len(parts) else '')
                             for i, c in enumerate(header or POINTFINDER_COLUMNS)})
        results[sample_name] = rows
    return results


def resfinder_reportable(row):
    """Exclude unnameable drugs and not-yet-validated classes from the report."""
    drug = row.get('Antimicrobial', '').strip().lower()
    cls = row.get('Class', '').strip().lower()
    if not drug or drug.startswith('unknown'):
        return False
    return cls not in RESFINDER_EXCLUDED_CLASSES


def build_resistance_marker_tables(amrfinder_path, mob_dirs_path):
    """Build plasmid and chromosome resistance marker tables by joining AMRFinderPlus hits with
    MOB-recon contig reports and mobtyper results."""
    plasmid_markers = {}
    chromo_markers = {}

    gene_hits, point_hits, amrfinder_rows = read_amrfinder_files(amrfinder_path)

    if not os.path.exists(mob_dirs_path):
        return plasmid_markers, chromo_markers, {}, gene_hits, point_hits, {}, amrfinder_rows

    # Per-sample: contig_id -> {molecule_type, cluster_id} from contig_report.txt
    sample_contig_map = {}
    # Per-sample: cluster_id -> plasmid metadata from mobtyper_results.txt
    sample_mobtyper = {}

    for d in os.listdir(mob_dirs_path):
        mob_dir = f"{mob_dirs_path}/{d}"
        if not os.path.isdir(mob_dir):
            continue
        sample_name = d.replace('_mob_recon', '')

        # Parse contig_report.txt - extract first word of contig_id, molecule_type, cluster_id
        contig_report = f"{mob_dir}/contig_report.txt"
        if os.path.exists(contig_report):
            sample_contig_map[sample_name] = {}
            with open(contig_report, 'r') as f:
                for row in csv.DictReader(f, delimiter='\t'):
                    cid_full = row.get('contig_id', '').strip()
                    cid = cid_full.split()[0] if cid_full else ''
                    if not cid:
                        continue
                    # Parse circularity from FASTA header: Unicycler (S/HS): circular=true.
                    cid_lower = cid_full.lower()
                    is_circular = ('circular=true' in cid_lower
                                   or 'circular=y' in cid_lower
                                   or 'rotated=true' in cid_lower)
                    circularity = 'circular' if is_circular else 'linear'
                    sample_contig_map[sample_name][cid] = {
                        'molecule_type': row.get('molecule_type', '').strip(),
                        'cluster_id': row.get('primary_cluster_id', '-').strip(),
                        'contig_length': row.get('contig_length', row.get('size', '0')).strip(),
                        'circularity': circularity,
                    }

        # Parse mobtyper_results.txt for per-cluster plasmid metadata
        mobtyper_file = f"{mob_dir}/mobtyper_results.txt"
        if os.path.exists(mobtyper_file):
            sample_mobtyper[sample_name] = {}
            with open(mobtyper_file, 'r') as f:
                for row in csv.DictReader(f, delimiter='\t'):
                    cluster_id = row.get('primary_cluster_id', '-').strip()
                    if not cluster_id or cluster_id == '-':
                        continue
                    sample_mobtyper[sample_name][cluster_id] = {
                        'num_contigs': row.get('num_contigs', '-').strip(),
                        'size': row.get('size', '-').strip(),
                        'secondary_cluster_id': row.get('secondary_cluster_id', '-').strip(),
                        'rep_type': row.get('rep_type(s)', '-').strip(),
                        'relaxase_type': row.get('relaxase_type(s)', '-').strip(),
                        'mpf_type': row.get('mpf_type', '-').strip(),
                        'orit_type': row.get('orit_type(s)', '-').strip(),
                        'predicted_mobility': row.get('predicted_mobility', '-').strip(),
                        'mash_nearest': row.get('mash_nearest_neighbor', '-').strip(),
                    }

    # Join: group by sample, then by cluster_id for plasmids or chromosome
    all_samples = sorted(set(list(sample_contig_map.keys()) + list(gene_hits.keys())))
    for sample_name in all_samples:
        contig_map = sample_contig_map.get(sample_name, {})
        res_contigs = gene_hits.get(sample_name, {})
        mobtyper = sample_mobtyper.get(sample_name, {})

        plasmid_clusters = defaultdict(lambda: {'contigs': [], 'genes': set(), 'circularity': []})
        chromo_genes = defaultdict(list)

        for cid, info in contig_map.items():
            genes = res_contigs.get(cid, [])
            if info['molecule_type'] == 'plasmid':
                cluster = info['cluster_id'] if info['cluster_id'] != '-' else 'unknown'
                plasmid_clusters[cluster]['contigs'].append(cid)
                plasmid_clusters[cluster]['genes'].update(genes)
                plasmid_clusters[cluster]['circularity'].append(info.get('circularity', '-'))
            elif info['molecule_type'] == 'chromosome' and genes:
                chromo_genes[cid] = genes

        # Build plasmid markers rows (metadata from mobtyper_results.txt)
        for cluster_id, cluster_info in sorted(plasmid_clusters.items()):
            row_key = f"{sample_name}:{cluster_id}"
            join_or_dash = lambda s: '; '.join(sorted(s)) if s else '-'
            meta = mobtyper.get(cluster_id, {})
            mash_val = meta.get('mash_nearest', '-')
            mash_display = ncbi_nuccore_link(mash_val) if mash_val and mash_val != '-' else '-'

            num_contigs_val = meta.get('num_contigs', str(len(cluster_info['contigs'])))
            try:
                num_contigs_int = int(num_contigs_val)
            except (ValueError, TypeError):
                num_contigs_int = len(cluster_info['contigs'])

            size_val = meta.get('size', '0')
            try:
                size_kb = int(size_val) / 1000.0
            except (ValueError, TypeError):
                size_kb = 0.0

            # Cluster is circular if any contig in the cluster is circular
            cluster_circular = 'circular' if any(c == 'circular' for c in cluster_info['circularity']) else 'linear'

            plasmid_markers[row_key] = {
                'Sample': sample_name,
                'Cluster_ID': cluster_id,
                'Num_Contigs': num_contigs_int,
                'Size_K': size_kb,
                'Circularity': cluster_circular,
                'Secondary_Cluster_ID': meta.get('secondary_cluster_id', '-'),
                'Num_Resistance_Genes': len(cluster_info['genes']),
                'Resistance_Genes': join_or_dash(cluster_info['genes']),
                'Rep_Type': meta.get('rep_type', '-'),
                'Relaxase_Type': meta.get('relaxase_type', '-'),
                'MPF_Type': meta.get('mpf_type', '-'),
                'OriT_Type': meta.get('orit_type', '-'),
                'Predicted_Mobility': meta.get('predicted_mobility', '-'),
                'Closest_Match': mash_display,
                'Contigs': '; '.join(sorted(cluster_info['contigs'])),
            }

        # Build chromosome markers rows
        if chromo_genes:
            all_genes = set()
            all_contigs = []
            total_chromo_size = 0
            circularity_statuses = []
            for cid, genes in sorted(chromo_genes.items()):
                all_contigs.append(cid)
                all_genes.update(genes)
                circularity_statuses.append(contig_map.get(cid, {}).get('circularity', '-'))
                try:
                    total_chromo_size += int(contig_map.get(cid, {}).get('contig_length', 0))
                except (ValueError, TypeError):
                    pass
            chromo_circular = 'circular' if any(c == 'circular' for c in circularity_statuses) else 'linear'
            chromo_markers[sample_name] = {
                'Contigs': '; '.join(all_contigs),
                'Size_K': total_chromo_size / 1000.0,
                'Circularity': chromo_circular,
                'Num_Resistance_Genes': len(all_genes),
                'Resistance_Genes': '; '.join(sorted(all_genes)),
            }

    return (plasmid_markers, chromo_markers, sample_contig_map,
            gene_hits, point_hits, sample_mobtyper, amrfinder_rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--title', default='Hybrid Assembly Pipeline Report')
    parser.add_argument('--subtitle', default='')
    parser.add_argument('--intro', default='Comprehensive genomic surveillance pipeline for bacterial isolates, covering read QC, hybrid/long-read/short-read assembly, taxonomic classification, MLST typing, antimicrobial resistance profiling, plasmid identification, horizontal gene transfer detection, and outbreak clustering.')
    parser.add_argument('--analysis-name', default='')
    parser.add_argument('--workflow-id', default='')
    parser.add_argument('--workflow-description', default='')
    parser.add_argument('--group-id', default='GP001')
    parser.add_argument('--group-description', default='Healthcare-Associated Infections')
    parser.add_argument('--workflow-timestamp', default='')
    parser.add_argument('--taxonomy-tool', default='gtdbtk', choices=['gtdbtk', 'sylph', 'both'],
                        help='Species engine used for this run (labels, versions, agreement columns)')
    args = parser.parse_args()

    # Generate default workflow_id if not provided
    from datetime import datetime
    if args.workflow_id:
        run_id = args.workflow_id
    else:
        run_id = 'W' + datetime.now().strftime('%m%d%y%H%M%S')

    group_id = args.group_id
    group_desc = args.group_description
    workflow_desc = args.workflow_description
    workflow_timestamp = args.workflow_timestamp or datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # ── Data collection ──
    samples = defaultdict(dict)
    mlst_data = {}
    gtdbtk_data = {}
    # 'both' mode: sylph's call per sample, compared against the authoritative GTDB-Tk call.
    sylph_compare = {}
    taxonomy_tool = args.taxonomy_tool
    species_source = {'gtdbtk': 'GTDB species (GTDB-Tk)',
                      'sylph': 'GTDB species (sylph, assembly containment ANI)',
                      'both': 'GTDB species (GTDB-Tk; sylph shown for comparison)'}[taxonomy_tool]
    mob_data = {}

    def parse_quality(row):
        sample = row['Sample']
        try: samples[sample]['Median_Phred'] = float(row.get('Median_Phred', 0))
        except (ValueError, TypeError): samples[sample]['Median_Phred'] = 0.0
        # Pick up filter stats from fastp quality files
        if 'Reads_Before_Filter' in row:
            try:
                samples[sample]['Reads_Before_Filter'] = int(row.get('Reads_Before_Filter', 0))
                samples[sample]['Reads_After_Filter'] = int(row.get('Reads_After_Filter', 0))
                samples[sample]['Reads_After_Filter_K'] = int(row.get('Reads_After_Filter', 0)) / 1000.0
                samples[sample]['Pct_Reads_Passed'] = float(row.get('Pct_Reads_Passed', 0))
            except (ValueError, TypeError):
                pass
        if 'Bases_Before_Filter' in row:
            try:
                samples[sample]['Bases_After_Filter_Mb'] = int(row.get('Bases_After_Filter', 0)) / 1000000.0
                samples[sample]['Pct_Bases_Passed'] = float(row.get('Pct_Bases_Passed', 0))
            except (ValueError, TypeError):
                pass
    parse_tsv_dir('quality', '_quality_mqc.tsv', parse_quality)

    # Parse fastplong filter stats
    def parse_filterstats(row):
        sample = row['Sample']
        try:
            samples[sample]['Reads_Before_Filter'] = int(row.get('Reads_Before_Filter', 0))
            samples[sample]['Reads_After_Filter'] = int(row.get('Reads_After_Filter', 0))
            samples[sample]['Reads_After_Filter_K'] = int(row.get('Reads_After_Filter', 0)) / 1000.0
            samples[sample]['Pct_Reads_Passed'] = float(row.get('Pct_Reads_Passed', 0))
        except (ValueError, TypeError):
            pass
        if 'Bases_Before_Filter' in row:
            try:
                samples[sample]['Bases_After_Filter_Mb'] = int(row.get('Bases_After_Filter', 0)) / 1000000.0
                samples[sample]['Pct_Bases_Passed'] = float(row.get('Pct_Bases_Passed', 0))
            except (ValueError, TypeError):
                pass
    parse_tsv_dir('filterstats', '_filterstats_mqc.tsv', parse_filterstats)

    def parse_seqkit(row):
        sample = row['Sample']
        samples[sample]['Total_Reads_K'] = int(row.get('Total_Reads', 0)) / 1000.0
        samples[sample]['Total_Bases_Mb'] = int(row.get('Total_Bases', 0)) / 1000000.0
        samples[sample]['Avg_Read_Len'] = float(row.get('Avg_Len', 0))
        samples[sample]['GC_pct'] = float(row.get('GC_pct', 0))
    parse_tsv_dir('seqkit', '_seqkit_mqc.tsv', parse_seqkit)

    def parse_coverage(row):
        sample = row['Sample']
        samples[sample]['Mean_Coverage'] = int(float(row.get('Mean_Coverage', 0)))
        samples[sample]['Median_Coverage'] = int(float(row.get('Median_Coverage', 0)))
        for t in [20, 30, 40, 50, 100]:
            samples[sample][f'Pct_ge_{t}x'] = float(row.get(f'Pct_ge_{t}x', 0))
    parse_tsv_dir('coverage', '_coverage_mqc.tsv', parse_coverage)

    def parse_mlst(row):
        sample = row['Sample']
        samples[sample]['MLST_Scheme'] = row.get('MLST_Scheme', '-')
        samples[sample]['MLST_ST'] = row.get('MLST_ST', '-')
        mlst_data[sample] = {
            'Scheme': row.get('MLST_Scheme', '-'),
            'ST': row.get('MLST_ST', '-'),
            **{f'Locus_{i}': row.get(f'Locus_{i}', '-') for i in range(1, 8)}
        }
    parse_tsv_dir('mlst', '_mlst_parsed.tsv', parse_mlst)

    def parse_gtdbtk(row):
        sample = row['Sample']
        species = row.get('species', '')
        genus = row.get('genus', '')
        samples[sample]['Isolate_Identity'] = species or genus or row.get('classification', '-')
        samples[sample]['Closest_Genome_Ref'] = row.get('closest_genome_reference', '-')
        samples[sample]['Closest_Genome_ANI'] = row.get('closest_genome_ani', '-')
        samples[sample]['Closest_Genome_AF'] = row.get('closest_genome_af', '-')
        gtdbtk_data[sample] = {
            'Domain': row.get('domain', '-'), 'Phylum': row.get('phylum', '-'),
            'Class': row.get('class', '-'), 'Order': row.get('order', '-'),
            'Family': row.get('family', '-'), 'Genus': row.get('genus', '-'),
            'Species': row.get('species', '-'),
            'Closest_Reference': row.get('closest_genome_reference', '-'),
            'ANI': row.get('closest_genome_ani', '-'),
            'AF': row.get('closest_genome_af', '-'),
            # Absent from GTDB-Tk files; present from parse_sylph.py.
            'Method': row.get('method') or 'gtdbtk',
            'Call': row.get('call_status') or '-',
            'Second_Species': row.get('second_species') or '-',
            'Second_ANI': row.get('second_ani') or '-',
        }
    parse_tsv_dir('gtdbtk', '_gtdbtk_parsed.tsv', parse_gtdbtk)

    def parse_sylph_compare(row):
        sylph_compare[row['Sample']] = {
            'species': row.get('species') or '',
            'genus': row.get('genus') or '',
            'ani': row.get('closest_genome_ani') or '',
            'call': row.get('call_status') or '',
        }
    parse_tsv_dir('gtdbtk', '_sylph_parsed.tsv', parse_sylph_compare)

    def taxonomy_agreement(sample_name):
        """'yes' / 'no' / 'genus only' / '-' for GTDB-Tk vs sylph in 'both' mode."""
        sy = sylph_compare.get(sample_name)
        gt = gtdbtk_data.get(sample_name, {})
        if not sy:
            return '-'
        g_sp = gt.get('Species') if gt.get('Species') not in (None, '', '-') else ''
        if g_sp and sy['species']:
            return 'yes' if g_sp == sy['species'] else 'no'
        g_ge = gt.get('Genus') if gt.get('Genus') not in (None, '', '-') else ''
        if g_ge and sy['genus']:
            return 'genus only' if g_ge == sy['genus'] else 'no'
        return '-'

    def parse_quast(row):
        sample = row['Sample']
        samples[sample]['L50'] = int(row.get('L50', 0))
        samples[sample]['Largest_Contig_Mbp'] = float(row.get('Largest_Contig_Mbp', 0))
        samples[sample]['Num_Contigs'] = int(row.get('Num_Contigs', row.get('L50', 0)))
        samples[sample]['N50_kbp'] = float(row.get('N50_kbp', 0))
        samples[sample]['Total_Length_Mbp'] = float(row.get('Total_Length_Mbp', 0))
    parse_tsv_dir('quast', '_quast_mqc.tsv', parse_quast)

    def parse_amrfinder(row):
        sample = row['Sample']
        # The report counts core genes only, so the generic AMR_Genes fields which feed the General
        # Statistics column and general_statistics.tsv carry the core set.
        samples[sample]['AMR_Genes'] = int(row.get('Num_Core_AMR_Genes', 0))
        samples[sample]['AMR_Gene_List'] = row.get('Core_AMR_Genes', '-')
        samples[sample]['Core_AMR_Genes'] = int(row.get('Num_Core_AMR_Genes', 0))
        samples[sample]['Core_AMR_Gene_List'] = row.get('Core_AMR_Genes', '-')
        samples[sample]['Plus_AMR_Genes'] = int(row.get('Num_Plus_AMR_Genes', 0))
        samples[sample]['Plus_AMR_Gene_List'] = row.get('Plus_AMR_Genes', '-')
        samples[sample]['Point_Mutations'] = int(row.get('Num_Point_Mutations', 0))
        samples[sample]['Point_Mutation_List'] = row.get('Point_Mutations', '-')
    parse_tsv_dir('amrfinder', '_amrfinder_parsed.tsv', parse_amrfinder)

    def parse_resfinder(row):
        sample = row['Sample']
        samples[sample]['Resistant_Antimicrobials'] = int(row.get('Num_Resistant', 0))
        samples[sample]['Resistant_Antimicrobial_List'] = row.get('Resistant_Antimicrobials', '-')
    parse_tsv_dir('resfinder', '_resfinder_parsed.tsv', parse_resfinder)

    def parse_mob(row):
        sample = row['Sample']
        num_plasmid_clusters = int(row.get('Num_Plasmids', 0))
        num_plasmid_contigs = int(row.get('Num_Plasmid_Contigs', 0))
        num_chromosomes = int(row.get('Num_Chromosomes', 0))
        total_size_bp = int(row.get('Total_Plasmid_Size_bp', 0))
        samples[sample]['Plasmid_Clusters'] = num_plasmid_clusters
        samples[sample]['Plasmid_Cluster_List'] = row.get('Primary_Cluster_IDs', '-')
        samples[sample]['Plasmid_Contigs'] = num_plasmid_contigs
        samples[sample]['Chromosomal_Contigs'] = num_chromosomes
        mob_data[sample] = {
            'Num_Plasmid_Clusters': num_plasmid_clusters,
            'Num_Plasmid_Contigs': num_plasmid_contigs,
            'Num_Chromosomal_Contigs': num_chromosomes,
            'Total_Plasmid_Size_K': total_size_bp / 1000.0,
            'Replicon_Types': row.get('Replicon_Types', '-'),
            'Relaxase_Types': row.get('Relaxase_Types', '-'),
            'Mash_Nearest': row.get('Mash_Nearest', '-'),
            'Primary_Cluster_IDs': row.get('Primary_Cluster_IDs', '-'),
            'Secondary_Cluster_IDs': row.get('Secondary_Cluster_IDs', '-')
        }
    parse_tsv_dir('mob', '_mob_parsed.tsv', parse_mob)

    # ── Build resistance marker tables ──
    (plasmid_markers, chromo_markers, sample_contig_map, sample_amr_genes,
     sample_point_mutations, sample_mobtyper, amrfinder_rows) = \
        build_resistance_marker_tables('amrfinder_files', 'mob_dirs')
    resfinder_pheno = read_resfinder_pheno('pheno_files')
    pointfinder_hits = read_pointfinder('pointfinder_files')
    resfinder_by_contig, resfinder_raw = read_resfinder_hits('resfinder_hit_files')

    # ── Merge the two AMR callers for the General Statistics view ──
    # The per-tool tables stay untouched; only the summary fields are combined, so a gene found by
    # exactly one tool is not silently lost.
    for sample_name in list(samples.keys()):
        amrfinder_genes = split_list_field(samples[sample_name].get('AMR_Gene_List'))
        resfinder_genes = []
        for row in resfinder_pheno.get(sample_name, []):
            if resfinder_reportable(row):
                resfinder_genes.extend(
                    split_list_field(parse_genetic_background(row.get('Genetic background', ''))))
        merged_genes = merge_preserving_order(amrfinder_genes, resfinder_genes)

        amrfinder_muts = split_list_field(samples[sample_name].get('Point_Mutation_List'))
        pointfinder_muts = [format_point_mutation(r.get('Mutation', ''))
                            for r in pointfinder_hits.get(sample_name, [])]
        merged_muts = merge_preserving_order(amrfinder_muts, pointfinder_muts)

        samples[sample_name]['AMRFinder_Gene_List'] = '; '.join(amrfinder_genes) or '-'
        samples[sample_name]['ResFinder_Gene_List'] = '; '.join(
            merge_preserving_order(resfinder_genes, [])) or '-'
        samples[sample_name]['Merged_AMR_Gene_List'] = '; '.join(merged_genes) or '-'
        samples[sample_name]['Merged_AMR_Genes'] = len(merged_genes)
        samples[sample_name]['Merged_Point_Mutation_List'] = '; '.join(merged_muts) or '-'
        samples[sample_name]['Merged_Point_Mutations'] = len(merged_muts)

    # ── Write plasmid HGT metadata TSV ──
    # Maps a plasmid FASTA name to its size, resistance genes and species.
    with open('plasmid_hgt_metadata.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['plasmid_name', 'size_bp', 'resistance_genes', 'species'])

        # Contig-level rows, from the per-contig data build_resistance_marker_tables already has but
        # aggregates away.
        for sample_asm, contig_map in sorted(sample_contig_map.items()):
            # localised_hits is assembled further down, so merge the two callers' per-contig maps here
            # rather than reordering the file.
            hits = defaultdict(list)
            for contig, genes in sample_amr_genes.get(sample_asm, {}).items():
                hits[contig].extend(genes)
            known = {g.lower() for genes in sample_amr_genes.get(sample_asm, {}).values()
                     for g in genes}
            for contig, genes in resfinder_by_contig.get(sample_asm, {}).items():
                for gene in genes:
                    if gene.lower() not in known and gene not in hits[contig]:
                        hits[contig].append(gene)
            species_name = samples.get(sample_asm, {}).get('Isolate_Identity', '')
            for contig_id, info in sorted(contig_map.items()):
                if info.get('molecule_type') != 'plasmid':
                    continue
                cluster_id = info.get('cluster_id', '-')
                if not cluster_id or cluster_id == '-':
                    cluster_id = 'unknown'
                token = contig_token(contig_id)
                try:
                    size_bp = int(info.get('contig_length', 0))
                except (ValueError, TypeError):
                    size_bp = 0
                genes = ';'.join(sorted(hits.get(contig_id, [])))
                w.writerow([f"{sample_asm}_plasmid_{cluster_id}_{token}",
                            size_bp, genes, species_name])

        for row_key, pm in plasmid_markers.items():
            # row_key = "SAMPLE-03_HL:AA102" -> plasmid_name = "SAMPLE-03_HL_plasmid_AA102"
            parts = row_key.split(':', 1)
            if len(parts) == 2:
                sample_asm, cluster_id = parts
                plasmid_name = f"{sample_asm}_plasmid_{cluster_id}"
                size_bp = int(pm.get('Size_K', 0) * 1000)
                genes = pm.get('Resistance_Genes', '-')
                if genes == '-':
                    genes = ''
                else:
                    genes = genes.replace('; ', ';')
                species = samples.get(sample_asm, {}).get('Isolate_Identity', '')
                w.writerow([plasmid_name, size_bp, genes, species])

    # ── Write chromosomal AMR TSV for outbreak tree annotation ──
    with open('chromosomal_amr.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Sample', 'Chromosomal_AMR_Genes'])
        for sample_name, cm in sorted(chromo_markers.items()):
            genes = cm.get('Resistance_Genes', '')
            if genes and genes != '-':
                w.writerow([sample_name, genes])

    # ── Store chromosomal and plasmid AMR gene lists in samples dict ──
    for sample_name, cm in chromo_markers.items():
        genes = cm.get('Resistance_Genes', '')
        if genes and genes != '-':
            samples[sample_name]['Chromosomal_AMR_Gene_List'] = genes

    # Aggregate plasmid AMR genes per sample from plasmid_markers
    plasmid_amr_by_sample = defaultdict(set)
    for row_key, pm in plasmid_markers.items():
        parts = row_key.split(':', 1)
        if len(parts) == 2:
            sample_asm = parts[0]
            genes = pm.get('Resistance_Genes', '-')
            if genes and genes != '-':
                for g in genes.replace('; ', ';').split(';'):
                    g = g.strip()
                    if g:
                        plasmid_amr_by_sample[sample_asm].add(g)
    for sample_name, genes in plasmid_amr_by_sample.items():
        samples[sample_name]['Plasmid_AMR_Gene_List'] = '; '.join(sorted(genes))

    # ── Build virulence factor data from VirulenceFinder output ──
    sample_virulence = {}  # sample_name -> list of dicts (full hit details)
    sample_virulence_factors = {}  # sample_name -> sorted unique factor names
    if os.path.exists('virfinder_dirs'):
        for d in sorted(os.listdir('virfinder_dirs')):
            vf_dir = f"virfinder_dirs/{d}"
            results_file = f"{vf_dir}/results_tab.tsv"
            if not os.path.isdir(vf_dir) or not os.path.exists(results_file):
                continue
            sample_name = d.replace('_virulencefinder', '')
            hits = []
            factor_names = []
            with open(results_file) as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    vf = row.get('Virulence factor', '').strip()
                    if vf:
                        factor_names.append(vf)
                        hits.append({
                            'database': row.get('Database', ''),
                            'virulence_factor': vf,
                            'identity': row.get('Identity', ''),
                            'query_template_length': row.get('Query / Template length', ''),
                            'contig': row.get('Contig', ''),
                            'position': row.get('Position in contig', ''),
                            'protein_function': row.get('Protein function', ''),
                            'accession': row.get('Accession number', ''),
                        })
            sample_virulence[sample_name] = hits
            sample_virulence_factors[sample_name] = sorted(set(factor_names))

    # Store virulence summary in samples dict for general_statistics
    for sample_name, factors in sample_virulence_factors.items():
        if sample_name in samples:
            samples[sample_name]['Virulence_Factor_Count'] = len(factors)
            samples[sample_name]['Virulence_Factor_List'] = json.dumps(factors) if factors else '[]'
            samples[sample_name]['Virulence_Factor_Display'] = '; '.join(factors) if factors else '-'

    # ── Propagate QC values from L->HL, S->HS for hybrid assemblies ──
    qc_fields = ['Total_Reads_K', 'Total_Bases_Mb', 'Avg_Read_Len', 'GC_pct', 'Median_Phred',
                 'Reads_After_Filter_K', 'Pct_Reads_Passed', 'Reads_Before_Filter', 'Reads_After_Filter',
                 'Bases_After_Filter_Mb', 'Pct_Bases_Passed']
    for sample_name in list(samples.keys()):
        # HL gets values from L, HS gets values from S
        if sample_name.endswith('_HL'):
            base = sample_name[:-3] + '_L'
            for field in qc_fields:
                if field not in samples[sample_name] and field in samples.get(base, {}):
                    samples[sample_name][field] = samples[base][field]
        elif sample_name.endswith('_HS'):
            base = sample_name[:-3] + '_S'
            for field in qc_fields:
                if field not in samples[sample_name] and field in samples.get(base, {}):
                    samples[sample_name][field] = samples[base][field]

    # ── Convert ANI/AF to numeric for coloring ──
    for sample_name in samples:
        for field in ['Closest_Genome_ANI', 'Closest_Genome_AF']:
            val = samples[sample_name].get(field, '-')
            if val and val != '-':
                try:
                    samples[sample_name][field] = float(val)
                except (ValueError, TypeError):
                    pass

    # ── General Stats YAML ──
    # General Statistics columns, in display order.
    general_stats_configs = [
        # id, title, description, format, scale, suffix, hidden
        ('Isolate_Identity', 'Species', 'Taxonomic classification', '{}', False, '', False),
        ('MLST_Scheme', 'MLST Scheme', 'MLST Scheme detected', '{}', False, '', False),
        ('MLST_ST', 'MLST', 'MLST Sequence Type', '{}', False, '', False),
        ('Merged_AMR_Genes', 'AMR Genes', 'Number of acquired AMR genes detected', '{:,}', 'OrRd', '', False),
        # Title is "Plasmids" so it does not collide with Plasmid_Cluster_List below, which keeps
        # "Plasmid Clusters": both are visible by default and two identically-titled columns would be
        # indistinguishable in the header.
        ('Plasmid_Clusters', 'Plasmids', 'Number of plasmid clusters detected', '{:,}', 'BuPu', '', False),
        ('Median_Phred', 'Median Q', 'Median read quality (Phred score)', '{:.1f}', 'RdYlGn', '', False),
        ('Mean_Coverage', 'Mean Cov', 'Mean coverage depth', '{:,}', 'RdYlGn', 'x', True),
        ('Median_Coverage', 'Median Cov', 'Median coverage depth', '{:,}', 'RdYlGn', 'x', False),
        ('Plasmid_Cluster_List', 'Plasmid Clusters', 'Primary plasmid cluster identifiers', '{}', False, '', False),
        ('Merged_AMR_Gene_List', 'AMR Genes List', 'Names of acquired AMR genes detected', '{}', False, '', False),
        ('Merged_Point_Mutation_List', 'AMR Point Mutations', 'Resistance-conferring point mutations detected', '{}', False, '', False),
        ('Total_Reads_K', 'Total Reads (K)', 'Total number of reads in thousands', '{:,.1f}', 'Blues', ' K', False),
        ('Reads_After_Filter_K', 'Filtered Reads (K)', 'Total reads after filtering in thousands', '{:,.1f}', 'Greens', ' K', False),
        ('Pct_Reads_Passed', '% Reads Passed', 'Percentage of reads passed filter', '{:.1f}', 'RdYlGn', '%', False),
        ('Total_Bases_Mb', 'Total Bases (Mb)', 'Total number of bases in megabases', '{:,.1f}', 'Greens', ' Mb', True),
        ('Bases_After_Filter_Mb', 'Filtered Bases (Mb)', 'Total bases after filtering in megabases', '{:,.1f}', 'Greens', ' Mb', True),
        ('Pct_Bases_Passed', '% Bases Passed', 'Percentage of bases passed filter', '{:.1f}', 'RdYlGn', '%', True),
        ('Avg_Read_Len', 'Avg Len', 'Average read length', '{:,.0f}', 'Purples', ' bp', False),
        ('Total_Length_Mbp', 'Assembly Length', 'Total assembly length in Mbp', '{:.3f}', 'Greens', ' Mbp', False),
        ('N50_kbp', 'N50', 'N50 in kilobases', '{:,.1f}', 'RdYlGn', ' Kb', False),
        ('L50', 'L50', 'L50 - number of contigs comprising 50% of assembly', '{:,}', 'Blues', '', False),
        ('Largest_Contig_Mbp', 'Largest Contig', 'Largest contig size in Mbp', '{:.3f}', 'Greens', ' Mbp', False),
        ('GC_pct', 'GC %', 'GC content percentage', '{:.1f}', 'Oranges', '%', True),
        ('Pct_ge_20x', '\u226520x', 'Percent bases \u226520x coverage', '{:.1f}', 'RdYlGn', '%', True),
        ('Pct_ge_30x', '\u226530x', 'Percent bases \u226530x coverage', '{:.1f}', 'RdYlGn', '%', True),
        ('Pct_ge_40x', '\u226540x', 'Percent bases \u226540x coverage', '{:.1f}', 'RdYlGn', '%', True),
        ('Pct_ge_50x', '\u226550x', 'Percent bases \u226550x coverage', '{:.1f}', 'RdYlGn', '%', True),
        ('Pct_ge_100x', '\u2265100x', 'Percent bases \u2265100x coverage', '{:.1f}', 'RdYlGn', '%', True),
        ('Num_Contigs', 'Total Contigs', 'Total number of contigs in assembly', '{:,}', 'Blues', '', True),
        ('Plasmid_Contigs', 'Plasmid Contigs', 'Number of plasmid contigs', '{:,}', 'BuPu', '', True),
        ('Chromosomal_Contigs', 'Chromosomal Contigs', 'Number of chromosomal contigs', '{:,}', 'Blues', '', True),
        ('Closest_Genome_Ref', 'Closest Ref', 'Closest genome reference - click to view on NCBI', '{}', False, '', False),
        ('Closest_Genome_ANI', 'ANI', 'Average Nucleotide Identity to closest reference', '{:.2f}', 'RdYlGn', '%', False),
        ('Closest_Genome_AF', 'AF', 'Alignment Fraction to closest reference', '{:.2f}', 'RdYlGn', '', True),
        ('Virulence_Factor_Count', 'Virulence Factors', 'Number of unique virulence factors detected', '{:,}', 'Purples', '', True),
        ('Virulence_Factor_Display', 'Virulence Factor List', 'Names of virulence factors detected', '{}', False, '', True),
        ('ASM_Type', 'ASM Type', 'Assembly type', '{}', False, '', True),
    ]

    pconfig_list = []
    for col_id, col_title, desc, fmt, scale, suffix, hidden in general_stats_configs:
        entry = {'title': col_title, 'description': desc, 'format': fmt}
        entry['scale'] = scale if scale else False
        if suffix:
            entry['suffix'] = suffix
        if hidden:
            entry['hidden'] = True
        pconfig_list.append({col_id: entry})

    general_stats_data = {}
    for sample_name in sorted(samples.keys()):
        data = dict(samples[sample_name])
        asm_type = extract_asm_type(sample_name)
        # Skip read-only samples (e.g. _L, _S when only hybrid HL/HS was assembled).
        has_assembly_data = any(k in data for k in ('Isolate_Identity', 'MLST_ST', 'N50_kbp', 'Total_Length_Mbp'))
        if not has_assembly_data:
            continue
        if 'Closest_Genome_Ref' in data:
            data['Closest_Genome_Ref'] = ncbi_genome_link(data['Closest_Genome_Ref'])
        data['ASM_Type'] = asm_type
        # Scope 'plus' hits are excluded from the report.
        data.pop('Plus_AMR_Genes', None)
        data.pop('Plus_AMR_Gene_List', None)
        general_stats_data[sample_name] = data

    write_yaml('general_stats_mqc.yaml', {
        'id': 'custom_general_stats',
        'section_name': 'Isolate Summary and QC Metrics',
        'description': 'Combined statistics from all pipeline modules',
        'plot_type': 'generalstats',
        'pconfig': pconfig_list,
        'data': general_stats_data
    })

    # ── Species lookup for adding Species column to tables ──
    species_lookup = {}
    import re as _re_sp
    for sample_name, row in gtdbtk_data.items():
        sp = row.get('Species', '')
        if sp and sp != '-':
            # Remove GTDB polyphyletic suffixes (e.g. _A, _B, _C, _D)
            sp = _re_sp.sub(r'_[A-Z]$', '', sp)
            sp = _re_sp.sub(r'_[A-Z] ', ' ', sp)
            species_lookup[sample_name] = sp

    # ── MLST Table YAML ──
    mlst_headers = {
        'Species': {'title': 'Species', 'description': species_source, 'scale': False},
        'Scheme': {'title': 'Scheme', 'description': 'MLST scheme detected', 'scale': False},
        'ST': {'title': 'ST', 'description': 'Sequence Type assigned', 'scale': False},
    }
    for i in range(1, 8):
        mlst_headers[f'Locus_{i}'] = {'title': f'Locus {i}', 'description': f'Allele at locus {i}', 'scale': False}
    mlst_headers['ASM_Type'] = {'title': 'ASM Type', 'description': 'Assembly type', 'scale': False}

    # Add ASM_Type and Species to mlst_data
    for sample_name in mlst_data:
        mlst_data[sample_name]['ASM_Type'] = extract_asm_type(sample_name)
        mlst_data[sample_name]['Species'] = species_lookup.get(sample_name, '-')

    write_yaml('mlst_table_mqc.yaml', {
        'id': 'mlst_results_table',
        'section_name': 'Multilocus Sequence Typing',
        'description': 'Multi-Locus Sequence Typing results showing scheme, sequence type, and individual allele values',
        'plot_type': 'table',
        'pconfig': {'id': 'mlst_results_table', 'title': 'MLST Typing Results', 'sortRows': True, 'col1_header': 'Sample'},
        'headers': mlst_headers,
        'data': dict(sorted(mlst_data.items()))
    })

    # ── Taxonomy Table YAML ──
    taxonomy_col_meta = {
        'Domain': 'Taxonomic domain', 'Phylum': 'Taxonomic phylum', 'Class': 'Taxonomic class',
        'Order': 'Taxonomic order', 'Family': 'Taxonomic family', 'Genus': 'Taxonomic genus',
        'Species': 'Taxonomic species', 'Closest_Reference': 'Closest genome reference',
        'ANI': 'Average Nucleotide Identity (%); containment ANI for sylph',
        'AF': 'Alignment Fraction (GTDB-Tk only; sylph does not compute one)'
    }
    taxonomy_headers = {}
    for col, desc in taxonomy_col_meta.items():
        entry = {'title': col.replace('_', ' '), 'description': desc, 'scale': False}
        if col == 'ANI':
            entry['suffix'] = '%'
            entry['scale'] = 'RdYlGn'
            entry['format'] = '{:.2f}'
        elif col == 'AF':
            entry['scale'] = 'RdYlGn'
            entry['format'] = '{:.2f}'
        taxonomy_headers[col] = entry
    taxonomy_headers['ASM_Type'] = {'title': 'ASM Type', 'description': 'Assembly type', 'scale': False}
    taxonomy_headers['Method'] = {'title': 'Method', 'description': 'Species engine', 'scale': False}
    if taxonomy_tool == 'sylph':
        taxonomy_headers['Call'] = {'title': 'Call', 'description':
            'species: ANI >= 95; borderline: 94-95 (species kept, review); genus: species left empty; '
            'ambiguous: another species within the margin; unclassified: no hit', 'scale': False}
        taxonomy_headers['Second_Species'] = {'title': 'Next Species', 'description': 'Best hit of a different species', 'scale': False}
        taxonomy_headers['Second_ANI'] = {'title': 'Next ANI', 'description': 'Its containment ANI (%)', 'scale': False}
    if taxonomy_tool == 'both':
        taxonomy_headers['Sylph_Species'] = {'title': 'Sylph Species', 'description': 'sylph call (comparison only)', 'scale': False}
        taxonomy_headers['Sylph_ANI'] = {'title': 'Sylph ANI', 'description': 'sylph containment ANI (%)', 'scale': False}
        taxonomy_headers['Agreement'] = {'title': 'Agreement', 'description': 'GTDB-Tk vs sylph species', 'scale': False}

    taxonomy_table_data = {}
    for sample_name in sorted(gtdbtk_data.keys()):
        row = dict(gtdbtk_data[sample_name])
        if 'Closest_Reference' in row:
            row['Closest_Reference'] = ncbi_genome_link(row['Closest_Reference'])
        # Convert ANI/AF to numeric for coloring
        for field in ['ANI', 'AF']:
            if field in row and row[field] != '-':
                try:
                    row[field] = float(row[field])
                except (ValueError, TypeError):
                    pass
        row['ASM_Type'] = extract_asm_type(sample_name)
        if taxonomy_tool != 'sylph':
            for k in ('Call', 'Second_Species', 'Second_ANI'):
                row.pop(k, None)
        if taxonomy_tool == 'both':
            sy = sylph_compare.get(sample_name, {})
            row['Sylph_Species'] = sy.get('species') or sy.get('genus') or '-'
            row['Sylph_ANI'] = sy.get('ani') or '-'
            row['Agreement'] = taxonomy_agreement(sample_name)
        taxonomy_table_data[sample_name] = row

    write_yaml('taxonomy_table_mqc.yaml', {
        'id': 'taxonomy_table',
        'section_name': 'Taxonomic Classification',
        'description': 'Taxonomic classification with closest genome reference and ANI values',
        'plot_type': 'table',
        'pconfig': {'id': 'taxonomy_table', 'title': 'Taxonomic Classification', 'sortRows': True, 'col1_header': 'Sample'},
        'headers': taxonomy_headers,
        'data': taxonomy_table_data
    })

    # ── Plasmid Table YAML ──
    mob_headers_def = [
        ('Num_Plasmid_Clusters', 'Plasmid Clusters', 'Number of plasmid clusters detected', 'BuPu'),
        ('Num_Plasmid_Contigs', 'Plasmid Contigs', 'Number of plasmid contigs', 'BuPu'),
        ('Num_Chromosomal_Contigs', 'Chromosomal Contigs', 'Number of chromosomal contigs', 'Blues'),
        ('Total_Plasmid_Size_K', 'Total Size (Kb)', 'Total plasmid size in kilobases', 'Greens'),
        ('Replicon_Types', 'Replicon Types', 'Plasmid replicon types identified', False),
        ('Relaxase_Types', 'Relaxase Types', 'Relaxase enzyme types detected', False),
        ('Mash_Nearest', 'Closest Match', 'Nearest plasmid match', False),
        ('Primary_Cluster_IDs', 'Primary Cluster IDs', 'Primary cluster identifiers', False),
        ('Secondary_Cluster_IDs', 'Secondary Cluster IDs', 'Secondary cluster identifiers', False),
    ]
    mob_headers = {}
    for col_id, col_title, desc, scale in mob_headers_def:
        entry = {'title': col_title, 'description': desc, 'scale': scale if scale else False}
        if col_id == 'Total_Plasmid_Size_K':
            entry['format'] = '{:,.1f}'
            entry['suffix'] = ' Kb'
        mob_headers[col_id] = entry
    mob_headers['ASM_Type'] = {'title': 'ASM Type', 'description': 'Assembly type', 'scale': False}

    # Hyperlink mash_nearest in mob_data
    for sample_name in mob_data:
        mash_val = mob_data[sample_name].get('Mash_Nearest', '-')
        if mash_val and mash_val != '-':
            parts = [v.strip() for v in mash_val.split(';')]
            mob_data[sample_name]['Mash_Nearest'] = '; '.join(ncbi_nuccore_link(p) for p in parts if p)
        mob_data[sample_name]['ASM_Type'] = extract_asm_type(sample_name)

    # Plasmid ID table removed from report: data available in General Statistics
    write_yaml('plasmid_table_mqc.yaml', {})

    # ── Plasmid Resistance Markers YAML ──
    plasmid_res_headers = {
        'Sample': {'title': 'Sample', 'description': 'Sample name', 'scale': False},
        'Species': {'title': 'Species', 'description': species_source, 'scale': False},
        'Cluster_ID': {'title': 'Cluster', 'description': 'Plasmid cluster identifier', 'scale': False},
        'Secondary_Cluster_ID': {'title': 'Secondary Cluster', 'description': 'Secondary cluster identifier', 'scale': False},
        'Num_Contigs': {'title': 'Num Contigs', 'description': 'Number of contigs in this plasmid', 'scale': 'Blues'},
        'Size_K': {'title': 'Size (Kb)', 'description': 'Total plasmid size in kilobases', 'scale': 'Greens', 'format': '{:,.1f}', 'suffix': ' Kb'},
        'Circularity': {'title': 'Circularity', 'description': 'Circular if any contig in cluster is circular', 'scale': False},
        'Num_Resistance_Genes': {'title': 'AMR Genes', 'description': 'Number of resistance genes', 'scale': 'OrRd'},
        'Resistance_Genes': {'title': 'Resistance Genes', 'description': 'Resistance genes found on this plasmid', 'scale': False},
        'Rep_Type': {'title': 'Rep Type', 'description': 'Replicon type', 'scale': False},
        'Relaxase_Type': {'title': 'Relaxase Type', 'description': 'Relaxase enzyme type', 'scale': False},
        'MPF_Type': {'title': 'MPF Type', 'description': 'Mating pair formation type', 'scale': False},
        'OriT_Type': {'title': 'OriT Type', 'description': 'Origin of transfer type', 'scale': False},
        'Predicted_Mobility': {'title': 'Predicted Mobility', 'description': 'Predicted plasmid mobility', 'scale': False},
        'Closest_Match': {'title': 'Closest Match', 'description': 'Nearest plasmid match', 'scale': False},
        'Contigs': {'title': 'Contigs', 'description': 'Contig IDs in this plasmid cluster', 'scale': False},
        'ASM_Type': {'title': 'ASM Type', 'description': 'Assembly type', 'scale': False},
    }

    # Build color mapping for cluster IDs (same cluster = same color)
    _cluster_palette = CATEGORICAL_PALETTE
    _all_cluster_vals = set()
    for _row in plasmid_markers.values():
        _cv = _row.get('Cluster_ID', '-')
        if _cv and _cv != '-':
            _all_cluster_vals.add(_cv)
        _sv = _row.get('Secondary_Cluster_ID', '-')
        if _sv and _sv != '-':
            _all_cluster_vals.add(_sv)
    _cluster_color_map = {}
    for _i, _cv in enumerate(sorted(_all_cluster_vals)):
        _cluster_color_map[_cv] = _cluster_palette[_i % len(_cluster_palette)]

    # Build cond_formatting for Cluster_ID and Secondary_Cluster_ID
    _cluster_fmt_rules = {}
    _cluster_fmt_colours = []
    for _cv, _col in _cluster_color_map.items():
        _rule_name = f'cluster_{_cv}'
        _cluster_fmt_rules[_rule_name] = [{'s_eq': _cv}]
        _cluster_fmt_colours.append({_rule_name: _col})

    if _cluster_fmt_rules:
        plasmid_res_headers['Cluster_ID']['cond_formatting_rules'] = dict(_cluster_fmt_rules)
        plasmid_res_headers['Cluster_ID']['cond_formatting_colours'] = list(_cluster_fmt_colours)
        plasmid_res_headers['Secondary_Cluster_ID']['cond_formatting_rules'] = dict(_cluster_fmt_rules)
        plasmid_res_headers['Secondary_Cluster_ID']['cond_formatting_colours'] = list(_cluster_fmt_colours)

    # Add ASM_Type and Species to plasmid_markers
    for row_key in plasmid_markers:
        sample_part = row_key.split(':')[0] if ':' in row_key else row_key
        plasmid_markers[row_key]['ASM_Type'] = extract_asm_type(sample_part)
        plasmid_markers[row_key]['Species'] = species_lookup.get(sample_part, '-')

    # Row order: cluster, then secondary cluster, then resistance genes, then the plasmid row key.
    def _plasmid_sort_key(item):
        row_key, row = item
        def rank(value):
            v = str(value or '').strip()
            return (1, '') if v in ('', '-', 'unknown') else (0, v)
        return (rank(row.get('Cluster_ID')),
                rank(row.get('Secondary_Cluster_ID')),
                rank(row.get('Resistance_Genes')),
                row_key)

    write_yaml('plasmid_resistance_mqc.yaml', {
        'id': 'plasmid_resistance_table',
        'section_name': 'Plasmid Clustering and Resistance Markers',
        'description': 'Resistance genes located on plasmid contigs, grouped by plasmid cluster. Note: If any contig within a cluster is circular, the entire cluster is marked as circular.',
        'plot_type': 'table',
        'pconfig': {'id': 'plasmid_resistance_table', 'title': 'Plasmid Clustering and Resistance Markers', 'sortRows': False, 'col1_header': 'Plasmid', 'no_violin': True},
        'headers': plasmid_res_headers,
        'data': dict(sorted(plasmid_markers.items(), key=_plasmid_sort_key))
    })

    # ── Chromosomal Resistance Markers ──
    # Section suppressed: the same information is in "Chromosomal AMR Gene Profile by Sample" in a
    # form that compares across samples.
    open('chromosome_resistance_mqc.yaml', 'w').close()

    # ── Antimicrobial Resistance Profile by Sample ──
    # Rows are antimicrobials (with their drug class), columns are samples.
    resistance_profile = defaultdict(dict)   # sample -> {antimicrobial: 'R'|'S'}
    resistance_genes = defaultdict(dict)     # sample -> {antimicrobial: 'geneA; geneB'}
    antimicrobial_class = {}                 # antimicrobial -> class (title-cased)
    all_antimicrobials = set()

    for sample_name, rows in resfinder_pheno.items():
        for row in rows:
            if not resfinder_reportable(row):
                continue
            drug = title_case(row['Antimicrobial'].strip())
            phenotype = row.get('WGS-predicted phenotype', '').strip().lower()
            genes = parse_genetic_background(row.get('Genetic background', ''))
            all_antimicrobials.add(drug)
            antimicrobial_class.setdefault(drug, title_case(row.get('Class', '').strip()) or '-')
            # A "Resistant" call with an empty Genetic background comes from PointFinder, not from an
            # acquired gene: the acquired-gene table has nothing to attribute it to, and the ResFinder web
            # tool does not report it as resistant here either.
            is_resistant = phenotype == 'resistant' and bool(genes)
            resistance_profile[sample_name][drug] = 'R' if is_resistant else 'S'
            resistance_genes[sample_name][drug] = genes

    # Only samples with an assembly appear; read-only rows would be blank columns.
    assembled_samples = [s for s in samples if any(k in samples[s] for k in ('Isolate_Identity', 'MLST_ST', 'N50_kbp', 'Total_Length_Mbp'))]
    profile_samples = sorted(set(list(resistance_profile.keys())) & set(assembled_samples)) or sorted(resistance_profile.keys())

    # Row order: every antimicrobial with at least one resistant sample first, then the rest;
    # within each block sort by class, then by antimicrobial.
    def row_sort_key(am):
        any_resistant = any(resistance_profile.get(s, {}).get(am) == 'R' for s in profile_samples)
        return (0 if any_resistant else 1, antimicrobial_class.get(am, '-'), am)

    ordered_antimicrobials = sorted(all_antimicrobials, key=row_sort_key)

    # One background colour per class so a class reads as a contiguous block.
    class_rules = {}
    class_colours = []
    for i, cls in enumerate(sorted({antimicrobial_class.get(am, '-') for am in ordered_antimicrobials})):
        rule = f'amr_class_{i}'
        class_rules[rule] = [{'s_eq': cls}]
        class_colours.append({rule: CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)]})

    resistance_profile_headers = {
        'Class': {
            'title': 'Class',
            'description': 'Antimicrobial class',
            'scale': False,
            'cond_formatting_rules': class_rules,
            'cond_formatting_colours': class_colours,
        },
    }
    for sample_name in profile_samples:
        # No cond_formatting here: the cell carries the badge markup itself so the causal genes can
        # sit outside the badge.
        resistance_profile_headers[sample_name] = {
            'title': sample_name,
            'description': f'{sample_name} predicted phenotype (R=Resistant, S=Susceptible, -=not assessed), with the genes responsible',
            'scale': False,
        }

    # MultiQC's own badge colours for the reserved pass/warn rules.
    PHENOTYPE_BADGE = {'R': '#f0ad4e', 'S': '#5cb85c'}

    def phenotype_cell(call, genes):
        """Badge for the R/S call, with the causal genes as plain text beside it."""
        if call not in PHENOTYPE_BADGE:
            return '-'
        badge = f'<span class="badge" style="background-color:{PHENOTYPE_BADGE[call]}">{call}</span>'
        if call == 'R' and genes:
            return f'{badge} <span class="amr-genes">{genes}</span>'
        return badge

    resistance_profile_data = {}
    for am in ordered_antimicrobials:
        row_data = {'Class': antimicrobial_class.get(am, '-')}
        for sample_name in profile_samples:
            call = resistance_profile.get(sample_name, {}).get(am, '-')
            row_data[sample_name] = phenotype_cell(
                call, resistance_genes.get(sample_name, {}).get(am, ''))
        resistance_profile_data[am] = row_data

    if resistance_profile_data:
        write_yaml('resistance_profile_mqc.yaml', {
            'id': 'resistance_profile_table',
            'section_name': 'Phenotypic AMR Profile based on Acquired Genes by Sample',
            'description': 'Predicted resistance phenotype attributable to acquired resistance genes.',
            'plot_type': 'table',
            'pconfig': {
                'id': 'resistance_profile_table',
                'title': 'Phenotypic AMR Profile based on Acquired Genes by Sample',
                'sortRows': False,
                'col1_header': 'Antimicrobial',
            },
            'headers': resistance_profile_headers,
            'data': resistance_profile_data,
        })
    else:
        # Empty file: MULTIQC skips *_mqc.yaml under 11 bytes, so the section disappears rather than
        # rendering as an empty table.
        open('resistance_profile_mqc.yaml', 'w').close()

    # ── Phenotypic AMR Profile based on Point Mutations by Sample ──
    # Rows are the antimicrobial/class PointFinder attributes the mutation to, columns are
    # samples, cells carry the R/S badge with the causal mutations beside it: same shape as the
    # acquired-gene table above.
    pointfinder_profile = defaultdict(dict)   # sample -> {resistance: [mutations]}
    all_point_resistances = set()

    for sample_name, rows in pointfinder_hits.items():
        for row in rows:
            resistance = title_case((row.get('Resistance') or '').strip())
            mutation = format_point_mutation(row.get('Mutation', ''))
            if not resistance or not mutation:
                continue
            all_point_resistances.add(resistance)
            muts = pointfinder_profile[sample_name].setdefault(resistance, [])
            if mutation not in muts:
                muts.append(mutation)

    point_samples = sorted(set(pointfinder_hits.keys()) & set(assembled_samples)) \
        or sorted(pointfinder_hits.keys())

    if all_point_resistances and point_samples:
        # Resistant rows first, then alphabetical, matching the acquired-gene table.
        ordered_point_resistances = sorted(
            all_point_resistances,
            key=lambda r: (0 if any(pointfinder_profile.get(sm, {}).get(r)
                                    for sm in point_samples) else 1, r))

        point_headers = {}
        for sample_name in point_samples:
            point_headers[sample_name] = {
                'title': sample_name,
                'description': f'{sample_name} predicted phenotype from point mutations (R=Resistant, S=Susceptible)',
                'scale': False,
            }

        point_data = {}
        for resistance in ordered_point_resistances:
            row_data = {}
            for sample_name in point_samples:
                muts = pointfinder_profile.get(sample_name, {}).get(resistance)
                if sample_name not in pointfinder_hits:
                    row_data[sample_name] = '-'
                elif muts:
                    row_data[sample_name] = phenotype_cell('R', '; '.join(muts))
                else:
                    row_data[sample_name] = phenotype_cell('S', '')
            point_data[resistance] = row_data

        write_yaml('pointfinder_profile_mqc.yaml', {
            'id': 'pointfinder_profile_table',
            'section_name': 'Phenotypic AMR Profile based on Point Mutations by Sample',
            'description': 'Predicted resistance phenotype attributable to chromosomal point mutations.',
            'plot_type': 'table',
            'pconfig': {
                'id': 'pointfinder_profile_table',
                'title': 'Phenotypic AMR Profile based on Point Mutations by Sample',
                'sortRows': False,
                'col1_header': 'Antimicrobial/Class',
            },
            'headers': point_headers,
            'data': point_data,
        })
    else:
        # No PointFinder-supported species in this run.
        open('pointfinder_profile_mqc.yaml', 'w').close()

    # ── Presence/absence profiles ──
    # Three tables built the same way: samples as rows, features as columns, with the species
    # prepended for context.
    pointfinder_all = {s: [format_point_mutation(r.get('Mutation', ''))
                           for r in pointfinder_hits.get(s, [])]
                       for s in samples}

    # Both AMR callers report a contig, so both are split on the MOB-recon assignment.
    _gene_profile_samples = sorted(set(list(sample_contig_map.keys())
                                       + list(sample_amr_genes.keys())
                                       + list(resfinder_by_contig.keys())))
    assembled_samples = [s for s in samples if any(k in samples[s] for k in ('Isolate_Identity', 'MLST_ST', 'N50_kbp', 'Total_Length_Mbp'))]

    def build_feature_profile(hits_by_sample, molecule_type=None):
        """sample -> {feature: '+'}, optionally restricted to one molecule type."""
        profile = defaultdict(dict)
        features = set()
        for sample_name in _gene_profile_samples:
            contig_map = sample_contig_map.get(sample_name, {})
            hits = hits_by_sample.get(sample_name, {})
            for cid, found in hits.items():
                if molecule_type is not None:
                    info = contig_map.get(cid)
                    if not info or info['molecule_type'] != molecule_type:
                        continue
                for feature in found:
                    profile[sample_name][feature] = '+'
                    features.add(feature)
        return profile, sorted(features)

    def write_feature_profile(profile, features, yaml_name, section_id, section_name,
                              description, col_description):
        """Emit a presence/absence table, or nothing when no feature was found."""
        if not features:
            return
        rows = sorted(set(list(profile.keys()) + assembled_samples))
        headers = {
            'Species': {'title': 'Species', 'description': species_source, 'scale': False},
        }
        for feature in features:
            headers[feature] = {
                'title': feature,
                'description': f'{feature} ({col_description})',
                'scale': False,
                'cond_formatting_rules': {
                    'pass': [{'s_eq': '-'}],
                    'warn': [{'s_eq': '+'}],
                },
                'cond_formatting_colours': [
                    {'pass': '#dcffe0'},
                    {'warn': '#ffd2d2'},
                ],
            }

        data = {}
        for sample_name in rows:
            row = {f: profile.get(sample_name, {}).get(f, '-') for f in features}
            row['Species'] = species_lookup.get(sample_name, '-')
            data[sample_name] = row

        write_yaml(yaml_name, {
            'id': section_id,
            'section_name': section_name,
            'description': description,
            'plot_type': 'table',
            'pconfig': {'id': section_id, 'title': section_name, 'sortRows': True, 'col1_header': 'Sample'},
            'headers': headers,
            'data': dict(sorted(data.items())),
        })

    amrfinder_all = {s: [g for genes in hits.values() for g in genes]
                     for s, hits in sample_amr_genes.items()}

    # Both callers report a contig, so both are placed on the MOB-recon assignment.
    localised_hits = defaultdict(lambda: defaultdict(list))
    for sample_name, hits in sample_amr_genes.items():
        for contig, genes in hits.items():
            localised_hits[sample_name][contig].extend(genes)
    for sample_name, hits in resfinder_by_contig.items():
        known = {g.lower() for g in amrfinder_all.get(sample_name, [])}
        for contig, genes in hits.items():
            for gene in genes:
                if gene.lower() in known:
                    continue
                if gene not in localised_hits[sample_name][contig]:
                    localised_hits[sample_name][contig].append(gene)

    chromo_profile, chromo_features = build_feature_profile(localised_hits, 'chromosome')
    write_feature_profile(
        chromo_profile, chromo_features,
        'chromo_gene_profile_mqc.yaml', 'chromo_gene_profile_table',
        'Chromosomal AMR Gene Profile by Sample',
        'Acquired AMR genes located on chromosomal contigs. '
        '<b>+</b> = Detected, <b>-</b> = Not detected.',
        '+ = Detected, - = Not detected')

    plasmid_profile, plasmid_features = build_feature_profile(localised_hits, 'plasmid')
    write_feature_profile(
        plasmid_profile, plasmid_features,
        'plasmid_gene_profile_mqc.yaml', 'plasmid_gene_profile_table',
        'Plasmid AMR Gene Profile by Sample',
        'Acquired AMR genes located on plasmid contigs. '
        '<b>+</b> = Detected, <b>-</b> = Not detected.',
        '+ = Detected, - = Not detected')

    point_profile, point_features = build_feature_profile(sample_point_mutations)
    # Point mutations are chromosomal by definition, so PointFinder's calls need no location
    # caveat and merge in as plain '+'.
    point_features = set(point_features)
    for sample_name, muts in pointfinder_all.items():
        for mutation in muts:
            if not mutation:
                continue
            existing = {f.lower() for f in point_profile[sample_name]}
            if mutation.lower() not in existing:
                point_profile[sample_name][mutation] = '+'
            point_features.add(mutation)
    point_features = sorted(point_features)

    write_feature_profile(
        point_profile, point_features,
        'point_mutation_profile_mqc.yaml', 'point_mutation_profile_table',
        'AMR Point Mutation Profile by Sample',
        'Resistance-conferring point mutations. '
        '<b>+</b> = Detected, <b>-</b> = Not detected.',
        '+ = Detected, - = Not detected')

    # ── Software Versions (written to file, not in report) ──
    versions = {
        'SeqKit': '2.12.0', 'NanoPlot': '1.46.2', 'fastp': '1.0.1',
        'fastplong': '0.4.1', 'Rasusa': '3.0.0', 'Flye': '2.9.6', 'Medaka': '2.0.0',
        'Unicycler': '0.5.1', 'dnaapler': '1.3.0', 'BWA-MEM2': '2.3',
        'Polypolish': '0.6.1', 'minimap2': '2.30', 'samtools': '1.21',
        'mosdepth': '0.3.3', 'QUAST': '5.3.0', 'BUSCO': '6.0.0',
        'MLST': '2.32.2',
        'AMRFinderPlus': '4.2.7-2026-08-07.1', 'ResFinder': '4.7.2',
        'MOB-suite': '3.1.9', 'MultiQC': '1.33'
    }
    if taxonomy_tool in ('gtdbtk', 'both'):
        versions['GTDB-Tk'] = '2.6.1'
    if taxonomy_tool in ('sylph', 'both'):
        versions['sylph'] = '1.0.0'
    write_versions('software_versions.tsv', versions)

    # ── Filtered Reads Bar Plot ──
    filtered_reads_plot_data = {}
    for sample_name in sorted(samples.keys()):
        before = samples[sample_name].get('Reads_Before_Filter', 0)
        after = samples[sample_name].get('Reads_After_Filter', 0)
        failed = before - after
        if before > 0:
            filtered_reads_plot_data[sample_name] = {
                'Reads Passed': after,
                'Reads Failed': failed,
            }

    write_yaml('filtered_reads_plot_mqc.yaml', {
        'id': 'filtered_reads_plot',
        'section_name': 'Read Filtering Summary',
        'description': 'Number of reads before and after quality filtering (fastp/fastplong)',
        'plot_type': 'bargraph',
        'pconfig': {
            'id': 'filtered_reads_bargraph',
            'title': 'Read Filtering Summary',
            'ylab': 'Number of Reads',
            'cpswitch_counts_label': 'Read Count',
            'cpswitch_c_active': False,
        },
        'data': filtered_reads_plot_data
    })

    # ── Cumulative Coverage Distribution Plot (from mosdepth global.dist.txt) ──
    coverage_plot_data = {}
    if os.path.exists('mosdepth_dist'):
        for f in os.listdir('mosdepth_dist'):
            if f.endswith('.mosdepth.global.dist.txt'):
                sample_name = f.replace('.mosdepth.global.dist.txt', '')
                coverage_plot_data[sample_name] = {}
                with open(f'mosdepth_dist/{f}', 'r') as fh:
                    for line in fh:
                        parts = line.strip().split('\t')
                        if len(parts) >= 3 and parts[0] == 'total':
                            depth = int(parts[1])
                            fraction = float(parts[2])
                            if depth <= 500:
                                coverage_plot_data[sample_name][depth] = fraction

    write_yaml('coverage_plot_mqc.yaml', {
        'id': 'coverage_depth_plot',
        'section_name': 'Cumulative Coverage Distribution',
        'description': 'Proportion of bases covered at a given depth or higher (mosdepth)',
        'plot_type': 'linegraph',
        'pconfig': {
            'id': 'cumulative_coverage_linegraph',
            'title': 'Cumulative Coverage Distribution',
            'xlab': 'Coverage Depth (x)',
            'ylab': 'Fraction of Bases',
            'ymin': 0,
            'ymax': 1,
            'xmin': 0,
            'xmax': 500,
            'tt_label': '{point.x}x: {point.y:.3f}',
            'smooth_points': 500,
        },
        'data': coverage_plot_data
    })

    # ── MultiQC config YAML ──
    report_header = [
        {'Workflow ID': run_id},
    ]
    if args.analysis_name:
        report_header.append({'Workflow Description': args.analysis_name})
    report_header.append({'Group ID': args.group_id})
    report_header.append({'Group Description': args.group_description})

     # ── AI summary prompts (used only when params.ai_summary / ai_summary_full is on) ──
    AI_PROMPT_FULL = (
        "You are an expert bacterial genomic-surveillance analyst writing for public-health and\n"
        "hospital infection-prevention staff. You are interpreting an EPIquence report that combines\n"
        "QC, assembly, species identification and typing, antimicrobial resistance, plasmid analysis,\n"
        "and sample-to-sample clustering.\n\n"
        "OUTPUT FORMAT - follow exactly:\n"
        "- Do NOT output an overall report title.\n"
        '- Open with a "### Key Findings" section: 3 to 5 bullets giving the most important run-level takeaways (highest-consequence resistance, any outbreak or transmission signal, notable species, overall data confidence). Then continue with the subsections below.\n'
        '- Use "### " markdown headings for every section.\n'
        '- Under every heading write proper markdown bullets: each bullet starts on its own new line with "- ". Never put a bullet marker mid-paragraph.\n'
        "- SUMMARISE within each subsection. Do NOT produce one bullet per sample; describe patterns across the run and always indicate the sample names in the patterns\n"
        "- Use **bold** for any sample ID you mention (**SAMPLE-04_HL**), species (**Klebsiella pneumoniae**), sequence type (**ST258**), serotype/serovar, resistance class (**carbapenem-resistant**), drug and gene/allele (**bla_KPC-3**, **mgrB**), plasmid replicon/cluster, cluster ID, and QC verdict (**PASS**, **WARN**, **FAIL**).\n"
        '- PLAIN TEXT ONLY - no LaTeX and no math delimiters. Write ">=" and "<=" as plain ASCII, Greek letters as words, and "%" as a bare percent sign.\n\n'
        "Interpret; do not restate table values. Give the biological, clinical, or epidemiological MEANING\n"
        "and the REASONING. Sample suffixes indicate assembly route: _L = ONT long-read, _S = Illumina\n"
        "short-read, _HL = hybrid (long-first), _HS = hybrid (short-first), _A = pre-assembled. Coverage\n"
        "depth and assembly quality (QUAST, BUSCO completeness) gate confidence in downstream calls.\n\n"
        "### Key Findings\n"
        "- The 3 to 5 most important run-level takeaways, most urgent first.\n\n"
        "### Species and Typing\n"
        "- Summarise the species (GTDB taxonomy) and MLST or other typing seen across the run; flag unexpected species, mixed or contaminated calls, or typing that disagrees with the taxonomy.\n\n"
        "### Antimicrobial Resistance\n"
        "- Summarise resistance by drug class, naming causal genes/mutations and whether they are chromosomal or plasmid-borne; explain treatment impact and call out high-consequence phenotypes (carbapenemase, ESBL, colistin); Indicate samples of concern.\n\n"
        "### Plasmids and Mobile Elements\n"
        "- Summarise notable plasmids (replicon, size, predicted mobility) and the resistance or virulence they carry; note plasmids shared across samples as possible horizontal transfer; Indicate samples of concern.\n\n"
        "### Clustering and Outbreak Signals\n"
        "- Interpret the outbreak (clonal) and plasmid-HGT groupings: which samples group together and whether that fits transmission or a shared plasmid rather than clonal spread; state which links are strong versus weak given coverage and assembly.\n\n"
        "### Quality and Confidence\n"
        "- Say which results are trustworthy and which are not, connecting taxonomic purity, coverage depth, and assembly completeness. Do NOT call a sample well covered when depth is high but breadth or genome fraction is low.\n\n"
        "### Recommended Next Steps\n"
        "- Concrete, practical actions for the surveillance team, one bullet each."
    )
    AI_PROMPT_SHORT = (
        "Write a concise interpretive assessment of this EPIquence bacterial genomic-surveillance report.\n"
        "- Do NOT output a title.\n"
        '- Open with a "### Key Findings" section of 2 to 4 bullets giving the most important run-level takeaways, then brief subsections.\n'
        '- Use "### " headings and markdown bullets ("- " at the start of each new line).\n'
        "- SUMMARISE per subsection; do NOT need to go sample by sample. Name individual samples when its needed.\n"
        "- Use **bold** for any sample ID, species, sequence type, resistance class, drug, gene/mutation, plasmid, cluster ID, and QC verdict you mention.\n"
        '- Plain text only: no LaTeX, no $...$ expressions; write ">=" as plain ASCII.\n'
        "Cover: species and typing; antimicrobial resistance with causal genes and chromosomal-versus-plasmid location and treatment impact; plasmids shared across samples; outbreak and plasmid-HGT clusters and whether they suggest transmission or a shared plasmid; and confidence given coverage, taxonomic purity, and assembly completeness. "
        "Sample suffixes: _L ONT, _S Illumina, _HL/_HS hybrid, _A pre-assembled. Interpret, do not restate table values."
    )

    write_yaml('multiqc_config.yaml', {
        'title': args.title,
        'subtitle': args.subtitle,
        'intro_text': args.intro,
        'thousandsSep_format': ',',
        'ai_prompt_short': AI_PROMPT_SHORT,
        'ai_prompt_full': AI_PROMPT_FULL,
        'show_analysis_paths': False,
        'show_analysis_time': True,
        'show_intro': False,
        'report_header_info': report_header,
        'module_order': ['custom_content', 'quast', 'busco'],
        'run_modules': ['custom_content', 'quast', 'busco'],
        'skip_versions_section': True,
        'custom_content': {
            'order': [
                'custom_general_stats', 'mlst_results_table', 'taxonomy_table',
                'plasmid_resistance_table', 'chromosome_resistance_table',
                'resistance_profile_table',
                'chromo_gene_profile_table', 'plasmid_gene_profile_table',
                'point_mutation_profile_table', 'pointfinder_profile_table',
                'outbreak_detection',
                'outbreak_cluster_summary',
                'plasmid_hgt_detection',
                'plasmid_hgt_2_plots_HL', 'plasmid_hgt_summary_HL', 'plasmid_hgt_1_table_HL',
                'plasmid_hgt_2_plots_HS', 'plasmid_hgt_summary_HS', 'plasmid_hgt_1_table_HS',
                'plasmid_hgt_2_plots_L', 'plasmid_hgt_summary_L', 'plasmid_hgt_1_table_L',
                'plasmid_hgt_2_plots_S', 'plasmid_hgt_summary_S', 'plasmid_hgt_1_table_S',
                'filtered_reads_plot', 'coverage_depth_plot'
            ]
        },
        'report_section_order': {
            'ectyper_typing': {'order': -100},
            'kleborate_kpsc_typing': {'order': -101},
            'kleborate_kosc_typing': {'order': -102},
            'emmtyper_typing': {'order': -103},
            'resistance_profile_table': {'order': -200},
            'chromo_gene_profile_table': {'order': -201},
            'plasmid_gene_profile_table': {'order': -202},
            'point_mutation_profile_table': {'order': -203},
            # Directly after the acquired-gene profile (-200): the two phenotype tables are read together.
            'pointfinder_profile_table': {'order': -199},
            'outbreak_detection': {'order': -500},
            'outbreak_cluster_summary': {'order': -501},
            'plasmid_hgt_detection': {'order': -505},
            # NOTE: sections render in ASCENDING order value (lowest first).
            'plasmid_hgt_summary_HL': {'order': -517},
            'plasmid_hgt_2_plots_HL': {'order': -516},
            'plasmid_hgt_1_table_HL': {'order': -515},
            'plasmid_hgt_summary_HS': {'order': -514},
            'plasmid_hgt_2_plots_HS': {'order': -513},
            'plasmid_hgt_1_table_HS': {'order': -512},
            'plasmid_hgt_summary_L': {'order': -511},
            'plasmid_hgt_2_plots_L': {'order': -510},
            'plasmid_hgt_1_table_L': {'order': -509},
            'plasmid_hgt_summary_S': {'order': -508},
            'plasmid_hgt_2_plots_S': {'order': -507},
            'plasmid_hgt_1_table_S': {'order': -506},
            'quast': {'order': -600},
            'busco': {'order': -700},
            'filtered_reads_plot': {'order': -800},
            'coverage_depth_plot': {'order': -900},
        },
        # Native-module columns that duplicate the custom General Statistics table (QUAST N50/length)
        # or are better read in the BUSCO section itself.
        'table_columns_visible': {
            'QUAST': {'N50': False, 'Total_length': False},
            'BUSCO': {
                'complete': False, 'complete_single_copy': False,
                'single_copy': False, 'multi_copy': False,
                'fragmented': False, 'missing': False,
            },
        },
        'table_columns_placement': {'quast': {'L50': 900}},
        'custom_table_header_config': {
            'quast_table': {
                'L50': {'modify': '', 'format': '{:,.0f}', 'suffix': ''}
            }
        },
        'sp': {
            'busco': {'fn': 'short_summary.specific.bacteria_odb12.*.txt'},
        },
        'extra_fn_clean_exts': ['.fastq.gz', '.fq.gz', '.fastq', '.fq', '.fasta', '.fa', '.fna', '_quast', '_busco'],
        'fn_clean_exts': [
            {'type': 'remove', 'pattern': 'short_summary.specific.bacteria_odb12.'},
            '.txt'
        ],
        'custom_css_files': ['custom_styles.css'],
        'theme': 'default'
    })

    # ── Combined Summary TSVs (for database ingestion) ──

    def split_sample(sample_name):
        """Split sample like SAMPLE-04_HL into (SAMPLE-04, HL)."""
        parts = sample_name.rsplit('_', 1)
        if len(parts) == 2 and parts[1] in ('L', 'S', 'HL', 'HS'):
            return parts[0], parts[1]
        return sample_name, ''

    # Which resistance genes sit on which plasmid cluster, per sample.
    _cluster_amr = {}
    for _rk, _row in plasmid_markers.items():
        _sample = _rk.rsplit(':', 1)[0] if ':' in _rk else _row.get('Sample', '')
        _cid = str(_row.get('Cluster_ID', '') or '').strip()
        _genes_raw = str(_row.get('Resistance_Genes', '') or '').strip()
        if not _sample or not _cid or _cid == '-' or not _genes_raw or _genes_raw == '-':
            # Clusters with no resistance genes are skipped: a sample with eight cryptic plasmids would
            # otherwise produce a column of "AA155:-".
            continue
        _genes = [g.strip() for g in _genes_raw.split(';') if g.strip() and g.strip() != '-']
        if not _genes:
            continue
        _cluster_amr.setdefault(_sample, []).append((_cid, _genes))

    cluster_amr_map = {}
    for _sample, _entries in _cluster_amr.items():
        # Most-loaded plasmid first, so the clinically interesting one reads first if the cell is
        # truncated; cluster id breaks ties for stable output.
        _entries.sort(key=lambda e: (-len(e[1]), e[0]))
        cluster_amr_map[_sample] = '; '.join(
            f"{_cid}:{','.join(_genes)}" for _cid, _genes in _entries)

    # general_statistics.tsv: raw values, no units/symbols
    gs_columns = [
        ('workflow_id', lambda s, d: run_id),
        ('workflow_description', lambda s, d: workflow_desc),
        ('workflow_timestamp', lambda s, d: workflow_timestamp),
        ('group_id', lambda s, d: group_id),
        ('group_description', lambda s, d: group_desc),
        ('sample', lambda s, d: s),
        ('isolate_id', lambda s, d: split_sample(s)[0]),
        ('assembly_type', lambda s, d: split_sample(s)[1]),
        ('species', lambda s, d: d.get('Isolate_Identity', '')),
        ('mlst_scheme', lambda s, d: d.get('MLST_Scheme', '')),
        ('mlst_st', lambda s, d: d.get('MLST_ST', '')),
        # amr_gene_count / amr_gene_list are the MERGED set.
        ('amr_gene_count', lambda s, d: d.get('Merged_AMR_Genes', '')),
        ('amr_gene_list', lambda s, d: d.get('Merged_AMR_Gene_List', '')),
        ('amrfinder_gene_list', lambda s, d: d.get('AMRFinder_Gene_List', '')),
        ('resfinder_gene_list', lambda s, d: d.get('ResFinder_Gene_List', '')),
        # amr_gene_count / amr_gene_list above are the core set, matching the report and the
        # chromosomal/plasmid splits.
        ('amr_plus_gene_count', lambda s, d: d.get('Plus_AMR_Genes', '')),
        ('amr_plus_gene_list', lambda s, d: d.get('Plus_AMR_Gene_List', '')),
        ('chromosomal_amr_gene_list', lambda s, d: d.get('Chromosomal_AMR_Gene_List', '')),
        ('plasmid_amr_gene_list', lambda s, d: d.get('Plasmid_AMR_Gene_List', '')),
        ('plasmid_cluster_count', lambda s, d: d.get('Plasmid_Clusters', '')),
        ('plasmid_cluster_list', lambda s, d: d.get('Plasmid_Cluster_List', '')),
        ('plasmid_cluster_amr', lambda s, d: cluster_amr_map.get(s, '')),
        ('point_mutation_count', lambda s, d: d.get('Merged_Point_Mutations', '')),
        ('point_mutation_list', lambda s, d: d.get('Merged_Point_Mutation_List', '')),
        ('amrfinder_point_mutation_list', lambda s, d: d.get('Point_Mutation_List', '')),
        ('resistant_antimicrobial_count', lambda s, d: d.get('Resistant_Antimicrobials', '')),
        ('resistant_antimicrobial_list', lambda s, d: d.get('Resistant_Antimicrobial_List', '')),
        ('virulence_factor_count', lambda s, d: d.get('Virulence_Factor_Count', '')),
        ('virulence_factor_list', lambda s, d: d.get('Virulence_Factor_List', '')),
        ('median_phred', lambda s, d: d.get('Median_Phred', '')),
        ('mean_coverage', lambda s, d: d.get('Mean_Coverage', '')),
        ('median_coverage', lambda s, d: d.get('Median_Coverage', '')),
        ('total_reads', lambda s, d: int(d.get('Total_Reads_K', 0) * 1000) if d.get('Total_Reads_K') else ''),
        ('total_bases', lambda s, d: int(d.get('Total_Bases_Mb', 0) * 1000000) if d.get('Total_Bases_Mb') else ''),
        ('avg_read_len', lambda s, d: d.get('Avg_Read_Len', '')),
        ('gc_pct', lambda s, d: d.get('GC_pct', '')),
        ('filtered_reads', lambda s, d: d.get('Reads_After_Filter', int(d.get('Reads_After_Filter_K', 0) * 1000) if d.get('Reads_After_Filter_K') else '')),
        ('pct_reads_passed', lambda s, d: d.get('Pct_Reads_Passed', '')),
        ('filtered_bases', lambda s, d: int(d.get('Bases_After_Filter_Mb', 0) * 1000000) if d.get('Bases_After_Filter_Mb') else ''),
        ('pct_bases_passed', lambda s, d: d.get('Pct_Bases_Passed', '')),
        ('pct_ge_20x', lambda s, d: d.get('Pct_ge_20x', '')),
        ('pct_ge_30x', lambda s, d: d.get('Pct_ge_30x', '')),
        ('pct_ge_40x', lambda s, d: d.get('Pct_ge_40x', '')),
        ('pct_ge_50x', lambda s, d: d.get('Pct_ge_50x', '')),
        ('pct_ge_100x', lambda s, d: d.get('Pct_ge_100x', '')),
        ('total_contigs', lambda s, d: d.get('Num_Contigs', '')),
        ('plasmid_contigs', lambda s, d: d.get('Plasmid_Contigs', '')),
        ('chromosomal_contigs', lambda s, d: d.get('Chromosomal_Contigs', '')),
        ('assembly_length_bp', lambda s, d: int(d.get('Total_Length_Mbp', 0) * 1000000) if d.get('Total_Length_Mbp') else ''),
        ('n50_bp', lambda s, d: int(d.get('N50_kbp', 0) * 1000) if d.get('N50_kbp') else ''),
        ('l50', lambda s, d: d.get('L50', '')),
        ('largest_contig_bp', lambda s, d: int(d.get('Largest_Contig_Mbp', 0) * 1000000) if d.get('Largest_Contig_Mbp') else ''),
        ('closest_genome_ref', lambda s, d: d.get('Closest_Genome_Ref', '')),
        ('ani', lambda s, d: d.get('Closest_Genome_ANI', '')),
        ('af', lambda s, d: d.get('Closest_Genome_AF', '')),
    ]
    with open('general_statistics.tsv', 'w') as f:
        f.write('\t'.join(col[0] for col in gs_columns) + '\n')
        for sample_name in sorted(samples.keys()):
            data = samples[sample_name]
            vals = []
            for col_name, fn in gs_columns:
                try:
                    v = fn(sample_name, data)
                    vals.append(str(v) if v != '' and v is not None else '')
                except:
                    vals.append('')
            f.write('\t'.join(vals) + '\n')

    with open('taxonomy.tsv', 'w') as f:
        tax_cols = ['workflow_id', 'workflow_description', 'workflow_timestamp', 'group_id', 'group_description', 'sample', 'isolate_id', 'assembly_type', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'closest_reference', 'ani', 'af',
                    'method', 'call_status', 'second_species', 'second_ani', 'sylph_species', 'sylph_ani', 'agreement']
        f.write('\t'.join(tax_cols) + '\n')
        for sample_name in sorted(gtdbtk_data.keys()):
            row = gtdbtk_data[sample_name]
            iso_id, asm_type = split_sample(sample_name)
            vals = [
                run_id, workflow_desc, workflow_timestamp, group_id, group_desc, sample_name, iso_id, asm_type,
                row.get('Domain', ''), row.get('Phylum', ''), row.get('Class', ''),
                row.get('Order', ''), row.get('Family', ''), row.get('Genus', ''),
                row.get('Species', ''), row.get('Closest_Reference', ''),
                str(row.get('ANI', '')), str(row.get('AF', '')),
                row.get('Method', ''),
                row.get('Call', '') if row.get('Call') != '-' else '',
                row.get('Second_Species', '') if row.get('Second_Species') != '-' else '',
                row.get('Second_ANI', '') if row.get('Second_ANI') != '-' else '',
                (sylph_compare.get(sample_name, {}).get('species') or sylph_compare.get(sample_name, {}).get('genus') or ''),
                sylph_compare.get(sample_name, {}).get('ani', ''),
                taxonomy_agreement(sample_name) if taxonomy_tool == 'both' else '',
            ]
            f.write('\t'.join(vals) + '\n')

    with open('mlst.tsv', 'w') as f:
        mlst_cols = ['workflow_id', 'workflow_description', 'workflow_timestamp', 'group_id', 'group_description', 'sample', 'isolate_id', 'assembly_type', 'scheme', 'st'] + [f'locus_{i}' for i in range(1, 8)]
        f.write('\t'.join(mlst_cols) + '\n')
        for sample_name in sorted(mlst_data.keys()):
            row = mlst_data[sample_name]
            iso_id, asm_type = split_sample(sample_name)
            vals = [run_id, workflow_desc, workflow_timestamp, group_id, group_desc, sample_name, iso_id, asm_type, row.get('Scheme', ''), row.get('ST', '')]
            for i in range(1, 8):
                vals.append(row.get(f'Locus_{i}', ''))
            f.write('\t'.join(vals) + '\n')

    # amrfinder.tsv: every AMRFinderPlus hit, unfiltered, one row per hit.
    def db_col(name):
        return (name.strip().lower().replace(' ', '_').replace('%', 'pct')
                .replace('/', '_').replace('-', '_'))

    provenance = ['workflow_id', 'workflow_description', 'workflow_timestamp',
                  'group_id', 'group_description', 'sample', 'isolate_id', 'assembly_type']

    def provenance_values(sample_name):
        iso_id, asm_type = split_sample(sample_name)
        return [run_id, workflow_desc, workflow_timestamp, group_id, group_desc,
                sample_name, iso_id, asm_type]

    amrfinder_cols = None
    with open('amrfinder.tsv', 'w') as f:
        for sample_name in sorted(amrfinder_rows.keys()):
            for row in amrfinder_rows[sample_name]:
                if amrfinder_cols is None:
                    # AMRFinderPlus already carries Name; sample/isolate_id/assembly_type replace it because the
                    # filename is the authoritative sample id.
                    amrfinder_cols = [c for c in row.keys() if c != 'Name']
                    f.write('\t'.join(provenance + [db_col(c) for c in amrfinder_cols]) + '\n')
                vals = provenance_values(sample_name) + [str(row.get(c, '') or '') for c in amrfinder_cols]
                f.write('\t'.join(vals) + '\n')
        if amrfinder_cols is None:
            f.write('\t'.join(provenance) + '\n')

    # resfinder_phenotype.tsv: the complete phenotype panel, every antimicrobial whether resistant
    # or not.
    with open('resfinder_phenotype.tsv', 'w') as f:
        f.write('\t'.join(provenance + ['antimicrobial', 'antimicrobial_class',
                                        'wgs_predicted_phenotype', 'match',
                                        'genetic_background']) + '\n')
        for sample_name in sorted(resfinder_pheno.keys()):
            for row in resfinder_pheno[sample_name]:
                f.write('\t'.join(provenance_values(sample_name) + [
                    row.get('Antimicrobial', ''),
                    row.get('Class', ''),
                    row.get('WGS-predicted phenotype', ''),
                    row.get('Match', ''),
                    row.get('Genetic background', ''),
                ]) + '\n')

    # resfinder_genes.tsv: every ResFinder acquired-gene hit, with the contig it sits on resolved
    # to a molecule type so the database can filter on location.
    with open('resfinder_genes.tsv', 'w') as f:
        f.write('\t'.join(provenance + ['resistance_gene', 'identity', 'coverage',
                                        'alignment_length_gene_length',
                                        'position_in_reference', 'contig',
                                        'molecule_type', 'position_in_contig',
                                        'phenotype', 'accession']) + '\n')
        for sample_name in sorted(resfinder_raw.keys()):
            contig_map = sample_contig_map.get(sample_name, {})
            for row in resfinder_raw[sample_name]:
                contig_full = (row.get('Contig') or '').strip()
                contig = contig_full.split()[0] if contig_full else ''
                f.write('\t'.join(provenance_values(sample_name) + [
                    row.get('Resistance gene', ''),
                    row.get('Identity', ''),
                    row.get('Coverage', ''),
                    row.get('Alignment Length/Gene Length', ''),
                    row.get('Position in reference', ''),
                    contig,
                    contig_map.get(contig, {}).get('molecule_type', ''),
                    row.get('Position in contig', ''),
                    row.get('Phenotype', ''),
                    row.get('Accession no.', ''),
                ]) + '\n')

    # pointfinder.tsv: every PointFinder hit, unfiltered.
    with open('pointfinder.tsv', 'w') as f:
        f.write('\t'.join(provenance + ['mutation', 'mutation_id', 'nucleotide_change',
                                        'amino_acid_change', 'resistance', 'pmid']) + '\n')
        for sample_name in sorted(pointfinder_hits.keys()):
            for row in pointfinder_hits[sample_name]:
                f.write('\t'.join(provenance_values(sample_name) + [
                    row.get('Mutation', ''),
                    format_point_mutation(row.get('Mutation', '')),
                    row.get('Nucleotide change', ''),
                    row.get('Amino acid change', ''),
                    row.get('Resistance', ''),
                    row.get('PMID', ''),
                ]) + '\n')

    # point_mutations.tsv: POINT rows resolved to their contig and molecule type, so a mutation
    # can be traced to chromosome or plasmid the same way genes are.
    with open('point_mutations.tsv', 'w') as f:
        f.write('\t'.join(provenance + ['mutation', 'contig', 'molecule_type',
                                        'amr_class', 'amr_subclass', 'method',
                                        'pct_identity', 'pct_coverage']) + '\n')
        for sample_name in sorted(amrfinder_rows.keys()):
            contig_map = sample_contig_map.get(sample_name, {})
            for row in amrfinder_rows[sample_name]:
                if (row.get('Subtype') or '').strip().upper() != AMRFINDER_SUBTYPE_POINT:
                    continue
                contig = (row.get('Contig id') or '').strip().split()[0] if row.get('Contig id') else ''
                mol = contig_map.get(contig, {}).get('molecule_type', '')
                f.write('\t'.join(provenance_values(sample_name) + [
                    row.get('Element symbol', ''),
                    contig,
                    mol,
                    row.get('Class', ''),
                    row.get('Subclass', ''),
                    row.get('Method', ''),
                    row.get('% Identity to reference', ''),
                    row.get('% Coverage of reference', ''),
                ]) + '\n')

    # plasmid_resistance_markers.tsv: contig-level with individual contig sizes
    plasmid_res_rows = []
    for sample_name in sorted(set(list(sample_contig_map.keys()) + list(sample_amr_genes.keys()))):
        contig_map = sample_contig_map.get(sample_name, {})
        res_contigs = sample_amr_genes.get(sample_name, {})
        mobtyper = sample_mobtyper.get(sample_name, {})

        for cid, info in sorted(contig_map.items()):
            if info['molecule_type'] != 'plasmid':
                continue
            genes = res_contigs.get(cid, [])
            cluster_id = info['cluster_id'] if info['cluster_id'] != '-' else 'unknown'
            meta = mobtyper.get(cluster_id, {})
            # Use individual contig size, not cluster size
            contig_size = info.get('contig_length', '0')
            iso_id, asm_type = split_sample(sample_name)
            plasmid_res_rows.append({
                'workflow_id': run_id,
                'workflow_description': workflow_desc,
                'workflow_timestamp': workflow_timestamp,
                'group_id': group_id,
                'group_description': group_desc,
                'sample': sample_name,
                'isolate_id': iso_id,
                'assembly_type': asm_type,
                'contig': cid,
                'cluster_id': cluster_id,
                'secondary_cluster_id': meta.get('secondary_cluster_id', '-'),
                'size': contig_size,
                'circularity': info.get('circularity', '-'),
                'num_resistance_genes': len(genes),
                'resistance_genes': '; '.join(sorted(genes)) if genes else '-',
                'rep_type': meta.get('rep_type', '-'),
                'relaxase_type': meta.get('relaxase_type', '-'),
                'mpf_type': meta.get('mpf_type', '-'),
                'orit_type': meta.get('orit_type', '-'),
                'predicted_mobility': meta.get('predicted_mobility', '-'),
                'closest_match': meta.get('mash_nearest', '-'),
            })

    with open('plasmid_resistance_markers.tsv', 'w') as f:
        pr_cols = ['workflow_id', 'workflow_description', 'workflow_timestamp', 'group_id', 'group_description', 'sample', 'isolate_id', 'assembly_type', 'contig', 'cluster_id', 'secondary_cluster_id', 'size', 'circularity', 'num_resistance_genes',
                   'resistance_genes', 'rep_type', 'relaxase_type', 'mpf_type', 'orit_type',
                   'predicted_mobility', 'closest_match']
        f.write('\t'.join(pr_cols) + '\n')
        for row in plasmid_res_rows:
            f.write('\t'.join(str(row.get(c, '')) for c in pr_cols) + '\n')

    # chromosomal_resistance_markers.tsv: contig-level with individual contig sizes
    chromo_res_rows = []
    for sample_name in sorted(set(list(sample_contig_map.keys()) + list(sample_amr_genes.keys()))):
        contig_map = sample_contig_map.get(sample_name, {})
        res_contigs = sample_amr_genes.get(sample_name, {})

        for cid, info in sorted(contig_map.items()):
            if info['molecule_type'] != 'chromosome':
                continue
            genes = res_contigs.get(cid, [])
            if genes:
                contig_size = info.get('contig_length', '0')
                iso_id, asm_type = split_sample(sample_name)
                chromo_res_rows.append({
                    'workflow_id': run_id,
                    'workflow_description': workflow_desc,
                    'workflow_timestamp': workflow_timestamp,
                    'group_id': group_id,
                    'group_description': group_desc,
                    'sample': sample_name,
                    'isolate_id': iso_id,
                    'assembly_type': asm_type,
                    'contig': cid,
                    'size': contig_size,
                    'circularity': info.get('circularity', '-'),
                    'num_resistance_genes': len(genes),
                    'resistance_genes': '; '.join(sorted(genes)),
                })

    with open('chromosomal_resistance_markers.tsv', 'w') as f:
        cr_cols = ['workflow_id', 'workflow_description', 'workflow_timestamp', 'group_id', 'group_description', 'sample', 'isolate_id', 'assembly_type', 'contig', 'size', 'circularity', 'num_resistance_genes', 'resistance_genes']
        f.write('\t'.join(cr_cols) + '\n')
        for row in chromo_res_rows:
            f.write('\t'.join(str(row.get(c, '')) for c in cr_cols) + '\n')

    # virulence_factors.tsv: per-hit detail from VirulenceFinder
    with open('virulence_factors.tsv', 'w') as f:
        vf_cols = ['workflow_id', 'workflow_description', 'workflow_timestamp', 'group_id', 'group_description',
                   'sample', 'isolate_id', 'assembly_type', 'database', 'virulence_factor', 'identity',
                   'query_template_length', 'contig', 'position', 'protein_function', 'accession']
        f.write('\t'.join(vf_cols) + '\n')
        for sample_name, hits in sorted(sample_virulence.items()):
            iso_id, asm_type = split_sample(sample_name)
            for hit in hits:
                vals = [run_id, workflow_desc, workflow_timestamp, group_id, group_desc,
                        sample_name, iso_id, asm_type,
                        hit['database'], hit['virulence_factor'], hit['identity'],
                        hit['query_template_length'], hit['contig'], hit['position'],
                        hit['protein_function'], hit['accession']]
                f.write('\t'.join(str(v) for v in vals) + '\n')

    # ── Species-Based Typing: Parse results and write TSVs + MultiQC sections ──
    from parse_species_typing import (parse_ectyper_dirs, parse_kleborate_dirs, parse_emmtyper_results,
                                       parse_sistr_dirs, parse_seroba_dirs, parse_meningotype_results,
                                       parse_elgato_dirs)

    ectyper_data = parse_ectyper_dirs('ectyper_dirs')
    kleborate_data = parse_kleborate_dirs('kleborate_dirs')
    emmtyper_data = parse_emmtyper_results('emmtyper_files')
    sistr_data = parse_sistr_dirs('sistr_dirs')
    seroba_data = parse_seroba_dirs('seroba_dirs')
    meningotype_data = parse_meningotype_results('meningotype_files')
    elgato_data = parse_elgato_dirs('elgato_dirs')

    # ECTyper TSV and MultiQC section
    ectyper_mqc_cols = ['Serotype', 'QC', 'Pathotype', 'PathotypeGenes', 'StxSubtypes']
    ectyper_mqc_display = {
        'Serotype': 'Serotype', 'QC': 'QC', 'Pathotype': 'Pathotype',
        'PathotypeGenes': 'Pathotype Genes', 'StxSubtypes': 'Stx Subtypes'
    }
    with open('ectyper.tsv', 'w') as f:
        # Get all columns from first result
        all_ec_cols = []
        for v in ectyper_data.values():
            all_ec_cols = list(v.keys())
            break
        header = ['workflow_id', 'workflow_description', 'workflow_timestamp', 'group_id', 'group_description',
                  'sample', 'isolate_id', 'assembly_type'] + (all_ec_cols if all_ec_cols else ectyper_mqc_cols)
        f.write('\t'.join(header) + '\n')
        for sample_name, row in sorted(ectyper_data.items()):
            iso_id, asm_type = split_sample(sample_name)
            vals = [run_id, workflow_desc, workflow_timestamp, group_id, group_desc,
                    sample_name, iso_id, asm_type]
            for col in (all_ec_cols if all_ec_cols else ectyper_mqc_cols):
                vals.append(str(row.get(col, '')))
            f.write('\t'.join(vals) + '\n')

    if ectyper_data:
        ec_table_data = {}
        for sample_name, row in ectyper_data.items():
            ec_table_data[sample_name] = {'Species': species_lookup.get(sample_name, '-')}
            ec_table_data[sample_name].update({ectyper_mqc_display.get(c, c): row.get(c, '') for c in ectyper_mqc_cols})
        write_yaml('ectyper_mqc.yaml', {
            'id': 'ectyper_typing',
            'section_name': 'Escherichia coli Typing',
            'description': 'Serotyping, pathotyping and Shiga toxin subtyping for E. coli isolates',
            'plot_type': 'table',
            'pconfig': {'id': 'ectyper_table', 'title': 'ECTyper Results', 'sortRows': True},
            'data': ec_table_data
        })
    else:
        open('ectyper_mqc.yaml', 'w').close()

    # Kleborate TSVs and MultiQC sections KPSC columns
    kpsc_mqc_cols = ['species', 'ST', 'virulence_score', 'resistance_score', 'K_type', 'K_locus',
                     'K_locus_confidence', 'O_locus', 'O_type', 'O_locus_confidence',
                     'Bla_Carb_acquired', 'Bla_ESBL_acquired',
                     'YbST', 'Yersiniabactin', 'AbST', 'Aerobactin',
                     'SmST', 'Salmochelin', 'RmST', 'RmpADC', 'rmpA2', 'rmpA', 'rmpD', 'rmpC']
    kpsc_mqc_display = {
        'species': 'Species', 'ST': 'ST', 'virulence_score': 'Virulence Score', 'resistance_score': 'Resistance Score',
        'K_type': 'K Type', 'K_locus': 'K Locus', 'K_locus_confidence': 'K Confidence',
        'O_locus': 'O Locus', 'O_type': 'O Type', 'O_locus_confidence': 'O Confidence',
        'Bla_Carb_acquired': 'Carbapenemases', 'Bla_ESBL_acquired': 'ESBLs',
        'YbST': 'YbST', 'Yersiniabactin': 'Yersiniabactin', 'AbST': 'AbST', 'Aerobactin': 'Aerobactin',
        'SmST': 'SmST', 'Salmochelin': 'Salmochelin',
        'RmST': 'RmST', 'RmpADC': 'RmpADC', 'rmpA2': 'rmpA2', 'rmpA': 'rmpA', 'rmpD': 'rmpD', 'rmpC': 'rmpC'
    }

    kosc_mqc_cols = ['species', 'ST', 'YbST', 'Yersiniabactin', 'CbST', 'Colibactin',
                     'AbST', 'Aerobactin', 'SmST', 'Salmochelin',
                     'RmST', 'RmpADC', 'rmpA2']
    kosc_mqc_display = {
        'species': 'Species', 'ST': 'ST',
        'YbST': 'YbST', 'Yersiniabactin': 'Yersiniabactin',
        'CbST': 'CbST', 'Colibactin': 'Colibactin',
        'AbST': 'AbST', 'Aerobactin': 'Aerobactin',
        'SmST': 'SmST', 'Salmochelin': 'Salmochelin',
        'RmST': 'RmST', 'RmpADC': 'RmpADC', 'rmpA2': 'rmpA2'
    }

    kpsc_data = {s: d for s, d in kleborate_data.items() if d.get('_preset') == 'kpsc'}
    kosc_data = {s: d for s, d in kleborate_data.items() if d.get('_preset') == 'kosc'}

    with open('kleborate.tsv', 'w') as f:
        # Get all columns from first result
        all_kleb_cols = []
        for v in kleborate_data.values():
            all_kleb_cols = [k for k in v.keys() if not k.startswith('_')]
            break
        header = ['workflow_id', 'workflow_description', 'workflow_timestamp', 'group_id', 'group_description',
                  'sample', 'isolate_id', 'assembly_type', 'preset'] + (all_kleb_cols if all_kleb_cols else kpsc_mqc_cols)
        f.write('\t'.join(header) + '\n')
        for sample_name, row in sorted(kleborate_data.items()):
            iso_id, asm_type = split_sample(sample_name)
            preset = row.get('_preset', '')
            vals = [run_id, workflow_desc, workflow_timestamp, group_id, group_desc,
                    sample_name, iso_id, asm_type, preset]
            for col in (all_kleb_cols if all_kleb_cols else kpsc_mqc_cols):
                vals.append(str(row.get(col, '')))
            f.write('\t'.join(vals) + '\n')

    if kpsc_data:
        kpsc_table_data = {}
        for sample_name, row in kpsc_data.items():
            kpsc_table_data[sample_name] = {kpsc_mqc_display.get(c, c): row.get(c, '') for c in kpsc_mqc_cols}
        write_yaml('kleborate_kpsc_mqc.yaml', {
            'id': 'kleborate_kpsc_typing',
            'section_name': 'Klebsiella pneumoniae Species Complex Typing',
            'description': '<p style="margin:0 0 8px 0;">KPSC serotyping, virulence and resistance scoring</p>'
                          '<details style="margin:0;"><summary><b>Scoring Reference</b> <small>(click to expand)</small></summary>'
                          '<table style="font-size:12px; border-collapse:collapse; margin:8px 0;">'
                          '<tr><th colspan="2" style="text-align:left; padding:4px 8px; border-bottom:2px solid #ccc;">Virulence Score</th></tr>'
                          '<tr><td style="padding:2px 8px; font-weight:bold;">0</td><td style="padding:2px 8px;">Negative for yersiniabactin (ybt), colibactin (clb), aerobactin (iuc)</td></tr>'
                          '<tr><td style="padding:2px 8px; font-weight:bold;">1</td><td style="padding:2px 8px;">Yersiniabactin only</td></tr>'
                          '<tr><td style="padding:2px 8px; font-weight:bold;">2</td><td style="padding:2px 8px;">Yersiniabactin and colibactin (or colibactin only)</td></tr>'
                          '<tr><td style="padding:2px 8px; font-weight:bold;">3</td><td style="padding:2px 8px;">Aerobactin without yersiniabactin or colibactin</td></tr>'
                          '<tr><td style="padding:2px 8px; font-weight:bold;">4</td><td style="padding:2px 8px;">Aerobactin with yersiniabactin (without colibactin)</td></tr>'
                          '<tr><td style="padding:2px 8px; font-weight:bold;">5</td><td style="padding:2px 8px;">Yersiniabactin, colibactin and aerobactin</td></tr>'
                          '<tr><th colspan="2" style="text-align:left; padding:8px 8px 2px; border-bottom:2px solid #ccc;">Resistance Score</th></tr>'
                          '<tr><td style="padding:2px 8px; font-weight:bold;">0</td><td style="padding:2px 8px;">No ESBL, no carbapenemase (regardless of colistin resistance)</td></tr>'
                          '<tr><td style="padding:2px 8px; font-weight:bold;">1</td><td style="padding:2px 8px;">ESBL, no carbapenemase (regardless of colistin resistance)</td></tr>'
                          '<tr><td style="padding:2px 8px; font-weight:bold;">2</td><td style="padding:2px 8px;">Carbapenemase without colistin resistance (regardless of ESBL/OmpK)</td></tr>'
                          '<tr><td style="padding:2px 8px; font-weight:bold;">3</td><td style="padding:2px 8px;">Carbapenemase with colistin resistance (regardless of ESBL/OmpK)</td></tr>'
                          '<tr><th colspan="2" style="text-align:left; padding:8px 8px 2px; border-bottom:2px solid #ccc;">Hypermucoidy Loci</th></tr>'
                          '<tr><td colspan="2" style="padding:2px 8px;">rmpADC (rmpA, rmpD, rmpC), rmpA2</td></tr>'
                          '</table></details>',
            'plot_type': 'table',
            'pconfig': {'id': 'kleborate_kpsc_table', 'title': 'Kleborate KPSC Results', 'sortRows': True},
            'data': kpsc_table_data
        })
    else:
        open('kleborate_kpsc_mqc.yaml', 'w').close()

    if kosc_data:
        kosc_table_data = {}
        for sample_name, row in kosc_data.items():
            kosc_table_data[sample_name] = {kosc_mqc_display.get(c, c): row.get(c, '') for c in kosc_mqc_cols}
        write_yaml('kleborate_kosc_mqc.yaml', {
            'id': 'kleborate_kosc_typing',
            'section_name': 'Klebsiella oxytoca Species Complex Typing',
            'description': 'KOSC virulence typing',
            'plot_type': 'table',
            'pconfig': {'id': 'kleborate_kosc_table', 'title': 'Kleborate KOSC Results', 'sortRows': True},
            'data': kosc_table_data
        })
    else:
        open('kleborate_kosc_mqc.yaml', 'w').close()

    # emmtyper TSV and MultiQC section
    with open('emmtyper.tsv', 'w') as f:
        header = ['workflow_id', 'workflow_description', 'workflow_timestamp', 'group_id', 'group_description',
                  'sample', 'isolate_id', 'assembly_type', 'species',
                  'emm_type', 'emm_cluster', 'num_hits', 'num_clusters',
                  'possible_emm_2', 'position_1', 'position_2']
        f.write('\t'.join(header) + '\n')
        for sample_name, row in sorted(emmtyper_data.items()):
            iso_id, asm_type = split_sample(sample_name)
            vals = [run_id, workflow_desc, workflow_timestamp, group_id, group_desc,
                    sample_name, iso_id, asm_type,
                    species_lookup.get(sample_name, '-'),
                    row.get('emm_type', ''), row.get('emm_cluster', ''),
                    row.get('num_hits', ''), row.get('num_clusters', ''),
                    row.get('possible_emm_2', ''), row.get('position_1', ''), row.get('position_2', '')]
            f.write('\t'.join(vals) + '\n')

    if emmtyper_data:
        emm_table_data = {}
        for sample_name, row in emmtyper_data.items():
            emm_table_data[sample_name] = {
                'Species': species_lookup.get(sample_name, '-'),
                'emm Type': row.get('emm_type', ''),
                'emm Cluster': row.get('emm_cluster', ''),
            }
        write_yaml('emmtyper_mqc.yaml', {
            'id': 'emmtyper_typing',
            'section_name': 'Streptococcus pyogenes emm Typing',
            'description': 'Typing and clustering for S. pyogenes isolates.',
            'plot_type': 'table',
            'pconfig': {'id': 'emmtyper_table', 'title': 'emmtyper Results', 'sortRows': True},
            'data': emm_table_data
        })
    else:
        open('emmtyper_mqc.yaml', 'w').close()

    # ── Additional species-based typing tools (SISTR, SeroBA, meningotype, el_gato) ──
    def _emit_typing(data, tsv_name, yaml_name, section_id, section_name, description, pconfig_id, table_title, col_display):
        """Write a full <tool>.tsv (workflow metadata + all raw columns) and a MultiQC custom-
        content table section (<tool>_mqc.yaml) using selected, renamed columns."""
        raw_cols = [c for c, _ in col_display]
        with open(tsv_name, 'w') as f:
            all_cols = []
            for v in data.values():
                all_cols = list(v.keys())
                break
            meta_cols = ['workflow_id', 'workflow_description', 'workflow_timestamp',
                         'group_id', 'group_description', 'sample', 'isolate_id',
                         'assembly_type']
            # Tool output often carries its own sample column ('Sample' in seroba, 'Name' elsewhere).
            _meta_lower = {c.lower() for c in meta_cols}
            tool_cols = [c for c in (all_cols if all_cols else raw_cols)
                         if c.strip().lower() not in _meta_lower]
            # ECTyper and Kleborate report a species themselves, so their tables carry one already.
            add_species = not any(c.strip().lower() == 'species' for c in tool_cols)
            header = meta_cols + (['species'] if add_species else []) + tool_cols
            f.write('\t'.join(header) + '\n')
            for sample_name, row in sorted(data.items()):
                iso_id, asm_type = split_sample(sample_name)
                vals = [run_id, workflow_desc, workflow_timestamp, group_id, group_desc,
                        sample_name, iso_id, asm_type]
                if add_species:
                    vals.append(species_lookup.get(sample_name, '-'))
                for col in tool_cols:
                    vals.append(str(row.get(col, '')))
                f.write('\t'.join(vals) + '\n')
        if data:
            table_data = {}
            for sample_name, row in data.items():
                entry = {'Species': species_lookup.get(sample_name, '-')}
                for raw, disp in col_display:
                    entry[disp] = row.get(raw, '')
                table_data[sample_name] = entry
            write_yaml(yaml_name, {
                'id': section_id,
                'section_name': section_name,
                'description': description,
                'plot_type': 'table',
                'pconfig': {'id': pconfig_id, 'title': table_title, 'sortRows': True},
                'data': table_data
            })
        else:
            open(yaml_name, 'w').close()

    _emit_typing(sistr_data, 'sistr.tsv', 'sistr_mqc.yaml', 'sistr_typing',
                 'Salmonella Serovar Prediction',
                 'cgMLST-based serovar and antigen prediction for Salmonella enterica isolates.',
                 'sistr_table', 'Salmonella Serovar Prediction',
                 [('cgmlst_subspecies', 'Subspecies'), ('serogroup', 'Serogroup'),
                  ('serovar_cgmlst', 'Serovar'), ('h1', 'H1 - Phase1 Flagellar Antigen'),
                  ('h2', 'H2 - Phase2 Flagellar Antigen'), ('o_antigen', 'O Antigen'),
                  ('qc_status', 'QC')])

    _emit_typing(seroba_data, 'seroba.tsv', 'seroba_mqc.yaml', 'seroba_typing',
                 'Streptococcus pneumoniae Serotyping',
                 'Capsular serotype prediction from Illumina reads for S. pneumoniae.',
                 'seroba_table', 'Streptococcus pneumoniae Serotyping',
                 [('Serotype', 'Serotype'), ('Genetic_Variant', 'Genetic Variant'),
                  ('Contamination_Status', 'Contamination Status')])

    _emit_typing(meningotype_data, 'meningotype.tsv', 'meningotype_mqc.yaml', 'meningotype_typing',
                 'Neisseria meningitidis Typing',
                 'Serogroup and fine typing (PorA, FetA, MLST, Bexsero antigens) for N. meningitidis.',
                 'meningotype_table', 'Neisseria meningitidis Typing',
                 [('SEROGROUP', 'Serogroup'), ('MLST', 'MLST'), ('ctrA', 'ctrA'), ('PorA', 'PorA'),
                  ('FetA', 'FetA'), ('PorB', 'PorB'), ('fHbp', 'fHbp'), ('NHBA', 'NHBA'),
                  ('NadA', 'NadA'), ('BAST', 'BAST')])

    _emit_typing(elgato_data, 'elgato.tsv', 'elgato_mqc.yaml', 'elgato_typing',
                 'Legionella pneumophila Sequence Typing',
                 'Sequence-based typing (SBT) for L. pneumophila.',
                 'elgato_table', 'Legionella pneumophila Sequence Typing',
                 [('ST', 'MLST'), ('flaA', 'flaA'), ('pilE', 'pilE'), ('asd', 'asd'),
                  ('mip', 'mip'), ('mompS', 'mompS'), ('proA', 'proA'), ('neuA_neuAH', 'neuA/neuAH')])

    # ── Custom CSS ──
    CSS_CONTENT = '''/* Custom Styles for Hybrid Assembly Pipeline */

/* ===== GLOBAL THEME ===== */
body { background-color: #fffefe; }
.side-nav { border-right: #ffffff; }
.side-nav-wrapper { background-color: #ffffff; }
a { color: rgb(115 119 127); }
div#mqc-nav-header { display: none !important; }
.side-nav .side-nav-title a { display: none !important; }
div#mqc-warning-many-samples { display: none !important; }
#mqc_loading_warning { display: none !important; }
.text-muted { display: none; }
/* ...but MultiQC greys the config-modal row of a HIDDEN column with .text-muted,
   so the blanket rule above erased exactly those rows: "Configure columns"
   listed only the visible columns, and a hidden one could not be re-enabled
   until "Show All" was pressed. Restore anything .text-muted inside a table or
   a config modal; the blanket rule still applies to the report chrome it was
   written for. */
tr.text-muted { display: table-row; }
td.text-muted, th.text-muted { display: table-cell; }
.mqc_config_modal .text-muted { display: revert; }

/* ===== LOGO REPLACEMENT ===== */
.side-nav-logo, .side-nav-logo a, .side-nav-logo svg, .side-nav h1.side-nav-logo { display: none !important; }
#page_title .multiqc-logo-wrapper, #page_title .multiqc-logo, #page_title a[href*="multiqc.info"] { display: none !important; }

#page_title::before {
    content: "";
    display: block;
    width: 220px;
    height: 50px;
    margin-bottom: 15px;
    background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 195 45'%3E%3Cg id='icon-group' stroke='none' fill='%23D85A3B' transform='matrix(0.65, 0, 0, 0.65, 0, -8)'%3E%3Cpath d='M72 22H28C12.56 22 0 34.564 0 50c0 13.236 10.764 24 24 24h48c13.236 0 24-10.764 24-24v-4c0-13.236-10.764-24-24-24zM57.916 46.752a379.337 379.337 0 01-6.132 4.1C59.16 55.24 67.176 58.84 76 59.716V42c0-1.1.9-2 2-2s2 .9 2 2v18c1.1 0 2 .9 2 2s-.9 2-2 2c-12.132 0-22.592-4.92-32-10.736C38.592 59.072 28.132 64 16 64c-1.104 0-2-.9-2-2s.896-2 2-2c4.192 0 8.188-.7 12.04-1.813-.008-.063-.04-.115-.04-.187v-8c0-1.1.896-2 2-2s2 .9 2 2v6.9c4.272-1.601 8.328-3.729 12.22-6.045a567.82 567.82 0 01-6.136-4.096C31.34 42.168 25.336 38.124 20 36.64V50c0 1.1-.896 2-2 2s-2-.9-2-2V36c-1.104 0-2-.9-2-2s.896-2 2-2c7.504 0 15.304 5.304 24.332 11.444 2.476 1.684 5.04 3.412 7.668 5.083 2.628-1.679 5.191-3.407 7.668-5.083C64.695 37.304 72.496 32 80 32c1.1 0 2 .9 2 2s-.9 2-2 2c-4.721 0-10.04 2.824-16 6.68V50c0 1.1-.9 2-2 2s-2-.9-2-2v-4.656c-.687.468-1.376.928-2.084 1.408z'%3E%3C/path%3E%3C/g%3E%3Cg id='text-group' fill='%23D85A3B' transform='matrix(1.8, 0, 0, 1.8, 68, 4)'%3E%3Cpath d='M6.57 13.87C5.52 14.10 3.92 14.32 2.60 14.35L2.59 10.30L5.94 10.30C5.99 10.08 6.01 9.81 5.94 9.56L2.59 9.56L2.58 5.99C3.89 6.02 5.12 6.17 6.29 6.43C6.65 6.27 6.80 5.31 6.43 5.31L2.23 5.31C1.68 5.31 1.41 5.38 1.12 5.64C1.27 6.09 1.32 6.59 1.32 7.15L1.32 14.32C1.32 14.77 1.62 15.02 2.06 15.02L6.50 15.02C6.93 15.02 6.90 14.08 6.57 13.87ZM10.77 5.12C9.63 5.12 9.03 5.28 8.68 5.45C8.83 5.91 8.88 6.71 8.88 7.01L8.88 14.49C8.88 14.92 9.18 15.19 9.62 15.19C9.98 15.19 10.25 15.01 10.36 14.90C10.19 14.48 10.16 13.72 10.16 13.40L10.16 11.17L10.53 11.17C13.36 11.17 14.83 10.09 14.83 7.99C14.83 6.01 13.48 5.12 10.77 5.12ZM10.54 10.64L10.16 10.64L10.16 5.70C10.29 5.67 10.49 5.67 10.74 5.67C12.67 5.67 13.43 6.52 13.43 8.08C13.43 9.72 12.59 10.64 10.54 10.64ZM17.51 15.19C17.88 15.19 18.13 15.01 18.26 14.90C18.09 14.45 18.06 13.89 18.06 13.34L18.06 5.85C18.06 5.42 17.74 5.17 17.30 5.17C16.95 5.17 16.69 5.33 16.58 5.45C16.73 5.91 16.77 6.47 16.77 7.01L16.77 14.49C16.77 14.92 17.08 15.19 17.51 15.19ZM25.17 18.19C25.65 18.19 25.93 17.92 25.93 17.50L25.93 8.64C25.37 8.11 24.29 7.70 23.06 7.70C20.68 7.70 19.89 9.73 19.89 11.66C19.89 14.15 20.83 15.22 22.33 15.22C23.46 15.22 24.32 14.50 24.65 13.50L24.65 15.20C24.65 16.14 24.67 17.19 24.44 17.89C24.60 18.05 24.82 18.20 25.17 18.19ZM22.65 14.48C21.70 14.48 21.25 13.57 21.25 11.70C21.24 10.14 21.50 8.30 23.11 8.30C23.77 8.30 24.39 8.55 24.67 8.79L24.67 11.24C24.67 12.52 24.09 14.48 22.65 14.48ZM34.31 14C34.22 14.17 34.05 14.20 33.94 14.20C33.67 14.20 33.56 14 33.56 13.48L33.54 8.46C33.54 7.98 33.22 7.70 32.82 7.70C32.45 7.70 32.20 7.88 32.09 7.98C32.30 8.57 32.30 9.06 32.30 11.23C32.30 13.61 31.35 14.55 30.46 14.55C29.47 14.55 29.15 13.85 29.15 12.12L29.15 8.40C29.15 7.95 28.84 7.70 28.42 7.70C28.08 7.70 27.86 7.85 27.69 7.99C27.85 8.44 27.89 9.07 27.89 9.55L27.89 12.29C27.89 14.43 28.59 15.22 30.09 15.22C31.21 15.22 32.05 14.59 32.42 13.48L32.42 13.66C32.42 14.67 32.84 15.16 33.59 15.16C34.15 15.16 34.47 14.88 34.47 14.46C34.47 14.31 34.43 14.15 34.31 14ZM41.09 13.17C40.80 13.93 40.01 14.36 39.00 14.36C37.79 14.36 36.85 13.72 36.72 11.79C37.21 11.84 37.70 11.87 38.21 11.87C40.14 11.87 41.36 11.14 41.36 9.74C41.36 8.40 40.25 7.70 38.74 7.70C36.54 7.70 35.35 9.16 35.35 11.48C35.35 13.90 36.64 15.20 38.68 15.20C39.89 15.20 40.85 14.76 41.20 13.89C41.40 13.43 41.27 13.24 41.09 13.17ZM38.68 8.23C39.33 8.23 40.00 8.51 40.00 9.56C40.00 11.00 38.88 11.33 37.65 11.33C37.32 11.31 37.00 11.30 36.71 11.27C36.72 9.28 37.49 8.23 38.68 8.23ZM48.99 13.34L48.99 10.60C48.99 8.69 48.27 7.70 46.73 7.70C45.58 7.70 44.74 8.32 44.35 9.39L44.35 9.23C44.35 8.22 43.93 7.71 43.19 7.71C42.63 7.71 42.29 8.01 42.29 8.43C42.29 8.60 42.35 8.74 42.46 8.89C42.55 8.72 42.73 8.69 42.83 8.69C43.09 8.69 43.20 8.89 43.22 9.41L43.22 14.43C43.22 14.91 43.53 15.19 43.96 15.19C44.31 15.19 44.58 15.02 44.69 14.91C44.48 14.34 44.48 13.83 44.48 11.66C44.48 9.32 45.46 8.41 46.35 8.41C47.39 8.41 47.71 9.24 47.73 10.77L47.73 14.49C47.73 14.92 48.05 15.19 48.47 15.19C48.79 15.19 49.03 15.05 49.20 14.90C49.03 14.45 49 13.82 48.99 13.34ZM53.93 15.20C55.19 15.20 56.15 14.67 56.46 13.68C56.60 13.20 56.49 13.01 56.31 12.92C56.08 13.86 55.27 14.36 54.28 14.36C52.98 14.36 51.98 13.54 51.98 11.28C51.98 9.28 52.77 8.23 53.90 8.23C54.60 8.23 55.17 8.55 55.17 9.38C55.17 9.83 55.03 10.16 54.73 10.30C54.89 10.51 55.17 10.60 55.44 10.60C56.03 10.60 56.45 10.15 56.45 9.48C56.45 8.47 55.61 7.70 54.03 7.70C51.84 7.70 50.62 9.14 50.62 11.48C50.62 13.90 51.86 15.20 53.93 15.20ZM63.46 13.17C63.17 13.93 62.38 14.36 61.38 14.36C60.16 14.36 59.22 13.72 59.09 11.79C59.58 11.84 60.07 11.87 60.58 11.87C62.51 11.87 63.73 11.14 63.73 9.74C63.73 8.40 62.62 7.70 61.11 7.70C58.91 7.70 57.72 9.16 57.72 11.48C57.72 13.90 59.01 15.20 61.05 15.20C62.26 15.20 63.22 14.76 63.57 13.89C63.77 13.43 63.64 13.24 63.46 13.17ZM61.05 8.23C61.70 8.23 62.37 8.51 62.37 9.56C62.37 11.00 61.25 11.33 60.02 11.33C59.70 11.31 59.37 11.30 59.08 11.27C59.09 9.28 59.86 8.23 61.05 8.23Z'%3E%3C/path%3E%3C/g%3E%3C/svg%3E");
    background-repeat: no-repeat;
    background-size: contain;
}

/* ===== HIDE AI BUTTONS & ELEMENTS ===== */
.ai-generate-button-wrapper { display: none !important; }
.ai-copy-button-wrapper { display: none !important; }
.ai-plot-buttons-container { display: none !important; }
[title="AI"] { display: none !important; }

/* ===== HIDE FOOTER ===== */
footer, .footer { display: none !important; }

/* ===== RENAME QUAST AND BUSCO SECTIONS ===== */
h2#quast { font-size: 0 !important; line-height: 0 !important; }
h2#quast::after { content: "Assembly Quality Overview"; font-size: 24px; line-height: normal; }
#mqc-module-section-quast > .text-muted { display: none !important; }

h2#busco { font-size: 0 !important; line-height: 0 !important; }
h2#busco::after { content: "Assembly Completeness"; font-size: 24px; line-height: normal; }

/* Sidebar nav renaming */
.side-nav a[href="#quast"] { font-size: 0; }
.side-nav a[href="#quast"]::after { content: "Assembly Quality Overview"; font-size: 13px; }
.side-nav a[href="#busco"] { font-size: 0; }
.side-nav a[href="#busco"]::after { content: "Assembly Completeness"; font-size: 13px; }

/* ===== GENERAL TABLE STYLING ===== */
.mqc_table { font-size: 12px; }
.mqc_table td { max-width: none; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; padding: 4px 8px; }
.mqc_table td:first-child { font-weight: 600; }
#general_stats_table td { max-width: 100px; }
#general_stats_table th { max-width: 100px; font-size: 11px; }
#mlst_results_table td { max-width: 120px; }
#taxonomy_table td { max-width: 150px; }

/* ===== WIDE TABLES - override inherited .mqc_table td constraints ===== */
#amr_overview_table td,
#plasmid_overview_table td,
#plasmid_resistance_table td,
#chromosome_resistance_table td {
    max-width: none !important;
    white-space: normal !important;
    word-wrap: break-word !important;
    overflow: visible !important;
    text-overflow: unset !important;
}
#amr_overview_table,
#plasmid_overview_table,
#plasmid_resistance_table,
#chromosome_resistance_table {
    table-layout: auto !important;
    width: 100% !important;
}
.mqc_section h1 { border-bottom: 3px solid #3498db; padding-bottom: 10px; }
.mqc_table thead th { background-color: #f5f5f5; font-weight: 600; border-bottom: 2px solid #ddd; }
.mqc_table tbody tr:hover { outline: 2px solid #1976d2; outline-offset: -1px; }
.mqc-welcome, #mqc_welcome, .intro-panel { display: none !important; }
.mqc_table a { color: rgb(115 119 127); text-decoration: none; }
.mqc_table a:hover { text-decoration: underline; }

/* ===== AI SUMMARY: hide interactive buttons (summary is pre-generated) ===== */
.ai-generate-button-wrapper { display: none !important; }
.ai-copy-button-wrapper { display: none !important; }
.ai-plot-buttons-container { display: none !important; }
button.ai-generate-button-more, .ai-generate-button-more { display: none !important; }
[title="AI"] { display: none !important; }

/* ===== AI SUMMARY: native <details> accordion (no JS) ===== */
/* bin/ai_collapse.py rewrites MultiQC's AI wrapper into a details/summary element. */
#global_ai_summary_wrapper {
    margin-bottom: 26px;   /* breathing room before the next section (e.g. General Statistics) */
}
#global_ai_summary_wrapper > summary.ai-summary-header {
    cursor: pointer;
    font-size: 15px;
    font-weight: 700;
    color: #C9530B;
    margin-bottom: 4px;
    list-style: none;
}
/* The outbreak section injects a document-wide rule that hides native details
   markers, so this accordion draws its own caret (the span in ai_collapse.py)
   and rotates it on open, exactly like the outbreak and HGT accordions. These
   rules are scoped here so the caret still works when the outbreak arm is off
   and that global rule is never emitted. */
#global_ai_summary_wrapper > summary.ai-summary-header::-webkit-details-marker {
    display: none;
}
#global_ai_summary_wrapper > summary.ai-summary-header .ai-caret {
    display: inline-block;
    margin-right: 8px;
    font-size: inherit;   /* track the 15px header rather than being fixed smaller */
    line-height: 1;
    transition: transform 0.2s;
}
#global_ai_summary_wrapper[open] > summary.ai-summary-header .ai-caret {
    transform: rotate(90deg);
}
/* smaller, tighter body text + line spacing */
#global_ai_summary_response,
#global_ai_summary_detailed_analysis_response {
    font-size: 12px;
    line-height: 1.3;
}
/* disclaimer banner (real HTML, selectable) */
#global_ai_summary_wrapper .ai-disclaimer {
    font-size: 12px;
    line-height: 1.45;
    color: #8a5a1b;
    background: #fff7ed;
    border-left: 4px solid #F5731F;
    padding: 8px 10px;
    margin: 8px 0 12px 0;
}
/* Genes responsible for an R call, shown beside (not inside) the badge in the
   Antimicrobial Resistance Profile table. Deliberately not .text-muted, which
   this stylesheet hides. */
.amr-genes {
    font-size: 11px;
    color: #444;
    margin-left: 4px;
    white-space: normal;
}

/* provider/endpoint hidden, model name kept */
.ai-summary-disclaimer { font-size: 0 !important; }
.ai-summary-disclaimer-provider { display: none !important; }
.ai-summary-disclaimer-model { font-size: 12px !important; }
.ai-summary-disclaimer-model::before { content: "Model: "; font-size: 12px; }
/* bulleted structure inside the summary */
#global_ai_summary_response ul { margin: 2px 0 6px 20px; list-style: disc; }
#global_ai_summary_response li { margin-bottom: 2px; }
#global_ai_summary_response h1,
#global_ai_summary_response h2,
#global_ai_summary_response h3,
#global_ai_summary_response h4 {
    font-size: 13px; font-weight: 700; color: #C9530B; margin: 8px 0 3px 0;
}
'''
    with open('custom_styles.css', 'w') as f:
        f.write(CSS_CONTENT)

    print(f"Generated MultiQC config with {len(samples)} samples")
    print(f"  - MLST data: {len(mlst_data)} samples")
    print(f"  - Taxonomy data: {len(gtdbtk_data)} samples")
    print(f"  - AMR gene data: {len(sample_amr_genes)} samples")
    print(f"  - ResFinder phenotype data: {len(resfinder_pheno)} samples")
    print(f"  - Plasmid data: {len(mob_data)} samples")
    print(f"  - Plasmid resistance markers: {len(plasmid_markers)} entries")
    print(f"  - Chromosome resistance markers: {len(chromo_markers)} entries")

if __name__ == '__main__':
    main()