#!/usr/bin/env python3
"""Select reference sample and group contigs for plasmid comparison maps."""
import argparse
import csv
import fnmatch
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict


# FASTA / AMRFinderPlus parsing

CIRC_RE = re.compile(
    r'circular[= _]true|topology[= _]circular|\[circularised\]'
    r'|circular_contig|rotated[= _][Tt]rue|\bcircular\b',
    re.IGNORECASE,
)


def parse_fasta_contigs(path):
    contigs = []
    cur_name = cur_extra = None
    cur_lines = []

    def _flush():
        if cur_name is not None:
            seq = "".join(cur_lines)
            contigs.append({
                'name': cur_name, 'extra': cur_extra, 'seq': seq,
                'is_circ': bool(CIRC_RE.search(cur_extra)),
                'length': len(seq),
            })

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                _flush()
                parts = line[1:].strip().split(None, 1)
                cur_name = parts[0]
                cur_extra = parts[1] if len(parts) > 1 else ""
                cur_lines = []
            else:
                cur_lines.append(line.strip())
    _flush()
    return contigs


def load_amrfinder(path):
    """Return ``{contig: [{gene, pident, coverage}]}`` from an AMRFinderPlus TSV."""
    by_contig = defaultdict(list)
    if not path or not os.path.exists(path):
        return by_contig
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if (row.get('Subtype') or '').strip().upper() != 'AMR':
                continue
            if (row.get('Scope') or '').strip().lower() == 'plus':
                continue
            gene = (row.get('Element symbol') or '').strip()
            if not gene:
                continue
            contig = (row.get('Contig id') or '').strip()
            contig = contig.split()[0] if contig else ''
            pident = float(row.get('% Identity to reference') or 0)
            coverage = float(row.get('% Coverage of reference') or 0)
            by_contig[contig].append({'gene': gene, 'pident': pident, 'coverage': coverage})
    return by_contig


def find_amrfinder(sample_name, search_dir):
    p = os.path.join(search_dir, f"{sample_name}_amrfinder.tsv")
    if os.path.exists(p):
        return p
    m = re.match(r'^(.+?)_plasmid_', sample_name)
    if m:
        p = os.path.join(search_dir, f"{m.group(1)}_amrfinder.tsv")
        if os.path.exists(p):
            return p
    for f in sorted(os.listdir(search_dir)):
        if f.endswith('_amrfinder.tsv'):
            prefix = f.replace('_amrfinder.tsv', '')
            if sample_name.startswith(prefix):
                return os.path.join(search_dir, f)
    return None


# Containment / similarity checks (minimap2)

# minimap2 used to be optional here: a missing binary was caught and redundancy detection
# silently switched off, which changed which contigs reached the comparison maps depending on
# the host.

def _minimap2_missing():
    sys.exit('ERROR: minimap2 not found on PATH. SELECT_PLASMID_REF must run in the '
             'minimap2 environment configured in nextflow.config.')


def pairwise_coverage(contigs, min_id=90):
    """Return dict: (query_name, target_name) -> pct_coverage_of_query."""
    if len(contigs) < 2:
        return {}

    cov = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta = os.path.join(tmpdir, "all.fasta")
        with open(fasta, 'w') as f:
            for cd in contigs:
                f.write(f">{cd['name']}\n{cd['seq']}\n")

        paf = os.path.join(tmpdir, "self.paf")
        try:
            subprocess.run(
                ["minimap2", "-x", "asm5", "-c", "--secondary=no", "-D", fasta, fasta],
                stdout=open(paf, 'w'), stderr=subprocess.DEVNULL,
                timeout=120, check=False,
            )
        except FileNotFoundError:
            _minimap2_missing()
        except subprocess.TimeoutExpired:
            print('WARN: minimap2 containment check timed out; redundancy detection '
                  'skipped for this sample', file=sys.stderr)
            return {}

        covered = defaultdict(lambda: defaultdict(set))  # query -> target -> set of positions
        qlens = {}
        with open(paf) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                qname, qlen, qstart, qend = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
                tname = parts[5]
                nmatch, alen = int(parts[9]), int(parts[10])
                if qname == tname:
                    continue
                aln_id = (nmatch / alen * 100) if alen > 0 else 0
                if aln_id < min_id:
                    continue
                qlens[qname] = qlen
                for pos in range(qstart, qend):
                    covered[qname][tname].add(pos)

        for qname in covered:
            for tname in covered[qname]:
                if qlens.get(qname, 0) > 0:
                    cov[(qname, tname)] = len(covered[qname][tname]) / qlens[qname] * 100

    return cov


# Classification + filtering + grouping

def annotate_contigs(contigs, amr_by_contig, priority_patterns):
    for cd in contigs:
        cd['amr_genes'] = set()
        for entry in amr_by_contig.get(cd['name'], []):
            cd['amr_genes'].add(entry['gene'])
        m = re.search(r'resistance_genes=(\S+)', cd['extra'])
        if m and m.group(1) != '-':
            for g in m.group(1).replace(';', ',').split(','):
                g = g.strip()
                if g:
                    cd['amr_genes'].add(g)

        cd['has_amr'] = len(cd['amr_genes']) > 0
        cd['has_priority'] = False
        if priority_patterns and cd['amr_genes']:
            for gene in cd['amr_genes']:
                for pat in priority_patterns:
                    if fnmatch.fnmatch(gene, pat) or fnmatch.fnmatch(gene.lower(), pat.lower()):
                        cd['has_priority'] = True
                        break
                if cd['has_priority']:
                    break


def classify_contigs(contigs, min_length, circular_min, fragment_min):
    for cd in contigs:
        if cd['length'] < min_length:
            cd['class'] = 'artifact'
        elif cd['has_priority']:
            cd['class'] = 'priority_amr'
        elif cd['has_amr']:
            cd['class'] = 'amr_carrier'
        elif cd['is_circ'] and cd['length'] >= circular_min:
            cd['class'] = 'primary_replicon'
        elif cd['is_circ']:
            cd['class'] = 'small_circular'
        elif cd['length'] >= fragment_min:
            cd['class'] = 'significant_fragment'
        else:
            cd['class'] = 'small_fragment'


def filter_contigs(contigs, cov_map, containment_cov, circular_min, fragment_min, min_length,
                   single_mode=False):
    """Mark contigs as included/excluded, detect redundancy."""
    has_circular = any(cd['is_circ'] and cd['length'] >= circular_min for cd in contigs
                       if cd['class'] != 'artifact')

    for cd in contigs:
        cd['redundant'] = False
        cd['redundant_of'] = ''
        cd['include'] = False
        cd['reason'] = ''

    for cd in contigs:
        if cd['is_circ'] or cd['class'] == 'artifact':
            continue
        for circ in contigs:
            if not circ['is_circ'] or circ['class'] == 'artifact':
                continue
            pct = cov_map.get((cd['name'], circ['name']), 0)
            if pct >= containment_cov:
                unique_amr = cd['amr_genes'] - circ['amr_genes']
                if unique_amr:
                    cd['class'] = 'amr_carrier'
                    print(f"      {cd['name']}: {pct:.0f}% contained in {circ['name']} "
                          f"but unique AMR: {','.join(sorted(unique_amr))}, keeping", file=sys.stderr)
                else:
                    cd['redundant'] = True
                    cd['redundant_of'] = circ['name']
                break

    for cd in contigs:
        if cd['redundant']:
            cd['include'] = False
            cd['reason'] = f"redundant ({cd['redundant_of']})"
        elif cd['class'] == 'artifact':
            cd['include'] = False
            cd['reason'] = f"artifact (< {min_length} bp)"
        elif single_mode:
            cd['include'] = True
            cd['reason'] = 'MOB cluster member'
        elif cd['class'] in ('priority_amr', 'amr_carrier'):
            cd['include'] = True
            cd['reason'] = 'AMR carrier'
        elif cd['class'] in ('primary_replicon', 'small_circular'):
            cd['include'] = True
            cd['reason'] = 'circular replicon'
        elif cd['class'] == 'significant_fragment':
            if has_circular:
                cd['include'] = False
                cd['reason'] = 'linear fragment (circular exists)'
            else:
                cd['include'] = True
                cd['reason'] = 'largest linear (no circular)'
        elif cd['class'] == 'small_fragment':
            cd['include'] = False
            cd['reason'] = f"small fragment (< {fragment_min} bp)"
        else:
            cd['include'] = False
            cd['reason'] = 'unclassified'

    included = [cd for cd in contigs if cd['include']]
    if not included and contigs:
        longest = max(contigs, key=lambda c: c['length'])
        longest['include'] = True
        longest['reason'] = 'fallback (longest)'
        included = [longest]

    return included


def group_contigs(included, cov_map, dissimilar_threshold=50):
    """Group included contigs into separate comparison map groups."""
    circulars = [cd for cd in included if cd['is_circ']]
    linears = [cd for cd in included if not cd['is_circ']]

    if not circulars:
        for cd in included:
            cd['group'] = 0
        return [included]

    groups = []
    assigned = set()

    for i, circ in enumerate(circulars):
        group = [circ]
        assigned.add(circ['name'])
        for lin in linears:
            if lin['name'] in assigned:
                continue
            pct = cov_map.get((lin['name'], circ['name']), 0)
            if pct >= dissimilar_threshold:
                group.append(lin)
                assigned.add(lin['name'])
        groups.append(group)

    for lin in linears:
        if lin['name'] not in assigned:
            groups.append([lin])
            assigned.add(lin['name'])

    for gi, group in enumerate(groups):
        for cd in group:
            cd['group'] = gi

    return groups


def refine_group_references(groups, best_sample, all_sample_results, search_dir,
                            priority_patterns, circular_min, containment_id):
    """For each group, check if another sample has a better matching contig (more AMR)."""
    refined = []

    for gi, group in enumerate(groups):
        # Representative contig from initially selected sample
        rep = max(group, key=lambda c: c.get('score', 0))
        rep_len = rep['length']

        # Score the initial reference for this group
        best_group_score = sum(c.get('score', 0) for c in group)
        best_ref_sample = best_sample['name']
        best_ref_path = best_sample['path']
        best_ref_contigs = group

        print(f"  Group {gi}: representative={rep['name']} ({rep_len:,} bp), "
              f"initial ref={best_ref_sample} (score={best_group_score:.1f})", file=sys.stderr)

        # Check each other sample for a better matching contig
        for sr in all_sample_results:
            if sr['name'] == best_sample['name']:
                continue

            # Parse this sample's contigs and AMR
            other_contigs = parse_fasta_contigs(sr['path'])
            amrfinder_path = find_amrfinder(sr['name'], search_dir)
            other_amr = load_amrfinder(amrfinder_path) if amrfinder_path else defaultdict(list)

            annotate_contigs(other_contigs, other_amr, priority_patterns)

            # Find best matching contig in this sample via minimap2
            matching = _find_matching_contigs(rep, other_contigs, min_cov=50, min_id=containment_id)

            if not matching:
                continue

            # Score the matching contigs
            candidate_score = 0
            for mc in matching:
                candidate_score += score_contig(mc, circular_min)

            if candidate_score > best_group_score:
                print(f"    Better ref found: {sr['name']} (score={candidate_score:.1f} vs {best_group_score:.1f})",
                      file=sys.stderr)
                for mc in matching:
                    flags = []
                    if mc.get('is_circ'):
                        flags.append('circular')
                    if mc.get('has_amr'):
                        flags.append(f"AMR: {','.join(sorted(mc['amr_genes']))}")
                    print(f"      {mc['name']}: {mc['length']:,} bp [{'; '.join(flags)}]", file=sys.stderr)
                best_group_score = candidate_score
                best_ref_sample = sr['name']
                best_ref_path = sr['path']
                best_ref_contigs = matching

        # Build comp list for this group (all samples except the chosen ref)
        comp_paths = [sr['path'] for sr in all_sample_results if sr['name'] != best_ref_sample]

        refined.append({
            'contigs': best_ref_contigs,
            'ref_sample': best_ref_sample,
            'ref_path': best_ref_path,
            'comp_paths': comp_paths,
            'group_index': gi,
        })

    return refined


def _find_matching_contigs(representative, candidate_contigs, min_cov=50, min_id=90):
    """Find contigs in candidate that match the representative contig."""
    if not candidate_contigs:
        return []

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write representative as target
        ref_fa = os.path.join(tmpdir, "ref.fasta")
        with open(ref_fa, 'w') as f:
            f.write(f">{representative['name']}\n{representative['seq']}\n")

        # Write candidates as query
        qry_fa = os.path.join(tmpdir, "query.fasta")
        with open(qry_fa, 'w') as f:
            for cd in candidate_contigs:
                f.write(f">{cd['name']}\n{cd['seq']}\n")

        paf = os.path.join(tmpdir, "align.paf")
        try:
            subprocess.run(
                ["minimap2", "-x", "asm5", "-c", "--secondary=no", ref_fa, qry_fa],
                stdout=open(paf, 'w'), stderr=subprocess.DEVNULL,
                timeout=60, check=False,
            )
        except FileNotFoundError:
            _minimap2_missing()
        except subprocess.TimeoutExpired:
            print(f"WARN: minimap2 matching against {representative['name']} timed out",
                  file=sys.stderr)
            return []

        # Parse PAF: find candidate contigs that cover the representative
        best_by_query = {}
        with open(paf) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                qname = parts[0]
                qlen = int(parts[1])
                qstart, qend = int(parts[2]), int(parts[3])
                tlen = int(parts[6])
                nmatch, alen = int(parts[9]), int(parts[10])
                aln_id = (nmatch / alen * 100) if alen > 0 else 0
                if aln_id < min_id:
                    continue
                # Coverage of the representative (target)
                tstart, tend = int(parts[7]), int(parts[8])
                t_cov = (tend - tstart) / tlen * 100 if tlen > 0 else 0

                if qname not in best_by_query or t_cov > best_by_query[qname]:
                    best_by_query[qname] = t_cov

        # Return candidate contigs that cover >= min_cov of the representative
        matching = []
        for cd in candidate_contigs:
            cov = best_by_query.get(cd['name'], 0)
            if cov >= min_cov:
                matching.append(cd)

        return matching


def score_contig(cd, circular_min):
    score = 0
    if cd.get('has_priority'):
        score += 10000
    if cd.get('has_amr'):
        score += 1000
        # Bonus per AMR gene: favors contigs with more AMR genes
        score += len(cd.get('amr_genes', set())) * 50
    if cd['is_circ'] and cd['length'] >= circular_min:
        score += 500
        # Extra bonus for circular + AMR: strongly prefer circular contigs carrying AMR
        if cd.get('has_amr'):
            score += 5000
    elif cd['is_circ']:
        score += 100
    if cd['length'] > 0:
        score += math.log10(cd['length'])
    cd['score'] = score
    return score


# Group mode + AMR-coverage reference selection (single mode)

def decide_group_mode(requested, asm_type, has_closed_replicon):
    """Resolve --group-mode auto into 'split' or 'single'."""
    if requested in ('split', 'single'):
        return requested
    if (asm_type or '').strip().upper() == 'S':
        return 'single'
    if not has_closed_replicon:
        return 'single'
    return 'split'


def matches_priority(gene, patterns):
    for pat in patterns or []:
        if fnmatch.fnmatch(gene, pat) or fnmatch.fnmatch(gene.lower(), pat.lower()):
            return True
    return False


def sample_amr_genes(sample_result):
    """Distinct AMR genes carried by a sample's *included* contigs."""
    genes = set()
    for cd in sample_result['included']:
        genes |= cd.get('amr_genes', set())
    return genes


def compute_amr_union(sample_results):
    """Union of AMR genes across every sample in the cluster."""
    union = set()
    for sr in sample_results:
        union |= sample_amr_genes(sr)
    return union


def _contiguity(sample_result):
    """Return (n_contigs, n50) for the included contigs."""
    lengths = sorted((cd['length'] for cd in sample_result['included']), reverse=True)
    if not lengths:
        return (0, 0)
    total = sum(lengths)
    acc = 0
    n50 = lengths[-1]
    for L in lengths:
        acc += L
        if acc >= total / 2:
            n50 = L
            break
    return (len(lengths), n50)


def amr_rank_key(sample_result, union_genes, priority_patterns):
    """Lexicographic ranking key for single-mode reference selection."""
    genes = sample_amr_genes(sample_result)
    n_priority = sum(1 for g in genes if matches_priority(g, priority_patterns))
    n_amr = len(genes)
    frac = (len(genes & union_genes) / len(union_genes)) if union_genes else 0.0
    amr_len = sum(cd['length'] for cd in sample_result['included'] if cd.get('has_amr'))
    n_contigs, n50 = _contiguity(sample_result)
    return (-n_priority, -n_amr, -frac, -amr_len, n_contigs, -n50, sample_result['name'])


def rank_samples_by_amr(sample_results, union_genes, priority_patterns):
    return sorted(sample_results, key=lambda sr: amr_rank_key(sr, union_genes, priority_patterns))


def write_amr_matrix(output_dir, sample_results, union_genes, ref_genes, priority_patterns):
    """Companion outputs: per-sample AMR presence/absence + per-gene summary."""
    genes = sorted(union_genes)
    matrix_path = os.path.join(output_dir, "amr_presence_matrix.tsv")
    with open(matrix_path, 'w') as out:
        out.write("sample\t" + "\t".join(genes) + "\n" if genes else "sample\n")
        for sr in sample_results:
            present = sample_amr_genes(sr)
            row = ["1" if g in present else "0" for g in genes]
            out.write(sr['name'] + ("\t" + "\t".join(row) if genes else "") + "\n")

    summary_path = os.path.join(output_dir, "amr_gene_summary.tsv")
    with open(summary_path, 'w') as out:
        out.write("gene\tn_samples\tis_priority\ton_reference\trepresentable_on_map\n")
        for g in genes:
            n = sum(1 for sr in sample_results if g in sample_amr_genes(sr))
            on_ref = g in ref_genes
            out.write(f"{g}\t{n}\t{'yes' if matches_priority(g, priority_patterns) else 'no'}\t"
                      f"{'yes' if on_ref else 'no'}\t{'yes' if on_ref else 'no'}\n")
    return matrix_path, summary_path



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastas", nargs="+", required=True)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--priority-genes", default=None)
    parser.add_argument("--min-length", type=int, default=1000)
    parser.add_argument("--circular-min", type=int, default=20000)
    parser.add_argument("--fragment-min", type=int, default=5000)
    parser.add_argument("--containment-cov", type=float, default=80)
    parser.add_argument("--containment-id", type=float, default=90)
    parser.add_argument("--dissimilar-threshold", type=float, default=50)
    parser.add_argument("--group-mode", choices=["auto", "split", "single"], default="auto",
                        help="Contig grouping strategy. 'split' = one map per replicon group "
                             "(good for closed assemblies); 'single' = keep the whole MOB-suite "
                             "cluster as one reference (good for fragmented short-read bins); "
                             "'auto' = single for short-read (S) or when no closed replicon exists.")
    parser.add_argument("--asm-type", default="",
                        help="Assembly type (L/S/HL/HS/A), used by --group-mode auto")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    priority_patterns = []
    if args.priority_genes:
        priority_patterns = [p.strip() for p in args.priority_genes.split(',') if p.strip()]

    amr_by_plasmid = {}
    if args.metadata and os.path.exists(args.metadata):
        with open(args.metadata) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                name = row.get('plasmid_name', '')
                genes = row.get('resistance_genes', '')
                if genes and genes != '-':
                    amr_by_plasmid[name] = set(g.strip() for g in genes.replace(';', ',').split(',') if g.strip())

    search_dir = os.path.dirname(args.fastas[0]) or '.'

    # --- Pass 1: parse, annotate, classify, containment ---
    prepared = []

    for fasta_path in args.fastas:
        sample_name = os.path.splitext(os.path.basename(fasta_path))[0]
        contigs = parse_fasta_contigs(fasta_path)

        amrfinder_path = find_amrfinder(sample_name, search_dir)
        amr_by_contig = load_amrfinder(amrfinder_path) if amrfinder_path else defaultdict(list)

        plasmid_amr = amr_by_plasmid.get(sample_name, set())
        if plasmid_amr and not any(amr_by_contig.values()):
            for cd in contigs:
                for gene in plasmid_amr:
                    amr_by_contig[cd['name']].append({'gene': gene, 'pident': 0, 'coverage': 0})

        annotate_contigs(contigs, amr_by_contig, priority_patterns)
        classify_contigs(contigs, args.min_length, args.circular_min, args.fragment_min)

        # Pairwise coverage for containment + grouping
        non_artifact = [cd for cd in contigs if cd['class'] != 'artifact']
        cov_map = pairwise_coverage(non_artifact, min_id=args.containment_id) if len(non_artifact) >= 2 else {}

        prepared.append({
            'path': fasta_path, 'name': sample_name,
            'contigs': contigs, 'cov_map': cov_map,
        })

    # --- Decide grouping mode (needs cluster-wide circularity) ---
    has_closed = any(
        cd['is_circ'] and cd['length'] >= args.circular_min
        for p in prepared for cd in p['contigs'] if cd['class'] != 'artifact'
    )
    group_mode = decide_group_mode(args.group_mode, args.asm_type, has_closed)
    single_mode = (group_mode == 'single')
    print(f"  Group mode: {group_mode} (requested={args.group_mode}, "
          f"asm_type={args.asm_type or '?'}, closed_replicon={'yes' if has_closed else 'no'})",
          file=sys.stderr)

    # --- Pass 2: filter, group, score ---
    sample_results = []

    for p in prepared:
        contigs = p['contigs']
        cov_map = p['cov_map']

        included = filter_contigs(contigs, cov_map, args.containment_cov,
                                  args.circular_min, args.fragment_min, args.min_length,
                                  single_mode=single_mode)

        if single_mode:
            # The MOB-suite bin IS the unit: never re-partition it.
            for cd in included:
                cd['group'] = 0
            groups = [included] if included else []
        else:
            groups = group_contigs(included, cov_map, args.dissimilar_threshold)

        total_score = 0
        for cd in included:
            total_score += score_contig(cd, args.circular_min)

        sample_results.append({
            'score': total_score, 'path': p['path'], 'name': p['name'],
            'included': included, 'all_contigs': contigs, 'groups': groups,
        })

        print(f"  {p['name']}: score={total_score:.1f} ({len(contigs)} contigs, "
              f"{len(included)} included, {len(groups)} group(s))", file=sys.stderr)
        for cd in contigs:
            flags = []
            if cd.get('has_priority'):
                flags.append('PRIORITY')
            if cd.get('is_circ'):
                flags.append('circular')
            if cd.get('has_amr'):
                flags.append(f"AMR: {','.join(sorted(cd['amr_genes']))}")
            if cd.get('redundant'):
                flags.append(f"REDUNDANT of {cd['redundant_of']}")
            flag_str = f" [{'; '.join(flags)}]" if flags else ""
            status = "+" if cd.get('include') else "-"
            grp = f" grp={cd['group']}" if cd.get('include') and 'group' in cd else ""
            print(f"    {status} {cd['name']}: {cd['length']:,} bp, class={cd['class']}"
                  f"{flag_str}{grp} | {cd.get('reason', '')}", file=sys.stderr)

    # --- Reference selection ---
    amr_union = compute_amr_union(sample_results)
    selection = {}

    if single_mode:
        # Whole-bin reference: rank by AMR coverage, not by contig score.
        ranked = rank_samples_by_amr(sample_results, amr_union, priority_patterns)
        best = ranked[0]
        best_genes = sample_amr_genes(best)
        n_priority = sum(1 for g in best_genes if matches_priority(g, priority_patterns))
        selection = {
            'basis': 'amr_coverage',
            'priority_genes_matched': sorted(g for g in best_genes
                                             if matches_priority(g, priority_patterns)),
            'n_priority_genes': n_priority,
            'n_amr_genes': len(best_genes),
            'union_size': len(amr_union),
            'union_coverage': (len(best_genes & amr_union) / len(amr_union)) if amr_union else 0.0,
            'missing_from_reference': sorted(amr_union - best_genes),
        }
        print(f"\n  Reference by AMR coverage: {best['name']} "
              f"({n_priority} priority, {len(best_genes)}/{len(amr_union)} union genes)",
              file=sys.stderr)
        if selection['missing_from_reference']:
            print(f"    NOT representable on map: "
                  f"{', '.join(selection['missing_from_reference'])}", file=sys.stderr)

        refined_groups = [{
            'contigs': best['included'],
            'ref_sample': best['name'],
            'ref_path': best['path'],
            'comp_paths': [sr['path'] for sr in sample_results if sr['name'] != best['name']],
            'group_index': 0,
        }]
    else:
        # --- Pick best sample (initial selection for grouping) ---
        sample_results.sort(key=lambda x: x['score'], reverse=True)
        best = sample_results[0]
        selection = {'basis': 'contig_score', 'score': best['score']}

        # --- Refine: per-group reference selection ---
        refined_groups = refine_group_references(
            best['groups'], best, sample_results, search_dir,
            priority_patterns, args.circular_min, args.containment_id,
        )

    # --- Write per-group outputs ---
    group_labels = []
    alphabet = 'abcdefghijklmnopqrstuvwxyz'

    for rg in refined_groups:
        gi = rg['group_index']
        suffix = alphabet[gi] if len(refined_groups) > 1 else ''
        group = rg['contigs']

        # Sort by score within group
        group.sort(key=lambda x: x.get('score', 0), reverse=True)

        # Check if group has AMR genes
        circ_in_group = [cd for cd in group if cd.get('is_circ')]
        amr_in_group = [cd for cd in group if cd.get('has_amr')]
        total_len = sum(cd['length'] for cd in group)
        n_amr = sum(len(cd.get('amr_genes', set())) for cd in amr_in_group)
        group_has_amr = n_amr > 0

        ref_path = os.path.join(args.output_dir, f"ref_{suffix}.fasta" if suffix else "ref.fasta")
        with open(ref_path, 'w') as out:
            for i, cd in enumerate(group, 1):
                header = f">contig_{i} length={cd['length']}"
                if cd.get('is_circ'):
                    header += " circular=true"
                if cd.get('has_amr'):
                    header += f" amr_genes={';'.join(sorted(cd['amr_genes']))}"
                header += f" original_name={cd['name']}"
                cleaned = cd.get('extra', '')
                cleaned = re.sub(r'\s*length=\d+', '', cleaned)
                cleaned = re.sub(r'\s*circular=[Tt]rue', '', cleaned)
                cleaned = re.sub(r'\s*rotated=[Tt]rue', '', cleaned)
                cleaned = re.sub(r'\s*resistance_genes=\S+', '', cleaned)
                cleaned = cleaned.strip()
                if cleaned:
                    header += f" {cleaned}"
                out.write(header + '\n')
                for j in range(0, len(cd['seq']), 80):
                    out.write(cd['seq'][j:j + 80] + '\n')

        ref_sample_name = rg['ref_sample']
        details = []
        if circ_in_group:
            details.append("circular")
        else:
            details.append("linear")
        if n_amr > 0:
            details.append(f"{n_amr} AMR genes")
        detail_str = f" ({', '.join(details)})" if details else ""
        display_name = f"{ref_sample_name}{detail_str}"

        # Collect all AMR gene names for this group
        all_amr_genes = set()
        for cd in amr_in_group:
            all_amr_genes.update(cd.get('amr_genes', set()))

        group_labels.append({
            'suffix': suffix,
            'ref_fasta': os.path.basename(ref_path),
            'display_name': display_name,
            'has_amr': group_has_amr,
            'amr_gene_names': sorted(all_amr_genes),
            'num_contigs': len(group),
            'total_length': total_len,
            'circular': len(circ_in_group),
            'amr': len(amr_in_group),
            'contig_names': [cd['name'] for cd in group],
            'ref_sample': ref_sample_name,
        })

    # Comp list (same for all groups - full sample FASTAs, exclude best initial sample)
    comp_files = [s['path'] for s in sample_results if s['path'] != best['path']]
    comp_path = os.path.join(args.output_dir, "comp_list.txt")
    with open(comp_path, 'w') as out:
        out.write(' '.join(comp_files))

    # AMRFinderPlus per-group (different groups may have different reference samples)
    for gl in group_labels:
        af_path = find_amrfinder(gl['ref_sample'], search_dir)
        gl['amrfinder_tsv'] = af_path or ''

    # Info file (uses initial best for backward compat)
    amrfinder_path = find_amrfinder(best['name'], search_dir)
    info_path = os.path.join(args.output_dir, "ref_info.txt")
    all_ref_contigs = []
    for rg in refined_groups:
        all_ref_contigs.extend(rg['contigs'])
    circ_all = [cd for cd in all_ref_contigs if cd.get('is_circ')]
    amr_all = [cd for cd in all_ref_contigs if cd.get('has_amr')]
    with open(info_path, 'w') as out:
        out.write(f"reference_sample={best['name']}\n")
        out.write(f"reference_score={best['score']:.1f}\n")
        out.write(f"group_mode={group_mode}\n")
        out.write(f"selection_basis={selection.get('basis', '')}\n")
        if single_mode:
            out.write(f"priority_genes_matched={';'.join(selection.get('priority_genes_matched', []))}\n")
            out.write(f"amr_union_size={selection.get('union_size', 0)}\n")
            out.write(f"amr_union_coverage={selection.get('union_coverage', 0.0):.3f}\n")
            missing = selection.get('missing_from_reference', [])
            out.write(f"amr_missing_from_reference={';'.join(missing)}\n")
        out.write(f"num_groups={len(refined_groups)}\n")
        out.write(f"num_contigs={len(all_ref_contigs)}\n")
        out.write(f"total_length={sum(cd['length'] for cd in all_ref_contigs)}\n")
        out.write(f"circular_contigs={len(circ_all)}\n")
        out.write(f"amr_contigs={len(amr_all)}\n")
        if amr_all:
            all_amr = set()
            for cd in amr_all:
                all_amr.update(cd.get('amr_genes', set()))
            out.write(f"amr_genes={';'.join(sorted(all_amr))}\n")
        if amrfinder_path:
            out.write(f"amrfinder_tsv={amrfinder_path}\n")

    ref_genes = set()
    for rg in refined_groups:
        for cd in rg['contigs']:
            ref_genes |= cd.get('amr_genes', set())

    manifest = {
        'reference_sample': best['name'],
        'score': best['score'],
        'group_mode': group_mode,
        'asm_type': args.asm_type,
        'selection': selection,
        'amr_union': sorted(amr_union),
        'amr_on_reference': sorted(ref_genes),
        'amr_missing_from_reference': sorted(amr_union - ref_genes),
        'comp_files': [os.path.basename(f) for f in comp_files],
        'groups': group_labels,
    }
    manifest_path = os.path.join(args.output_dir, "groups_manifest.json")
    with open(manifest_path, 'w') as out:
        json.dump(manifest, out, indent=2)

    # Companion outputs: AMR presence/absence matrix + per-gene summary.
    write_amr_matrix(args.output_dir, sample_results, amr_union, ref_genes, priority_patterns)

    print(f"\nInitial best sample: {best['name']} (score={best['score']:.1f})", file=sys.stderr)
    print(f"  {len(refined_groups)} group(s):", file=sys.stderr)
    for gl in group_labels:
        ref_note = f" [ref: {gl['ref_sample']}]" if gl['ref_sample'] != best['name'] else ""
        print(f"  Group {gl['suffix'] or '-'}: {gl['num_contigs']} contigs, "
              f"{gl['total_length']:,} bp, {gl['display_name']}{ref_note}", file=sys.stderr)
    print(f"  Comparisons: {len(comp_files)} samples", file=sys.stderr)


if __name__ == '__main__':
    main()