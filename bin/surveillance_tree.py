#!/usr/bin/env python3
"""Render one cumulative interactive MST per species from SKA distances."""

import argparse
import csv
import heapq
import json
import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from epiquence_lib.surveillance_scope import make_filter, scope_label

TEMPLATE = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'templates', 'surveillance_tree.html')


def now_in(tz_name):
    """Current time in the configured IANA zone, e.g. America/Chicago."""
    if tz_name:
        try:
            from zoneinfo import ZoneInfo
            return datetime.now(ZoneInfo(tz_name))
        except Exception as exc:
            print(f'WARN: timezone {tz_name!r} unavailable ({exc}); '
                  f'falling back to system local time', file=sys.stderr)
    return datetime.now()


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True


def resolve_columns(headers):
    """Locate sample/sample/distance columns case-insensitively."""
    hm = {h.lower().replace(' ', '_').replace('(', '').replace(')', ''): h
          for h in headers}
    s1 = hm.get('sample1') or hm.get('sample_1') or hm.get('isolate_1')
    s2 = hm.get('sample2') or hm.get('sample_2') or hm.get('isolate_2')
    d = (hm.get('distance') or hm.get('snp_distance') or hm.get('snps')
         or hm.get('snpdistance'))
    return s1, s2, d


def stream_pairs(path, keep=None):
    """Yield (sample1, sample2, distance) from a ska distance TSV."""
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


def first_pass(path, cap, max_pairs, keep=None):
    """Index samples and retain the closest pairs."""
    index = {}
    heap = []          # (-dist, i, j): max-heap by distance
    pruned = False

    for a, b, d in stream_pairs(path, keep):
        if a not in index:
            index[a] = len(index)
        if b not in index:
            index[b] = len(index)
        if d > cap:
            continue
        i, j = index[a], index[b]
        if i > j:
            i, j = j, i
        if len(heap) < max_pairs:
            heapq.heappush(heap, (-d, i, j))
        elif -heap[0][0] > d:
            heapq.heapreplace(heap, (-d, i, j))
            pruned = True

    retained = [(-nd, i, j) for nd, i, j in heap]

    if not retained:
        return index, [], cap

    if not pruned:
        # Everything at or below the requested cap is present, whether or not any pair actually
        # reaches it.
        return index, retained, cap

    # Pruning can split a tie: some pairs at the boundary distance kept, others dropped.
    boundary = max(d for d, _, _ in retained)
    retained = [t for t in retained if t[0] < boundary]
    eff_cap = boundary - 1
    if not retained:
        return index, [], 0

    return index, retained, eff_cap


def second_pass(path, index, uf, keep=None):
    """Minimum-weight edge between each pair of distinct components."""
    bridges = {}
    for a, b, d in stream_pairs(path, keep):
        i, j = index.get(a), index.get(b)
        if i is None or j is None:
            continue
        ci, cj = uf.find(i), uf.find(j)
        if ci == cj:
            continue
        key = (ci, cj) if ci < cj else (cj, ci)
        cur = bridges.get(key)
        if cur is None or d < cur[0]:
            bridges[key] = (d, i, j)
    return bridges


def build_mst(path, index, retained, keep=None):
    """Exact global MST."""
    n = len(index)
    uf = UnionFind(n)
    edges = []

    # Stage 1: all retained (cheap) edges, ascending.
    for d, i, j in sorted(retained):
        if uf.union(i, j):
            edges.append((i, j, d))

    # Stage 2: bridge the resulting components.
    n_components = len({uf.find(i) for i in range(n)})
    if n_components > 1:
        bridges = second_pass(path, index, uf, keep)
        for (_, _), (d, i, j) in sorted(bridges.items(), key=lambda kv: kv[1][0]):
            if uf.union(i, j):
                edges.append((i, j, d))

    return edges, n_components


def load_isolates(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline='', encoding='utf-8') as f:
        return {(r.get('sample') or '').strip(): r
                for r in csv.DictReader(f, delimiter='\t')
                if (r.get('sample') or '').strip()}


def load_partitions(path):
    """ReporTree partitions.tsv -> {sample: {column: value}}."""
    if not path or not os.path.exists(path) or path == 'NO_PARTITIONS':
        return {}, []
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        try:
            header = next(reader)
        except StopIteration:
            return {}, []
        if len(header) < 2:
            return {}, []
        cols = header[1:]
        out = {}
        for row in reader:
            if not row:
                continue
            out[row[0].strip()] = {
                c: (row[k + 1] if k + 1 < len(row) else '') for k, c in enumerate(cols)
            }
    return out, cols


def main():
    p = argparse.ArgumentParser(description='Render cumulative surveillance MST for one species')
    p.add_argument('--distances', required=True, help='ska distance TSV (long format)')
    p.add_argument('--isolates', required=True, help='Cumulative isolates TSV')
    p.add_argument('--partitions', default=None, help='ReporTree partitions.tsv')
    p.add_argument('--species', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--snp-cap', type=int, default=100)
    p.add_argument('--max-pairs', type=int, default=500000)
    p.add_argument('--height', type=int, default=900)
    p.add_argument('--timezone', default=None,
                   help='IANA timezone for the generated timestamp (e.g. America/Chicago)')
    p.add_argument('--group', default=None,
                   help='Restrict to isolates whose cumulative group_id equals this '
                        '(an empty string selects isolates without a group). '
                        'Omit for the species-wide view.')
    args = p.parse_args()

    if not os.path.exists(args.distances):
        print(f'[{args.species}] no distances file; skipping tree', file=sys.stderr)
        return

    isolates = load_isolates(args.isolates)
    keep = make_filter(isolates, args.group)
    label = scope_label(isolates, args.group)
    tag = args.species + (f' / group {args.group or "(none)"}' if args.group is not None else '')

    assert_no_duplicate_pairs(args.distances, keep=keep)
    index, retained, eff_cap = first_pass(args.distances, args.snp_cap, args.max_pairs, keep)
    n = len(index)

    if n < 2:
        print(f'[{tag}] {n} isolate(s) in scope; skipping tree', file=sys.stderr)
        return

    if eff_cap < args.snp_cap:
        print(f'[{args.species}] WARN: pairs within {args.snp_cap} exceeded '
              f'--max-pairs ({args.max_pairs}); cap lowered to {eff_cap}. '
              f'This population is highly clonal.', file=sys.stderr)

    mst_edges, n_components = build_mst(args.distances, index, retained, keep)

    # A spanning tree over n nodes has exactly n-1 edges.
    if len(mst_edges) != n - 1:
        sys.exit(f'ERROR: MST has {len(mst_edges)} edges for {n} isolates '
                 f'(expected {n - 1}); distances file may be incomplete')

    partitions, part_cols = load_partitions(args.partitions)

    rev = [''] * n
    for name, i in index.items():
        rev[i] = name

    meta_cols = ['Cluster (' + c + ')' for c in part_cols] + [
        'MLST scheme', 'AMR gene count', 'First seen', 'Last run', 'Group', 'Assembly type',
    ]

    samples = []
    runs = set()
    for i, name in enumerate(rev):
        iso = isolates.get(name, {})
        part = partitions.get(name, {})
        first_seen = iso.get('first_seen_workflow', '')
        if first_seen:
            runs.add(first_seen)
        mf = {'Cluster (' + c + ')': part.get(c, '') for c in part_cols}
        mf.update({
            'MLST scheme': iso.get('mlst_scheme', ''),
            'AMR gene count': iso.get('chromosomal_amr_gene_count', ''),
            'First seen': first_seen,
            'Last run': iso.get('workflow_id', ''),
            'Group': iso.get('group_id', ''),
            'Assembly type': iso.get('assembly_type', ''),
        })
        samples.append({
            'id': i,
            'name': name,
            'st': iso.get('mlst_st', ''),
            'amr': iso.get('chromosomal_amr_genes', ''),
            'metaFields': mf,
        })

    # Only MLST and AMR are offered as colour dimensions for now.
    color_fields = ['MLST']
    if any(s['amr'] for s in samples):
        color_fields.append('Chromosomal AMR')

    pairs_flat = []
    for d, i, j in retained:
        pairs_flat.extend((i, j, d))

    edges_json = [{'source': i, 'target': j, 'snps': d} for i, j, d in mst_edges]

    species_safe = args.species.replace(' ', '_')
    init_threshold = min(20, eff_cap) if eff_cap > 0 else 0
    # Only surfaced in the page when the cap actually had to be lowered, so a silently reduced
    # collapse range never goes unnoticed.
    cap_warning = ''
    if eff_cap < args.snp_cap:
        cap_warning = (
            f'\n    <div class="epiq-warn">Collapse range reduced from {args.snp_cap} to '
            f'{eff_cap}: this population is closely related enough that carrying every '
            f'pair within {args.snp_cap} would have exceeded the size budget.</div>')

    with open(TEMPLATE, encoding='utf-8') as f:
        html = f.read()

    for token, value in (
        ('@@EPIQ_div_id@@', 'surv_' + species_safe.replace('-', '_').replace('.', '_')),
        ('@@EPIQ_cluster_name@@', species_safe + '_surveillance'),
        # Scope (group) sits on its own line under the species, styled like the distance matrix page,
        # and joins the browser title with pipes.
        ('@@EPIQ_species@@', args.species),
        ('@@EPIQ_title_scope@@', f' | {label}' if label else ''),
        ('@@EPIQ_scope_html@@', f'\n    <div class="epiq-scope">{label}</div>' if label else ''),
        ('@@EPIQ_h@@', args.height),
        ('@@EPIQ_n_isolates@@', n),
        ('@@EPIQ_n_runs@@', len(runs) if runs else 'n/a'),
        ('@@EPIQ_snp_cap@@', eff_cap),
        ('@@EPIQ_init_threshold@@', init_threshold),
        ('@@EPIQ_generated@@', now_in(args.timezone).strftime('%Y-%m-%d %H:%M %Z').strip()),
        ('@@EPIQ_cap_warning@@', cap_warning),
        ('@@EPIQ_samples_json@@', json.dumps(samples)),
        ('@@EPIQ_edges_json@@', json.dumps(edges_json)),
        ('@@EPIQ_pairs_json@@', json.dumps(pairs_flat)),
        ('@@EPIQ_color_fields_json@@', json.dumps(color_fields)),
    ):
        html = html.replace(token, str(value))

    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(html)

    size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(f'[{tag}] tree: {n} isolates, {n_components} lineage(s) at '
          f'<={eff_cap} distance, {len(retained)} pairs shipped, {size_mb:.1f} MB',
          file=sys.stderr)


if __name__ == '__main__':
    main()