#!/usr/bin/env python3
"""Load ReporTree partitions and SKA2 distances into surveillance database."""
import argparse
import csv
import json
import os
import sys
from datetime import datetime

from sqlalchemy import (
    create_engine, MetaData, Table, Column, Integer, BigInteger, Float,
    String, Text, DateTime, Boolean, JSON, Index, inspect, text
)
from sqlalchemy.exc import SQLAlchemyError

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib.utils import add_db_arguments, resolve_db, mask_url

metadata_obj = MetaData()

# ── Schema definition (database-agnostic) ──

surveillance_species = Table('surveillance_species', metadata_obj,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('species_name', Text, unique=True, nullable=False),
    Column('sketch_path', Text),
    Column('created_at', DateTime, default=datetime.utcnow),
    Column('updated_at', DateTime, default=datetime.utcnow),
)

surveillance_isolates = Table('surveillance_isolates', metadata_obj,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('isolate_id', Text, nullable=False),
    Column('sample', Text, nullable=False),
    Column('species_id', Integer),
    Column('species_name', Text, nullable=False),
    Column('workflow_id', Text, nullable=False),
    Column('workflow_description', Text),
    Column('group_id', Text),
    Column('group_description', Text),
    Column('assembly_type', Text),
    Column('mlst_scheme', Text),
    Column('mlst_st', Text),
    Column('sketch_path', Text),
    Column('fasta_path', Text),
    Column('collection_date', DateTime),
    Column('metadata_json', JSON),
    Column('created_at', DateTime, default=datetime.utcnow),
    Index('idx_isolates_species', 'species_name'),
    Index('idx_isolates_workflow', 'workflow_id'),
)

surveillance_pairwise = Table('surveillance_pairwise_distances', metadata_obj,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('species_name', Text, nullable=False),
    Column('isolate_1', Text, nullable=False),
    Column('isolate_2', Text, nullable=False),
    Column('snp_distance', Integer, nullable=False),
    Column('run_id', Text, nullable=False),
    Column('created_at', DateTime, default=datetime.utcnow),
    Index('idx_pairwise_species', 'species_name'),
    Index('idx_pairwise_dist', 'snp_distance'),
)

surveillance_clusters = Table('surveillance_clusters', metadata_obj,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('species_name', Text, nullable=False),
    Column('cluster_name', Text, nullable=False),
    Column('nomenclature_code', Text),
    Column('threshold_snp', Integer, nullable=False),
    Column('member_count', Integer, default=0),
    Column('first_seen', DateTime),
    Column('last_updated', DateTime, default=datetime.utcnow),
    Column('is_active', Boolean, default=True),
    Index('idx_clusters_species', 'species_name'),
    Index('idx_clusters_threshold', 'threshold_snp'),
)

surveillance_members = Table('surveillance_cluster_members', metadata_obj,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('cluster_id', Integer),
    Column('isolate_id', Text, nullable=False),
    Column('species_name', Text, nullable=False),
    Column('threshold_snp', Integer, nullable=False),
    Column('nomenclature_code', Text),
    Column('assigned_at', DateTime, default=datetime.utcnow),
    Column('run_id', Text, nullable=False),
    Index('idx_members_cluster', 'cluster_id'),
    Index('idx_members_isolate', 'isolate_id'),
)

surveillance_history = Table('surveillance_cluster_history', metadata_obj,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('species_name', Text, nullable=False),
    Column('cluster_name', Text, nullable=False),
    Column('threshold_snp', Integer, nullable=False),
    Column('event_type', Text, nullable=False),
    Column('description', Text),
    Column('old_name', Text),
    Column('new_name', Text),
    Column('affected_isolates', JSON),  # stored as JSON array
    Column('run_id', Text, nullable=False),
    Column('event_time', DateTime, default=datetime.utcnow),
    Index('idx_history_species', 'species_name'),
    Index('idx_history_event', 'event_type'),
)

surveillance_runs = Table('surveillance_reportree_runs', metadata_obj,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('run_id', Text, unique=True, nullable=False),
    Column('species_name', Text, nullable=False),
    Column('workflow_id', Text, nullable=False),
    Column('n_isolates', Integer),
    Column('n_new_isolates', Integer),
    Column('thresholds', Text, default='10,20,50,100'),
    Column('partitions_path', Text),
    Column('status', Text, default='completed'),
    Column('started_at', DateTime, default=datetime.utcnow),
    Column('completed_at', DateTime),
)

# ── Denormalized summary tables (populated per-run, database-agnostic) ──

surveillance_sample_overview = Table('surveillance_sample_overview', metadata_obj,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('isolate_id', Text, nullable=False),
    Column('sample', Text, nullable=False),
    Column('species_name', Text, nullable=False),
    Column('assembly_type', Text),
    Column('nomenclature_code', Text),
    Column('cluster_hc100', Text),
    Column('cluster_hc50', Text),
    Column('cluster_hc20', Text),
    Column('cluster_hc10', Text),
    Column('in_outbreak_cluster', Boolean, default=False),
    Column('tightest_cluster_level', Text),
    Column('mlst_scheme', Text),
    Column('mlst_st', Text),
    Column('workflow_id', Text),
    Column('workflow_description', Text),
    Column('group_id', Text),
    Column('group_description', Text),
    Column('first_seen', DateTime),
    Index('idx_overview_species', 'species_name'),
    Index('idx_overview_isolate', 'isolate_id'),
    Index('idx_overview_sample', 'sample'),
    Index('idx_overview_code', 'nomenclature_code'),
)

surveillance_cluster_summary = Table('surveillance_cluster_summary', metadata_obj,
    Column('id', Integer, primary_key=True, autoincrement=True),
    Column('species_name', Text, nullable=False),
    Column('cluster_name', Text, nullable=False),
    Column('nomenclature_code', Text),
    Column('threshold_snp', Integer, nullable=False),
    Column('member_count', Integer),
    Column('first_seen', DateTime),
    Column('last_updated', DateTime),
    Column('is_active', Boolean, default=True),
    Column('isolate_ids', Text),  # comma-separated list (works on all DBs)
    Column('new_members_30d', Integer, default=0),
    Column('new_members_7d', Integer, default=0),
    Index('idx_summary_species', 'species_name'),
    Index('idx_summary_threshold', 'threshold_snp'),
)


def init_schema(engine):
    """Create all tables if they don't exist."""
    metadata_obj.create_all(engine)
    print("  Schema initialized (tables created if not present)")


def upsert_row(conn, table, match_cols, data):
    """Delete-then-insert upsert that works on both SQLite and PostgreSQL."""
    where_clause = None
    for col_name in match_cols:
        val = data.get(col_name)
        if val is None:
            continue
        cond = (table.c[col_name] == val)
        where_clause = cond if where_clause is None else (where_clause & cond)

    if where_clause is not None:
        conn.execute(table.delete().where(where_clause))
    conn.execute(table.insert().values(**data))


def ensure_species(conn, species_name, sketch_path=None):
    """Insert species if not exists, return species id."""
    result = conn.execute(
        surveillance_species.select().where(surveillance_species.c.species_name == species_name)
    ).fetchone()

    if result:
        if sketch_path:
            conn.execute(
                surveillance_species.update()
                .where(surveillance_species.c.species_name == species_name)
                .values(sketch_path=sketch_path, updated_at=datetime.utcnow())
            )
        return result[0]
    else:
        result = conn.execute(
            surveillance_species.insert().values(
                species_name=species_name, sketch_path=sketch_path,
                created_at=datetime.utcnow(), updated_at=datetime.utcnow()
            )
        )
        return result.inserted_primary_key[0]


def load_isolates(conn, isolates_tsv, species_name, species_id, workflow_id, workflow_description, group_id, group_description, run_id):
    """Load isolate metadata."""
    if not isolates_tsv or not os.path.exists(isolates_tsv):
        return 0

    count = 0
    with open(isolates_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            isolate_id = row.get('isolate_id', row.get('sample', ''))
            sample = row.get('sample', isolate_id)
            asm_type = row.get('assembly_type', '')
            if not isolate_id:
                continue

            upsert_row(conn, surveillance_isolates,
                ['sample'],
                {
                    'isolate_id': isolate_id,
                    'sample': sample,
                    'species_id': species_id,
                    'species_name': species_name,
                    'workflow_id': workflow_id,
                    'workflow_description': workflow_description,
                    'group_id': group_id,
                    'group_description': group_description,
                    'assembly_type': asm_type,
                    'mlst_scheme': row.get('mlst_scheme', ''),
                    'mlst_st': row.get('mlst_st', ''),
                    'sketch_path': row.get('sketch_path', ''),
                    'fasta_path': row.get('fasta_path', ''),
                    'created_at': datetime.utcnow(),
                })
            count += 1
    return count


def load_distances(conn, distances_tsv, species_name, run_id):
    """Load SKA2 pairwise distances."""
    if not distances_tsv or not os.path.exists(distances_tsv):
        return 0

    count = 0
    with open(distances_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        # Normalize headers: SKA2 outputs 'Sample 1', 'Sample 2', 'SNP distance'.
        headers = reader.fieldnames or []
        header_map = {h.lower().replace(' ', '_'): h for h in headers}
        print(f"  Distance file headers: {headers}")

        # Find the right column names
        iso1_col = header_map.get('sample_1') or header_map.get('sample1') or header_map.get('isolate_1')
        iso2_col = header_map.get('sample_2') or header_map.get('sample2') or header_map.get('isolate_2')
        dist_col = (header_map.get('snp_distance') or header_map.get('snps')
                    or header_map.get('distance') or header_map.get('snpdistance'))

        if not iso1_col or not iso2_col or not dist_col:
            print(f"  WARNING: Could not find required columns. Found: {headers}")
            print(f"  Looking for sample_1/sample_2/snp_distance (case-insensitive)")
            return 0

        print(f"  Using columns: {iso1_col}, {iso2_col}, {dist_col}")

        for row in reader:
            iso1 = row.get(iso1_col, '')
            iso2 = row.get(iso2_col, '')
            dist = row.get(dist_col, '')

            if not iso1 or not iso2 or dist == '':
                continue

            if iso1 == iso2:
                continue

            try:
                snp_dist = int(float(dist))
            except (ValueError, TypeError):
                continue

            if iso1 > iso2:
                iso1, iso2 = iso2, iso1

            upsert_row(conn, surveillance_pairwise,
                ['species_name', 'isolate_1', 'isolate_2'],
                {
                    'species_name': species_name,
                    'isolate_1': iso1,
                    'isolate_2': iso2,
                    'snp_distance': snp_dist,
                    'run_id': run_id,
                    'created_at': datetime.utcnow(),
                })
            count += 1
    print(f"  Distance file: {distances_tsv} ({count} pairs loaded, self-pairs skipped)")
    return count


def _cluster_from_distances(conn, species_name, run_id, sample_names, thresholds):
    """Build clusters from pairwise distances when ReporTree doesn't produce partition columns."""
    import re

    # Load all pairwise distances for this species from the database
    dist_rows = conn.execute(
        surveillance_pairwise.select().where(
            surveillance_pairwise.c.species_name == species_name
        )
    ).fetchall()

    distances = {}
    for r in dist_rows:
        distances[(r.isolate_1, r.isolate_2)] = r.snp_distance
        distances[(r.isolate_2, r.isolate_1)] = r.snp_distance

    cluster_count = 0
    member_count = 0

    # For each threshold, do single-linkage clustering
    isolate_assignments = {s: {} for s in sample_names}

    for threshold in sorted(thresholds):
        # Build adjacency: samples within threshold distance
        adj = {s: set() for s in sample_names}
        for i, s1 in enumerate(sample_names):
            for s2 in sample_names[i+1:]:
                d = distances.get((s1, s2))
                if d is not None and d <= threshold:
                    adj[s1].add(s2)
                    adj[s2].add(s1)

        # Find connected components (single-linkage clusters)
        visited = set()
        cluster_num = 0
        for s in sorted(sample_names):
            if s in visited:
                continue
            # BFS to find all connected samples
            component = []
            queue = [s]
            while queue:
                node = queue.pop(0)
                if node in visited:
                    continue
                visited.add(node)
                component.append(node)
                for neighbor in sorted(adj[node]):
                    if neighbor not in visited:
                        queue.append(neighbor)

            if len(component) > 1:
                cluster_num += 1
                cluster_name = f'cluster_{cluster_num}'
                for sample in component:
                    isolate_assignments[sample][threshold] = cluster_name
            else:
                singleton_name = f'singleton_{s}'
                isolate_assignments[s][threshold] = singleton_name

    # Propagate upward (fine → coarse)
    for sample in sample_names:
        last_cluster = ''
        for t in sorted(thresholds):
            val = isolate_assignments[sample].get(t, '')
            if val and not val.startswith('singleton'):
                last_cluster = val
            elif not val and last_cluster:
                isolate_assignments[sample][t] = last_cluster

    def cluster_to_code(cluster_name):
        if not cluster_name or cluster_name.startswith('singleton'):
            return 'C0'
        m = re.match(r'^cluster_(\d+)$', cluster_name)
        if m:
            return f'C{m.group(1)}'
        return cluster_name

    def build_nomenclature_code(sample):
        parts = []
        for t in sorted(thresholds, reverse=True):
            cn = isolate_assignments.get(sample, {}).get(t, '')
            parts.append(cluster_to_code(cn))
        return '-'.join(parts)

    # Insert clusters and members
    for threshold in sorted(thresholds):
        cluster_members_map = {}
        for sample in sample_names:
            cn = isolate_assignments[sample].get(threshold, '')
            if not cn:
                continue
            if cn not in cluster_members_map:
                cluster_members_map[cn] = []
            cluster_members_map[cn].append(sample)

        for cluster_name, members in cluster_members_map.items():
            nom_code = build_nomenclature_code(members[0])
            is_multi = len(members) > 1

            upsert_row(conn, surveillance_clusters,
                ['species_name', 'cluster_name', 'threshold_snp'],
                {
                    'species_name': species_name,
                    'cluster_name': cluster_name,
                    'threshold_snp': threshold,
                    'nomenclature_code': nom_code,
                    'member_count': len(members),
                    'first_seen': datetime.utcnow(),
                    'last_updated': datetime.utcnow(),
                    'is_active': is_multi,
                })
            cluster_count += 1

            result = conn.execute(
                surveillance_clusters.select().where(
                    (surveillance_clusters.c.species_name == species_name) &
                    (surveillance_clusters.c.cluster_name == cluster_name) &
                    (surveillance_clusters.c.threshold_snp == threshold)
                )
            ).fetchone()
            cluster_id = result[0] if result else None

            for sample in members:
                iso_code = build_nomenclature_code(sample)
                upsert_row(conn, surveillance_members,
                    ['species_name', 'isolate_id', 'threshold_snp'],
                    {
                        'cluster_id': cluster_id,
                        'isolate_id': sample,
                        'species_name': species_name,
                        'threshold_snp': threshold,
                        'nomenclature_code': iso_code,
                        'run_id': run_id,
                        'assigned_at': datetime.utcnow(),
                    })
                member_count += 1

    print(f"  Distance-based clustering: {cluster_count} clusters, {member_count} memberships")
    return cluster_count, member_count


def load_partitions(conn, partitions_tsv, species_name, run_id):
    """Load ReporTree partitions into clusters and cluster_members."""
    if not partitions_tsv or not os.path.exists(partitions_tsv):
        return 0, 0

    import re

    cluster_count = 0
    member_count = 0
    thresholds = [100, 50, 20, 10]

    with open(partitions_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames

        # Detect threshold columns: support multiple ReporTree output formats.
        threshold_cols = {}
        col_threshold_map = {}
        for h in headers:
            # single-Nx1.0 format (HC output)
            m = re.match(r'^single-(\d+)x[\d.]+$', h)
            if m:
                col_threshold_map[int(m.group(1))] = h
                continue
            m = re.match(r'^HC-single-(\d+)$', h, re.IGNORECASE)
            if m:
                col_threshold_map[int(m.group(1))] = h
                continue
            m = re.match(r'^HC(\d+)$', h, re.IGNORECASE)
            if m:
                col_threshold_map[int(m.group(1))] = h
                continue

        for t in thresholds:
            if t in col_threshold_map:
                threshold_cols[t] = col_threshold_map[t]

        print(f"  Total columns: {len(headers)}")
        print(f"  Detected threshold columns: {threshold_cols}")
        if not threshold_cols:
            print(f"  First 10 headers: {headers[:10]}")

        rows = list(reader)

    if not threshold_cols:
        print(f"  WARNING: No threshold columns found in {partitions_tsv}")
        # ReporTree may omit partition columns when all samples cluster together, or the cached
        # partitions may only have a subset of samples.
        sample_set = set()
        for row in rows:
            name = row.get('sequence', row.get('sample', row.get(headers[0], '')))
            if name:
                sample_set.add(name)

        # Also add samples from pairwise distances (covers historical samples in merged sketch)
        for r in conn.execute(
            surveillance_pairwise.select().where(
                surveillance_pairwise.c.species_name == species_name
            )
        ):
            sample_set.add(r.isolate_1)
            sample_set.add(r.isolate_2)

        sample_names = sorted(sample_set)
        if len(sample_names) > 1:
            print(f"  Falling back to distance-based clustering for {len(sample_names)} samples")
            cluster_count, member_count = _cluster_from_distances(
                conn, species_name, run_id, sample_names, thresholds
            )
            return cluster_count, member_count
        return 0, 0

    # --- Build nomenclature codes: C{N} format --- C0 = singleton (not clustered), C1 =
    # cluster_1, C2 = cluster_2, etc. Code is 4-part: C{100}-C{50}-C{20}-C{10}.
    all_isolates = [row.get('sequence', row.get('sample', row.get(headers[0], '')))
                    for row in rows]
    all_isolates = [i for i in all_isolates if i]
    missing = [t for t in thresholds if t not in threshold_cols]
    if missing and len(all_isolates) > 1:
        print(f"  Thresholds not computed by ReporTree (above max pairwise distance): "
              f"{sorted(missing)}")
        print(f"  -> assigning all {len(all_isolates)} isolate(s) to a single cluster "
              f"at each; every pair is closer than the threshold.")
    elif missing:
        print(f"  Thresholds not computed by ReporTree: {sorted(missing)} "
              f"(single isolate, so singleton at every level)")

    isolate_assignments = {}  # isolate -> {threshold: cluster_name}
    for row in rows:
        isolate = row.get('sequence', row.get('sample', row.get(headers[0], '')))
        if not isolate:
            continue
        isolate_assignments[isolate] = {}
        for t in thresholds:
            col = threshold_cols.get(t)
            if col:
                isolate_assignments[isolate][t] = row.get(col, '')
            elif len(all_isolates) > 1:
                isolate_assignments[isolate][t] = 'cluster_1'
            else:
                isolate_assignments[isolate][t] = f'singleton_{isolate}'

    # Propagate cluster assignments upward: if clustered at HC10, must also be clustered at HC20,
    # HC50, HC100.
    for isolate in isolate_assignments:
        last_cluster = ''
        for t in sorted(thresholds):  # 10, 20, 50, 100
            val = isolate_assignments[isolate][t]
            if val and not val.startswith('singleton'):
                last_cluster = val
            elif not val and last_cluster:
                # Missing column at coarser threshold: propagate finer cluster up
                isolate_assignments[isolate][t] = last_cluster

    def cluster_to_code(cluster_name):
        """Convert cluster_name to C{N} format."""
        if not cluster_name or cluster_name.startswith('singleton'):
            return 'C0'
        m = re.match(r'^cluster_(\d+)$', cluster_name)
        if m:
            return f'C{m.group(1)}'
        return cluster_name  # fallback

    def build_nomenclature_code(isolate):
        """Build 4-part C{N} nomenclature code for an isolate (HC100-HC50-HC20-HC10)."""
        parts = []
        for t in sorted(thresholds, reverse=True):  # 100, 50, 20, 10
            cn = isolate_assignments.get(isolate, {}).get(t, '')
            parts.append(cluster_to_code(cn))
        return '-'.join(parts)

    # --- Load each threshold --- Ensure we load all 4 thresholds even if some columns were
    # missing (propagated assignments fill the gaps)
    for threshold in sorted(thresholds):
        cluster_members_map = {}

        col_name = threshold_cols.get(threshold)
        for isolate, assignments in isolate_assignments.items():
            cluster_name = assignments.get(threshold, '')
            if not cluster_name:
                continue
            if cluster_name not in cluster_members_map:
                cluster_members_map[cluster_name] = []
            cluster_members_map[cluster_name].append(isolate)

        for cluster_name, members in cluster_members_map.items():
            nom_code = build_nomenclature_code(members[0])

            upsert_row(conn, surveillance_clusters,
                ['species_name', 'cluster_name', 'threshold_snp'],
                {
                    'species_name': species_name,
                    'cluster_name': cluster_name,
                    'threshold_snp': threshold,
                    'nomenclature_code': nom_code,
                    'member_count': len(members),
                    'first_seen': datetime.utcnow(),
                    'last_updated': datetime.utcnow(),
                    'is_active': len(members) > 1,
                })
            cluster_count += 1

            result = conn.execute(
                surveillance_clusters.select().where(
                    (surveillance_clusters.c.species_name == species_name) &
                    (surveillance_clusters.c.cluster_name == cluster_name) &
                    (surveillance_clusters.c.threshold_snp == threshold)
                )
            ).fetchone()
            cluster_id = result[0] if result else None

            for isolate in members:
                iso_code = build_nomenclature_code(isolate)
                upsert_row(conn, surveillance_members,
                    ['species_name', 'isolate_id', 'threshold_snp'],
                    {
                        'cluster_id': cluster_id,
                        'isolate_id': isolate,
                        'species_name': species_name,
                        'threshold_snp': threshold,
                        'nomenclature_code': iso_code,
                        'run_id': run_id,
                        'assigned_at': datetime.utcnow(),
                    })
                member_count += 1

    # Check if there are samples in pairwise distances not covered by partitions.
    partitioned_samples = set(isolate_assignments.keys())
    distance_samples = set()
    for r in conn.execute(
        surveillance_pairwise.select().where(
            surveillance_pairwise.c.species_name == species_name
        )
    ):
        distance_samples.add(r.isolate_1)
        distance_samples.add(r.isolate_2)

    missing_samples = distance_samples - partitioned_samples
    if missing_samples:
        print(f"  {len(missing_samples)} samples in distances but not in partitions; re-clustering all from distances")
        # Wipe what we just loaded and redo from distances for consistency
        conn.execute(surveillance_members.delete().where(
            surveillance_members.c.species_name == species_name))
        conn.execute(surveillance_clusters.delete().where(
            surveillance_clusters.c.species_name == species_name))
        all_samples = sorted(partitioned_samples | distance_samples)
        cluster_count, member_count = _cluster_from_distances(
            conn, species_name, run_id, all_samples, thresholds
        )

    return cluster_count, member_count


def load_nomenclature_changes(conn, changes_tsv, species_name, run_id):
    """Load nomenclature changes into cluster_history."""
    if not changes_tsv or not os.path.exists(changes_tsv):
        return 0

    count = 0
    with open(changes_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            event_type = row.get('event', row.get('type', 'unknown'))
            old_name = row.get('old_cluster', row.get('previous', ''))
            new_name = row.get('new_cluster', row.get('current', ''))
            threshold = row.get('threshold', '')
            isolates_str = row.get('isolates', row.get('samples', ''))
            description = row.get('description', row.get('details', ''))

            try:
                threshold_snp = int(threshold) if threshold else 0
            except (ValueError, TypeError):
                threshold_snp = 0

            affected = [x.strip() for x in isolates_str.split(',') if x.strip()] if isolates_str else []

            conn.execute(surveillance_history.insert().values(
                species_name=species_name,
                cluster_name=new_name or old_name,
                threshold_snp=threshold_snp,
                event_type=event_type,
                description=description,
                old_name=old_name,
                new_name=new_name,
                affected_isolates=affected,
                run_id=run_id,
                event_time=datetime.utcnow(),
            ))
            count += 1
    return count


def rebuild_summary_tables(conn, species_name):
    """Rebuild denormalized summary tables for this species."""
    now = datetime.utcnow()
    thirty_days_ago = datetime(now.year, now.month, now.day) 
    seven_days_ago = datetime(now.year, now.month, now.day)
    # Approximate: just use timedelta
    from datetime import timedelta
    thirty_days_ago = now - timedelta(days=30)
    seven_days_ago = now - timedelta(days=7)

    conn.execute(
        surveillance_sample_overview.delete().where(
            surveillance_sample_overview.c.species_name == species_name
        )
    )

    # Get all isolates at HC10 (base threshold)
    members_10 = conn.execute(
        surveillance_members.select().where(
            (surveillance_members.c.species_name == species_name) &
            (surveillance_members.c.threshold_snp == 10)
        )
    ).fetchall()

    for m10 in members_10:
        sample_name = m10.isolate_id  # This is the sample name with _L/_S from SKA
        nom_code = m10.nomenclature_code

        # Derive clean isolate_id and assembly_type from sample name
        parts = sample_name.rsplit('_', 1)
        if len(parts) == 2 and parts[1] in ('L', 'S', 'HL', 'HS'):
            clean_isolate_id = parts[0]
            asm_type = parts[1]
        else:
            clean_isolate_id = sample_name
            asm_type = ''

        # Get cluster assignments at each threshold
        cluster_at = {}
        member_count_at = {}
        for t in [10, 20, 50, 100]:
            cm_row = conn.execute(
                surveillance_members.select().where(
                    (surveillance_members.c.species_name == species_name) &
                    (surveillance_members.c.isolate_id == sample_name) &
                    (surveillance_members.c.threshold_snp == t)
                )
            ).fetchone()
            if cm_row and cm_row.cluster_id:
                c_row = conn.execute(
                    surveillance_clusters.select().where(
                        surveillance_clusters.c.id == cm_row.cluster_id
                    )
                ).fetchone()
                cluster_at[t] = c_row.cluster_name if c_row else 'Singleton'
                member_count_at[t] = c_row.member_count if c_row else 1
            else:
                cluster_at[t] = 'Singleton'
                member_count_at[t] = 1

        # Determine tightest cluster level in_outbreak_cluster = TRUE if clustered at HC20 or tighter
        # (≤20 SNPs)
        if member_count_at.get(10, 1) > 1:
            tightest = 'HC10'
            in_outbreak = True
        elif member_count_at.get(20, 1) > 1:
            tightest = 'HC20'
            in_outbreak = True
        elif member_count_at.get(50, 1) > 1:
            tightest = 'HC50'
            in_outbreak = False
        elif member_count_at.get(100, 1) > 1:
            tightest = 'HC100'
            in_outbreak = False
        else:
            tightest = 'Singleton'
            in_outbreak = False

        # Get isolate metadata from surveillance_isolates (matched by sample name)
        iso_row = conn.execute(
            surveillance_isolates.select().where(
                (surveillance_isolates.c.species_name == species_name) &
                (surveillance_isolates.c.sample == sample_name)
            )
        ).fetchone()

        conn.execute(surveillance_sample_overview.insert().values(
            isolate_id=iso_row.isolate_id if iso_row else clean_isolate_id,
            sample=iso_row.sample if iso_row else sample_name,
            species_name=species_name,
            assembly_type=iso_row.assembly_type if iso_row else asm_type,
            nomenclature_code=nom_code,
            cluster_hc100=cluster_at[100],
            cluster_hc50=cluster_at[50],
            cluster_hc20=cluster_at[20],
            cluster_hc10=cluster_at[10],
            in_outbreak_cluster=in_outbreak,
            tightest_cluster_level=tightest,
            mlst_scheme=iso_row.mlst_scheme if iso_row else '',
            mlst_st=iso_row.mlst_st if iso_row else '',
            workflow_id=iso_row.workflow_id if iso_row else '',
            workflow_description=iso_row.workflow_description if iso_row else '',
            group_id=iso_row.group_id if iso_row else '',
            group_description=iso_row.group_description if iso_row else '',
            first_seen=iso_row.created_at if iso_row else now,
        ))

    n_overview = len(members_10)

    conn.execute(
        surveillance_cluster_summary.delete().where(
            surveillance_cluster_summary.c.species_name == species_name
        )
    )

    clusters = conn.execute(
        surveillance_clusters.select().where(
            (surveillance_clusters.c.species_name == species_name) &
            (surveillance_clusters.c.is_active == True) &
            (surveillance_clusters.c.member_count > 1)
        )
    ).fetchall()

    for c in clusters:
        # Get member isolate IDs
        members = conn.execute(
            surveillance_members.select().where(
                surveillance_members.c.cluster_id == c.id
            )
        ).fetchall()

        isolate_ids = ','.join(sorted(set(m.isolate_id for m in members)))

        new_30d = sum(1 for m in members if m.assigned_at and m.assigned_at > thirty_days_ago)
        new_7d = sum(1 for m in members if m.assigned_at and m.assigned_at > seven_days_ago)

        conn.execute(surveillance_cluster_summary.insert().values(
            species_name=c.species_name,
            cluster_name=c.cluster_name,
            nomenclature_code=c.nomenclature_code,
            threshold_snp=c.threshold_snp,
            member_count=c.member_count,
            first_seen=c.first_seen,
            last_updated=c.last_updated,
            is_active=c.is_active,
            isolate_ids=isolate_ids,
            new_members_30d=new_30d,
            new_members_7d=new_7d,
        ))

    n_summary = len(clusters)
    print(f"  Summary tables rebuilt: {n_overview} sample_overview rows, {n_summary} cluster_summary rows")


def main():
    parser = argparse.ArgumentParser(description='Load surveillance results into database')
    add_db_arguments(parser)
    parser.add_argument('--species', required=True, help='Species name')
    parser.add_argument('--workflow-id', required=True, help='Workflow ID')
    parser.add_argument('--workflow-description', default='', help='Workflow description')
    parser.add_argument('--group-id', default='', help='Group ID')
    parser.add_argument('--group-description', default='', help='Group description')
    parser.add_argument('--run-id', required=True, help='Unique run identifier')
    parser.add_argument('--isolates-tsv', default=None, help='Isolate metadata TSV')
    parser.add_argument('--distances-tsv', default=None, help='SKA2 pairwise distances TSV')
    parser.add_argument('--partitions-tsv', default=None, help='ReporTree partitions.tsv')
    parser.add_argument('--nomenclature-changes', default=None, help='ReporTree nomenclature_changes.tsv')
    parser.add_argument('--merged-sketch-path', default=None, help='Path to merged species .skf')
    parser.add_argument('--partitions-output', default=None, help='Copy partitions.tsv here for next run')
    args = parser.parse_args()
    args.db_url = resolve_db(args)

    print(f"Loading surveillance data for {args.species} (run: {args.run_id})")
    print(f"  Database: {mask_url(args.db_url)}")

    engine = create_engine(args.db_url, echo=False)

    # Initialize schema via SQLAlchemy ORM (works for both SQLite and PostgreSQL)
    init_schema(engine)

    with engine.begin() as conn:
        species_id = ensure_species(conn, args.species, args.merged_sketch_path)
        print(f"  Species: {args.species} (id={species_id})")

        # Clean clustering data for this species before loading fresh results.
        for tbl in [surveillance_sample_overview, surveillance_cluster_summary,
                     surveillance_members, surveillance_clusters, surveillance_pairwise]:
            deleted = conn.execute(
                tbl.delete().where(tbl.c.species_name == args.species)
            )
            print(f"  Cleared {deleted.rowcount} rows from {tbl.name} for {args.species}")

        n_isolates = load_isolates(conn, args.isolates_tsv, args.species, species_id,
                                    args.workflow_id, args.workflow_description,
                                    args.group_id, args.group_description, args.run_id)
        print(f"  Isolates: {n_isolates} upserted from current run")

        n_distances = load_distances(conn, args.distances_tsv, args.species, args.run_id)
        print(f"  Distances: {n_distances} upserted")

        n_clusters, n_members = load_partitions(conn, args.partitions_tsv, args.species, args.run_id)
        print(f"  Clusters: {n_clusters} upserted, {n_members} memberships")

        # Backfill isolate records for samples in partitions/distances but not in isolates table.
        existing_isolates = set()
        for r in conn.execute(
            surveillance_isolates.select().where(
                surveillance_isolates.c.species_name == args.species
            )
        ):
            existing_isolates.add(r.sample)

        # Collect all sample names from cluster_members (which came from partitions)
        all_clustered = set()
        for r in conn.execute(
            surveillance_members.select().where(
                surveillance_members.c.species_name == args.species
            )
        ):
            all_clustered.add(r.isolate_id)

        missing_isolates = all_clustered - existing_isolates
        if missing_isolates:
            print(f"  Backfilling {len(missing_isolates)} isolate records (historical samples from merged sketch)")
            for sample_name in sorted(missing_isolates):
                # Double-check it truly doesn't exist (race condition safety)
                already = conn.execute(
                    surveillance_isolates.select().where(
                        surveillance_isolates.c.sample == sample_name
                    )
                ).fetchone()
                if already:
                    continue

                # Parse isolate_id by stripping assembly type suffix
                parts = sample_name.rsplit('_', 1)
                if len(parts) == 2 and parts[1] in ('L', 'S', 'HL', 'HS'):
                    iso_id = parts[0]
                    asm_type = parts[1]
                else:
                    iso_id = sample_name
                    asm_type = ''

                # Try to find metadata from ov_general_statistics (if table exists)
                gs_row = None
                try:
                    gs_row = conn.execute(
                        text("SELECT workflow_id, workflow_description, group_id, group_description, "
                             "mlst_scheme, mlst_st FROM ov_general_statistics WHERE sample = :s LIMIT 1"),
                        {'s': sample_name}
                    ).fetchone()
                except Exception:
                    pass

                # Fallback: look up another isolate with the same isolate_id in this species
                if not gs_row:
                    existing = conn.execute(
                        surveillance_isolates.select().where(
                            (surveillance_isolates.c.isolate_id == iso_id) &
                            (surveillance_isolates.c.species_name == args.species)
                        )
                    ).fetchone()
                else:
                    existing = None

                wf_id = (gs_row[0] if gs_row else
                         existing.workflow_id if existing else args.workflow_id)
                wf_desc = (gs_row[1] if gs_row else
                           existing.workflow_description if existing else args.workflow_description)
                grp_id = (gs_row[2] if gs_row else
                          existing.group_id if existing else args.group_id)
                grp_desc = (gs_row[3] if gs_row else
                            existing.group_description if existing else args.group_description)
                mlst_sch = (gs_row[4] if gs_row else
                            existing.mlst_scheme if existing else '')
                mlst_s = (gs_row[5] if gs_row else
                          existing.mlst_st if existing else '')

                conn.execute(surveillance_isolates.insert().values(
                    isolate_id=iso_id,
                    sample=sample_name,
                    species_id=species_id,
                    species_name=args.species,
                    workflow_id=wf_id,
                    workflow_description=wf_desc,
                    group_id=grp_id,
                    group_description=grp_desc,
                    assembly_type=asm_type,
                    mlst_scheme=mlst_sch,
                    mlst_st=mlst_s,
                    created_at=datetime.utcnow(),
                ))
            n_isolates += len(missing_isolates)

        # Backfill singletons: isolates in this run but not in any cluster_member (e.g. single-isolate
        # species where ReporTree had no data)
        if n_isolates > 0 and n_members == 0:
            print(f"  No partitions loaded; inserting all {n_isolates} isolates as singletons")
            isolate_rows = conn.execute(
                surveillance_isolates.select().where(
                    (surveillance_isolates.c.species_name == args.species) &
                    (surveillance_isolates.c.workflow_id == args.workflow_id)
                )
            ).fetchall()

            for idx, iso_row in enumerate(isolate_rows, start=1):
                sample_name = iso_row[2]  # 'sample' column
                for threshold in [100, 50, 20, 10]:
                    singleton_name = f'singleton_{idx}'
                    nom_code = 'C0-C0-C0-C0'

                    upsert_row(conn, surveillance_clusters,
                        ['species_name', 'cluster_name', 'threshold_snp'],
                        {
                            'species_name': args.species,
                            'cluster_name': singleton_name,
                            'threshold_snp': threshold,
                            'nomenclature_code': nom_code,
                            'member_count': 1,
                            'first_seen': datetime.utcnow(),
                            'last_updated': datetime.utcnow(),
                            'is_active': False,
                        })

                    result = conn.execute(
                        surveillance_clusters.select().where(
                            (surveillance_clusters.c.species_name == args.species) &
                            (surveillance_clusters.c.cluster_name == singleton_name) &
                            (surveillance_clusters.c.threshold_snp == threshold)
                        )
                    ).fetchone()
                    cluster_id = result[0] if result else None

                    upsert_row(conn, surveillance_members,
                        ['species_name', 'isolate_id', 'threshold_snp'],
                        {
                            'cluster_id': cluster_id,
                            'isolate_id': sample_name,
                            'species_name': args.species,
                            'threshold_snp': threshold,
                            'nomenclature_code': nom_code,
                            'run_id': args.run_id,
                            'assigned_at': datetime.utcnow(),
                        })
                    n_members += 1
            print(f"  Singleton backfill: {n_members} memberships added")

        n_changes = load_nomenclature_changes(conn, args.nomenclature_changes, args.species, args.run_id)
        print(f"  History events: {n_changes}")

        upsert_row(conn, surveillance_runs, ['run_id'], {
            'run_id': args.run_id,
            'species_name': args.species,
            'workflow_id': args.workflow_id,
            'n_isolates': n_isolates,
            'n_new_isolates': n_isolates,
            'partitions_path': args.partitions_tsv,
            'status': 'completed',
            'started_at': datetime.utcnow(),
            'completed_at': datetime.utcnow(),
        })

        rebuild_summary_tables(conn, args.species)

    # Copy partitions for next run
    if args.partitions_tsv and args.partitions_output and os.path.exists(args.partitions_tsv):
        import shutil
        shutil.copy2(args.partitions_tsv, args.partitions_output)
        print(f"  Partitions saved to: {args.partitions_output}")

    engine.dispose()
    print("Done")


if __name__ == '__main__':
    main()