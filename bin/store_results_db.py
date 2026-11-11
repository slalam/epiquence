#!/usr/bin/env python3
"""Store EPIquence pipeline summary results into a database."""
import argparse
import csv
import json
import os
import re
import sys
from datetime import datetime

from sqlalchemy import (
    create_engine, MetaData, Table, Column, Integer, BigInteger, Float, String, Text,
    DateTime, JSON, inspect, text, select, func
)
from sqlalchemy.exc import SQLAlchemyError

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib.utils import add_db_arguments, resolve_db, mask_url


def _add_missing_columns(engine, table_name, columns):
    """ALTER in any of `columns` that the existing table lacks."""
    inspector = inspect(engine)
    if not inspector.has_table(table_name):
        return
    existing = {c['name'] for c in inspector.get_columns(table_name)}
    for col, sql_type in columns.items():
        if col not in existing:
            with engine.begin() as conn:
                conn.execute(text(f'ALTER TABLE {table_name} ADD COLUMN {col} {sql_type}'))
            print(f'  Added column {table_name}.{col}')


def local_now():
    """Current time in the task's zone."""
    return datetime.now()


def parse_nf_datetime(value):
    """Parse Nextflow's workflow.start / workflow.complete string."""
    if not value:
        return None
    v = str(value).strip()
    if not v:
        return None
    try:
        dt = datetime.fromisoformat(v)
        return dt.replace(tzinfo=None) if dt.tzinfo is None else dt.astimezone().replace(tzinfo=None)
    except ValueError:
        pass
    for fmt in ('%Y-%m-%d %H:%M:%S', '%Y-%m-%dT%H:%M:%S'):
        try:
            return datetime.strptime(v, fmt)
        except ValueError:
            continue
    # "Mon Aug 03 15:22:31 CDT 2026": strptime cannot parse a named zone abbreviation with %Z, so
    # drop the token and read the rest as local.
    m = re.match(r'^(\w{3} \w{3} +\d+ \d{2}:\d{2}:\d{2}) +[A-Z]{2,5} +(\d{4})$', v)
    if m:
        try:
            return datetime.strptime(f'{m.group(1)} {m.group(2)}', '%a %b %d %H:%M:%S %Y')
        except ValueError:
            pass
    print(f'WARN: could not parse workflow start {value!r}', file=sys.stderr)
    return None


# Column type mappings for each table.
JSON_ARRAY_COLUMNS = {
    'amr_gene_list', 'chromosomal_amr_gene_list', 'plasmid_amr_gene_list',
    'plasmid_cluster_list', 'resistance_genes', 'rep_type', 'relaxase_type',
    'orit_type', 'point_mutation_list', 'resistant_antimicrobial_list',
    'amr_plus_gene_list', 'amrfinder_gene_list', 'resfinder_gene_list',
    'amrfinder_point_mutation_list',
}

BIGINT_COLUMNS = {
    'total_reads', 'total_bases', 'filtered_reads', 'filtered_bases',
    'assembly_length_bp', 'n50_bp', 'largest_contig_bp', 'size',
}

FLOAT_COLUMNS = {
    'median_phred', 'mean_coverage', 'median_coverage', 'avg_read_len',
    'gc_pct', 'pct_reads_passed', 'pct_bases_passed',
    'pct_ge_20x', 'pct_ge_30x', 'pct_ge_40x', 'pct_ge_50x', 'pct_ge_100x',
    'ani', 'af',
    'pct_coverage_of_reference', 'pct_identity_to_reference',
}

INTEGER_COLUMNS = {
    'amr_gene_count', 'plasmid_cluster_count',
    'total_contigs', 'plasmid_contigs',
    'chromosomal_contigs', 'l50',
    'num_resistance_genes',
    'point_mutation_count', 'resistant_antimicrobial_count',
    'amr_plus_gene_count',
    'match', 'start', 'stop', 'target_length', 'reference_sequence_length',
    'alignment_length', 'pmid',
}


# Columns to skip entirely (not stored in database)
SKIP_COLUMNS = {
    'sequence',
}


def parse_value(col_name, raw_value):
    """Convert a raw TSV value to the appropriate Python type."""
    if raw_value is None or raw_value == '':
        return None

    if col_name in JSON_ARRAY_COLUMNS:
        # Split on semicolons or commas and store as Python list SQLAlchemy JSON type handles
        # serialization
        items = []
        for delim in [';', ',']:
            if delim in raw_value:
                items = [x.strip() for x in raw_value.split(delim) if x.strip()]
                break
        if not items and raw_value.strip() and raw_value.strip() != '-':
            items = [raw_value.strip()]
        return items

    if col_name in FLOAT_COLUMNS:
        try:
            return float(raw_value)
        except (ValueError, TypeError):
            return None

    if col_name in INTEGER_COLUMNS:
        try:
            return int(float(raw_value))
        except (ValueError, TypeError):
            return None

    if col_name in BIGINT_COLUMNS:
        try:
            return int(float(raw_value))
        except (ValueError, TypeError):
            return None

    if col_name == 'workflow_timestamp' and raw_value:
        try:
            return datetime.strptime(raw_value, '%Y-%m-%d %H:%M:%S')
        except (ValueError, TypeError):
            return raw_value

    return raw_value if raw_value != '-' else None


def get_column_type(col_name):
    """Return SQLAlchemy column type for a given column name."""
    if col_name in JSON_ARRAY_COLUMNS:
        return JSON
    if col_name in FLOAT_COLUMNS:
        return Float
    if col_name in BIGINT_COLUMNS:
        return BigInteger
    if col_name in INTEGER_COLUMNS:
        return Integer
    if col_name == 'workflow_timestamp':
        return DateTime
    return Text()


def load_tsv(filepath):
    """Load a TSV file and return (headers, rows)."""
    if not filepath or not os.path.exists(filepath):
        return [], []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames or []
        rows = list(reader)
    return headers, rows


def ensure_table(engine, metadata, table_name, columns):
    """Create table if it doesn't exist."""
    inspector = inspect(engine)

    if not inspector.has_table(table_name):
        cols = [Column('id', Integer, primary_key=True, autoincrement=True)]
        for col_name in columns:
            col_type = get_column_type(col_name)
            cols.append(Column(col_name, col_type, nullable=True))
        table = Table(table_name, metadata, *cols)
        metadata.create_all(engine, tables=[table])
        print(f"  Created table: {table_name} ({len(columns)} columns)")
        return table
    else:
        # Table exists: check for missing columns and add them.
        existing_cols = {c['name'].lower() for c in inspector.get_columns(table_name)}
        table = Table(table_name, metadata, autoload_with=engine)

        with engine.begin() as conn:
            for col_name in columns:
                if col_name.lower() not in existing_cols:
                    col_type = get_column_type(col_name)
                    type_str = 'TEXT'
                    if col_type == Float:
                        type_str = 'REAL'
                    elif col_type == BigInteger:
                        type_str = 'BIGINT'
                    elif col_type == Integer:
                        type_str = 'INTEGER'
                    elif col_type == JSON:
                        type_str = 'JSON'
                    elif col_type == DateTime:
                        type_str = 'TIMESTAMP'
                    conn.execute(text(
                        f'ALTER TABLE {table_name} ADD COLUMN "{col_name}" {type_str}'))
                    existing_cols.add(col_name.lower())
                    print(f"  Added column {col_name} ({type_str}) to {table_name}")

        # Reload table metadata after alterations
        metadata.clear()
        table = Table(table_name, metadata, autoload_with=engine)
        return table


def insert_rows(engine, table, headers, rows):
    """Insert rows into a table with upsert behavior."""
    if not rows:
        print(f"  No rows to insert into {table.name}")
        return 0

    records = []
    for row in rows:
        record = {}
        for col_name in headers:
            raw = row.get(col_name, '')
            record[col_name] = parse_value(col_name, raw)
        records.append(record)

    has_workflow_id = 'workflow_id' in headers or 'Workflow_ID' in headers
    has_sample = 'sample' in headers
    wid_col_name = 'workflow_id' if 'workflow_id' in headers else 'Workflow_ID' if 'Workflow_ID' in headers else None
    sample_col_name = 'sample' if 'sample' in headers else None

    with engine.begin() as conn:
        if wid_col_name and hasattr(table.c, wid_col_name):
            wid_column = getattr(table.c, wid_col_name)

            # Check for pairwise tables with isolate_1, isolate_2
            has_isolate_1 = 'isolate_1' in headers and hasattr(table.c, 'isolate_1')
            has_isolate_2 = 'isolate_2' in headers and hasattr(table.c, 'isolate_2')

            if has_isolate_1 and has_isolate_2:
                # Pairwise tables: delete by (workflow_id, isolate_1, isolate_2)
                triples = set()
                for rec in records:
                    wid = rec.get(wid_col_name)
                    iso1 = rec.get('isolate_1')
                    iso2 = rec.get('isolate_2')
                    if wid and iso1 and iso2 and (wid, iso1, iso2) not in triples:
                        triples.add((wid, iso1, iso2))
                        conn.execute(
                            table.delete().where(
                                (wid_column == wid) &
                                (table.c.isolate_1 == iso1) &
                                (table.c.isolate_2 == iso2)
                            )
                        )
                if triples:
                    print(f"  Cleared {len(triples)} existing (workflow_id, isolate_1, isolate_2) triples in {table.name}")

            elif sample_col_name and hasattr(table.c, sample_col_name):
                sample_column = getattr(table.c, sample_col_name)
                # Delete existing rows for each (workflow_id, sample) pair
                pairs = set()
                for rec in records:
                    wid = rec.get(wid_col_name)
                    smp = rec.get(sample_col_name)
                    if wid and smp and (wid, smp) not in pairs:
                        pairs.add((wid, smp))
                        conn.execute(
                            table.delete().where(
                                (wid_column == wid) & (sample_column == smp)
                            )
                        )
                if pairs:
                    print(f"  Cleared {len(pairs)} existing (workflow_id, sample) pairs in {table.name}")
            else:
                # Fallback: delete all rows for this workflow_id
                workflow_ids = set(rec.get(wid_col_name) for rec in records if rec.get(wid_col_name))
                for wid in workflow_ids:
                    conn.execute(table.delete().where(wid_column == wid))
                if workflow_ids:
                    print(f"  Cleared existing rows for workflow(s) {', '.join(workflow_ids)} in {table.name}")

        conn.execute(table.insert(), records)

    print(f"  Inserted {len(records)} rows into {table.name}")
    return len(records)


def canonicalize_headers(headers, rows, table_name):
    """Lower-case every column name and resolve names that then collide."""
    canon, seen, renames = [], {}, []
    for h in headers:
        base = h.strip().lower()
        if base not in seen:
            seen[base] = 1
            canon.append(base)
        else:
            seen[base] += 1
            alt = f'{base}_{seen[base]}'
            while alt in seen:
                seen[base] += 1
                alt = f'{base}_{seen[base]}'
            seen[alt] = 1
            canon.append(alt)
            renames.append((h, alt))

    if renames:
        pretty = ', '.join(f'{o} -> {n}' for o, n in renames)
        print(f"  {table_name}: renamed duplicate column(s) after case folding: {pretty}")

    if canon != headers:
        remapped = []
        for row in rows:
            remapped.append({new: row.get(old, '') for old, new in zip(headers, canon)})
        rows = remapped
    return canon, rows


def store_tsv(engine, metadata, table_name, filepath):
    """Load a TSV and store it in the database."""
    headers, rows = load_tsv(filepath)
    if not headers:
        print(f"  Skipping {table_name}: file not found or empty ({filepath})")
        return 0

    # Filter out columns we don't want to store
    skip_lower = {c.lower() for c in SKIP_COLUMNS}
    keep = [h for h in headers if h.lower() not in skip_lower]
    headers, rows = canonicalize_headers(keep, rows, table_name)

    table = ensure_table(engine, metadata, table_name, headers)
    return insert_rows(engine, table, headers, rows)


def rebuild_meta_tables(engine, metadata):
    """Rebuild metadata lookup tables from general_statistics."""
    inspector = inspect(engine)

    meta_tables = {
        'meta_species': {
            'value_col': 'species_name',
            'source_col': 'species',
        },
        'meta_isolates': {
            'value_col': 'isolate_id',
            'source_col': 'isolate_id',
        },
        'meta_amr_genes': {
            'value_col': 'amr_gene',
            'source_col': 'amr_gene_list',  # JSON array: needs unnesting
        },
        'meta_point_mutations': {
            'value_col': 'point_mutation',
            'source_col': 'point_mutation_list',  # JSON array: needs unnesting
        },
        'meta_plasmid_clusters': {
            'value_col': 'plasmid_cluster',
            'source_col': 'plasmid_cluster_list',  # JSON array: needs unnesting
        },
        'meta_virulence_factors': {
            'value_col': 'virulence_factor',
            'source_col': 'virulence_factor_list',  # JSON array: needs unnesting
        },
    }

    for table_name, config in meta_tables.items():
        value_col = config['value_col']

        # Create table if it doesn't exist
        if not inspector.has_table(table_name):
            Table(table_name, metadata,
                Column('id', Integer, primary_key=True, autoincrement=True),
                Column(value_col, Text, unique=True, nullable=False),
                Column('first_seen', DateTime, default=local_now),
            )
            metadata.create_all(engine, tables=[metadata.tables[table_name]])
            print(f"  Created metadata table: {table_name}")
        elif table_name not in metadata.tables:
            Table(table_name, metadata, autoload_with=engine)

    # Collect distinct values from ov_general_statistics
    gs_table_name = 'ov_general_statistics'
    if not inspector.has_table(gs_table_name):
        print("  Skipping meta tables: general_statistics not found")
        return

    with engine.begin() as conn:
        # Read all rows from general_statistics
        rows = conn.execute(text(f"SELECT * FROM {gs_table_name}")).fetchall()
        col_names = [c['name'] for c in inspector.get_columns(gs_table_name)]

        # Build sets of distinct values
        species_set = set()
        isolate_set = set()
        amr_gene_set = set()
        point_mutation_set = set()
        plasmid_cluster_set = set()
        virulence_factor_set = set()

        for row in rows:
            row_dict = dict(zip(col_names, row))

            val = row_dict.get('species')
            if val and val != '-':
                species_set.add(val)

            val = row_dict.get('isolate_id')
            if val and val != '-':
                isolate_set.add(val)

            # AMR genes (JSON array)
            val = row_dict.get('amr_gene_list')
            if val:
                genes = _extract_list(val)
                amr_gene_set.update(genes)

            # Point mutations (JSON array)
            val = row_dict.get('point_mutation_list')
            if val:
                mutations = _extract_list(val)
                point_mutation_set.update(mutations)

            # Plasmid clusters (JSON array)
            val = row_dict.get('plasmid_cluster_list')
            if val:
                clusters = _extract_list(val)
                plasmid_cluster_set.update(clusters)

            # Virulence factors (JSON array)
            val = row_dict.get('virulence_factor_list')
            if val:
                factors = _extract_list(val)
                virulence_factor_set.update(factors)

        # Insert new values into each meta table
        now = local_now()

        for table_name, values, value_col in [
            ('meta_species', species_set, 'species_name'),
            ('meta_isolates', isolate_set, 'isolate_id'),
            ('meta_amr_genes', amr_gene_set, 'amr_gene'),
            ('meta_point_mutations', point_mutation_set, 'point_mutation'),
            ('meta_plasmid_clusters', plasmid_cluster_set, 'plasmid_cluster'),
            ('meta_virulence_factors', virulence_factor_set, 'virulence_factor'),
        ]:
            table = metadata.tables[table_name]
            existing = set()
            for r in conn.execute(table.select()):
                existing.add(getattr(r, value_col))

            new_values = values - existing
            for v in sorted(new_values):
                conn.execute(table.insert().values(**{value_col: v, 'first_seen': now}))

            print(f"  {table_name}: {len(new_values)} new, {len(existing)} existing, {len(values)} total distinct")


def _extract_list(val):
    """Extract a list of strings from a JSON array or delimited string."""
    if isinstance(val, list):
        return [x.strip() for x in val if x and str(x).strip() and str(x).strip() != '-']
    if isinstance(val, str):
        val = val.strip()
        if val.startswith('['):
            try:
                parsed = json.loads(val)
                return [x.strip() for x in parsed if x and str(x).strip() and str(x).strip() != '-']
            except (json.JSONDecodeError, TypeError):
                pass
        # Fallback: split on semicolons or commas
        for delim in [';', ',']:
            if delim in val:
                return [x.strip() for x in val.split(delim) if x.strip() and x.strip() != '-']
        if val and val != '-':
            return [val]
    return []


def main():
    parser = argparse.ArgumentParser(description='Store EPIquence summary results in database')
    add_db_arguments(parser)
    parser.add_argument('--general-statistics', default=None, help='general_statistics.tsv')
    parser.add_argument('--taxonomy', default=None, help='taxonomy.tsv')
    parser.add_argument('--mlst', default=None, help='mlst.tsv')
    parser.add_argument('--amrfinder', default=None, help='amrfinder.tsv')
    parser.add_argument('--resfinder-phenotype', default=None, help='resfinder_phenotype.tsv')
    parser.add_argument('--point-mutations', default=None, help='point_mutations.tsv')
    parser.add_argument('--pointfinder', default=None, help='pointfinder.tsv')
    parser.add_argument('--resfinder-genes', default=None, help='resfinder_genes.tsv')
    parser.add_argument('--plasmid-resistance', default=None, help='plasmid_resistance_markers.tsv')
    parser.add_argument('--chromosomal-resistance', default=None, help='chromosomal_resistance_markers.tsv')
    parser.add_argument('--pairwise-distances', default=None, help='pairwise_distances.tsv')
    parser.add_argument('--cga-pairwise-distances', default=None, help='cga_pairwise_distances.tsv')
    parser.add_argument('--virulence-factors', default=None, help='virulence_factors.tsv')
    parser.add_argument('--ectyper', default=None, help='ectyper.tsv')
    parser.add_argument('--kleborate', default=None, help='kleborate.tsv')
    parser.add_argument('--emmtyper', default=None, help='emmtyper.tsv')
    parser.add_argument('--sistr', default=None, help='sistr.tsv')
    parser.add_argument('--seroba', default=None, help='seroba.tsv')
    parser.add_argument('--meningotype', default=None, help='meningotype.tsv')
    parser.add_argument('--elgato', default=None, help='elgato.tsv')
    parser.add_argument('--workflow-id', default='', help='Workflow ID')
    parser.add_argument('--workflow-description', default='', help='Workflow description')
    parser.add_argument('--group-id', default='', help='Group ID')
    parser.add_argument('--group-description', default='', help='Group description')
    parser.add_argument('--sample-count', type=int, default=0, help='Number of samples in this run')
    parser.add_argument('--assembly-types', default='', help='Assembly types used (comma-separated)')
    parser.add_argument('--nextflow-run-name', default='', help='Nextflow run name')
    parser.add_argument('--workflow-start', default='',
                        help="Nextflow's workflow.start, the true start of this "
                             "execution. Without it created_at would be the moment "
                             "this loader began, which is the very end of the run.")
    args = parser.parse_args()
    args.db_url = resolve_db(args)

    print(f"Connecting to database: {mask_url(args.db_url)}")
    # Loader start, used only for the load-duration message below.
    loader_start = local_now()

    # created_at is the workflow's start, not this script's.
    start_time = parse_nf_datetime(args.workflow_start) or loader_start

    try:
        engine = create_engine(args.db_url, echo=False)
        with engine.connect() as conn:
            conn.execute(text('SELECT 1'))
        print("  Database connection successful")
    except SQLAlchemyError as e:
        print(f"  ERROR: Failed to connect to database: {e}", file=sys.stderr)
        sys.exit(1)

    # Create workflow_log table if not exists
    wf_metadata = MetaData()
    workflow_log = Table('workflow_log', wf_metadata,
        Column('id', Integer, primary_key=True, autoincrement=True),
        Column('workflow_id', Text, nullable=False),
        Column('workflow_description', Text),
        Column('group_id', Text),
        Column('group_description', Text),
        Column('status', Text, nullable=False, default='Running'),
        Column('sample_count', Integer),
        Column('assembly_types', Text),
        Column('nextflow_run_name', Text),
        Column('tables_loaded', Text),
        Column('total_rows_inserted', Integer),
        Column('error_message', Text),
        Column('created_at', DateTime),
        Column('updated_at', DateTime),
        Column('completed_at', DateTime),
        Column('duration_seconds', Float),
        # How long the database load itself took, distinct from the workflow duration.
        Column('db_load_seconds', Float),
        # Attempt ordinal for this row: 1 for the original run, 2 for the first resume, and so on.
        Column('execution_count', Integer),
    )
    # create_all does not ALTER existing tables, so new columns are added by hand for databases
    # created before these were introduced.
    wf_metadata.create_all(engine)
    _add_missing_columns(engine, 'workflow_log', {
        'db_load_seconds': 'DOUBLE PRECISION',
        'execution_count': 'INTEGER',
    })

    # Count samples from general_statistics file
    sample_count = 0
    if args.general_statistics and os.path.exists(args.general_statistics):
        with open(args.general_statistics) as f:
            sample_count = sum(1 for _ in f) - 1  # subtract header

    # One row per execution, not per workflow.
    with engine.begin() as conn:
        # Reconcile any earlier attempt at this workflow that never closed out.
        stale = conn.execute(
            select(workflow_log.c.id, workflow_log.c.created_at, workflow_log.c.updated_at)
            .where(workflow_log.c.workflow_id == args.workflow_id)
            .where(workflow_log.c.status.in_(['Running', 'Loaded']))
        ).all()
        for stale_id, stale_created, stale_updated in stale:
            closed = stale_updated or stale_created
            vals = {'status': 'Failed',
                    'error_message': 'Run did not complete; superseded by a later execution',
                    'completed_at': closed}
            if closed and stale_created:
                vals['duration_seconds'] = round((closed - stale_created).total_seconds(), 2)
            conn.execute(workflow_log.update()
                         .where(workflow_log.c.id == stale_id).values(**vals))
        if stale:
            print(f"  Marked {len(stale)} incomplete earlier execution(s) as Failed")

        prior = conn.execute(
            select(func.count())
            .select_from(workflow_log)
            .where(workflow_log.c.workflow_id == args.workflow_id)
        ).scalar() or 0
        n_exec = prior + 1
        result = conn.execute(workflow_log.insert().values(
            workflow_id=args.workflow_id,
            workflow_description=args.workflow_description,
            group_id=args.group_id,
            group_description=args.group_description,
            status='Running',
            sample_count=sample_count,
            assembly_types=args.assembly_types,
            nextflow_run_name=args.nextflow_run_name,
            created_at=start_time,
            updated_at=start_time,
            execution_count=n_exec,
        ))
        log_id = result.inserted_primary_key[0]
    label = 'status=Running' if n_exec == 1 else f'execution #{n_exec}, resume'
    print(f"  Workflow log entry created (id={log_id}, {label})")

    metadata = MetaData()
    total = 0
    tables_loaded = []
    error_msg = None

    try:
        table_map = [
            ('ov_general_statistics', args.general_statistics),
            ('ov_taxonomy', args.taxonomy),
            ('ov_mlst', args.mlst),
            ('ov_amrfinder', args.amrfinder),
            ('ov_resfinder_phenotype', args.resfinder_phenotype),
            ('ov_point_mutations', args.point_mutations),
            ('ov_pointfinder', args.pointfinder),
            ('ov_resfinder_genes', args.resfinder_genes),
            ('ov_plasmid_resistance_markers', args.plasmid_resistance),
            ('ov_chromosomal_resistance_markers', args.chromosomal_resistance),
            ('ov_ska2_pairwise_distances', args.pairwise_distances),
            ('ov_cga_pairwise_distances', args.cga_pairwise_distances),
            ('ov_virulence_factors', args.virulence_factors),
            ('ov_ectyper', args.ectyper),
            ('ov_kleborate', args.kleborate),
            ('ov_emmtyper', args.emmtyper),
            ('ov_sistr', args.sistr),
            ('ov_seroba', args.seroba),
            ('ov_meningotype', args.meningotype),
            ('ov_elgato', args.elgato),
        ]

        for table_name, filepath in table_map:
            if filepath and os.path.exists(filepath):
                n = store_tsv(engine, metadata, table_name, filepath)
                total += n
                if n > 0:
                    tables_loaded.append(f"{table_name}:{n}")

        # Rebuild metadata lookup tables from general_statistics
        if args.general_statistics and os.path.exists(args.general_statistics):
            rebuild_meta_tables(engine, metadata)

        status = 'Success'
    except Exception as e:
        status = 'Failed'
        error_msg = str(e)
        print(f"  ERROR: {e}", file=sys.stderr)

    # Update workflow_log with load results.
    end_time = local_now()
    load_seconds = (end_time - loader_start).total_seconds()
    with engine.begin() as conn:
        conn.execute(
            workflow_log.update().where(workflow_log.c.id == log_id).values(
                status='Loaded' if status == 'Success' else status,
                total_rows_inserted=total,
                tables_loaded=', '.join(tables_loaded),
                error_message=error_msg,
                updated_at=end_time,
                db_load_seconds=round(load_seconds, 2),
            )
        )
    print(f"  Workflow log updated (status={status}, db load={load_seconds:.1f}s, rows={total})")

    print(f"\nTotal: {total} rows inserted across all tables")
    engine.dispose()


if __name__ == '__main__':
    main()