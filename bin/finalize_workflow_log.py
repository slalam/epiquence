#!/usr/bin/env python3
"""Close out the workflow_log row for a run."""

import argparse
import re
import sys
from datetime import datetime

from sqlalchemy import create_engine, MetaData, Table, select, inspect

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib.utils import add_db_arguments, resolve_db


def parse_nf_datetime(value):
    """Parse Nextflow's workflow.start rendering."""
    if not value:
        return None
    v = str(value).strip()
    if not v or v.lower() in ('null', 'none'):
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
    # strptime cannot read a named zone abbreviation with %Z, so drop that token.
    m = re.match(r'^(\w{3} \w{3} +\d+ \d{2}:\d{2}:\d{2}) +[A-Z]{2,5} +(\d{4})$', v)
    if m:
        try:
            return datetime.strptime(f'{m.group(1)} {m.group(2)}', '%a %b %d %H:%M:%S %Y')
        except ValueError:
            pass
    print(f'WARN: could not parse workflow start {value!r}', file=sys.stderr)
    return None


def main():
    p = argparse.ArgumentParser(description='Finalize the workflow_log row for a run')
    add_db_arguments(p)
    p.add_argument('--workflow-id', required=True)
    p.add_argument('--workflow-start', default='',
                   help="Nextflow's workflow.start. Falls back to the row's own "
                        "created_at, which makes a manual re-run of this script "
                        "possible with nothing but --workflow-id.")
    p.add_argument('--status', default='Success')
    args = p.parse_args()
    args.db_url = resolve_db(args)

    if not args.workflow_id:
        print('finalize_workflow_log: no workflow id; nothing to update', file=sys.stderr)
        return

    # This process is last in the DAG, so "now" is the end of the run to within the time it took
    # to schedule this task.
    completed = datetime.now()
    started = parse_nf_datetime(args.workflow_start)

    try:
        engine = create_engine(args.db_url, echo=False)
        if not inspect(engine).has_table('workflow_log'):
            print('finalize_workflow_log: workflow_log table does not exist',
                  file=sys.stderr)
            return

        meta = MetaData()
        wf = Table('workflow_log', meta, autoload_with=engine)

        with engine.connect() as conn:
            row = conn.execute(
                select(wf.c.id, wf.c.created_at)
                .where(wf.c.workflow_id == args.workflow_id)
                .order_by(wf.c.id.desc())
            ).first()
        if not row:
            print(f'finalize_workflow_log: no workflow_log row for '
                  f'{args.workflow_id}; nothing to update', file=sys.stderr)
            return
        log_id, row_created = row

        if started is None and row_created:
            started = row_created
        duration = (completed - started).total_seconds() if started else None

        values = {
            'status': args.status,
            'completed_at': completed,
            'updated_at': completed,
        }
        if duration is not None:
            values['duration_seconds'] = round(duration, 2)
        if started:
            values['created_at'] = started

        with engine.begin() as conn:
            conn.execute(wf.update().where(wf.c.id == log_id).values(**values))

        dur = f'{duration:.1f}s' if duration is not None else 'unknown'
        print(f'  Workflow log finalized ({args.workflow_id}: {args.status}, '
              f'duration={dur})')
        engine.dispose()
    except Exception as exc:
        # Deliberately swallowed: everything else has already succeeded, and failing here would turn a
        # good run into a failed one.
        print(f'finalize_workflow_log: could not update log ({exc})', file=sys.stderr)


if __name__ == '__main__':
    main()