#!/usr/bin/env python3
"""Refresh Metabase's view of the analysis database."""

import argparse
import sys
import urllib.error
import urllib.request

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.realpath(__file__)))
from epiquence_lib.utils import add_key_arguments, resolve_key


ENDPOINTS = [
    ('sync_schema', 'schema sync'),
    ('rescan_values', 'field value rescan'),
]


def post(url, api_key, timeout, label):
    """POST to one Metabase endpoint."""
    # data=b'' forces a POST; these endpoints expect no body.
    req = urllib.request.Request(
        url, data=b'', method='POST',
        headers={'x-api-key': api_key, 'Content-Type': 'application/json'})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            print(f'  {label}: triggered (HTTP {resp.status})')
            return True
    except urllib.error.HTTPError as e:
        # Reached Metabase, but it refused: wrong key, wrong db id, no permission.
        body = ''
        try:
            body = e.read().decode('utf-8', 'replace')[:500]
        except Exception:
            pass
        print(f'  WARN: {label} failed (HTTP {e.code} {e.reason})', file=sys.stderr)
        if body:
            print(f'    {body}', file=sys.stderr)
    except urllib.error.URLError as e:
        # Never reached Metabase: DNS, connection refused, TLS, timeout.
        print(f'  WARN: {label} could not reach {url} ({e.reason})', file=sys.stderr)
    except Exception as e:
        print(f'  WARN: {label} failed ({e})', file=sys.stderr)
    return False


def main():
    p = argparse.ArgumentParser(
        description="Trigger Metabase schema sync and field-value rescan")
    p.add_argument('--base-url', default='', help='e.g. https://metabase.example.org')
    p.add_argument('--db-id', default='', help='Metabase database id')
    add_key_arguments(p)
    p.add_argument('--timeout', type=int, default=120,
                   help='Seconds before giving up on each call')
    args = p.parse_args()

    base_url = args.base_url.strip().rstrip('/')
    db_id = args.db_id.strip()
    api_key = resolve_key(args)

    if not (base_url and db_id and api_key):
        print('Metabase refresh skipped: metabase_base_url, metabase_db_id and '
              'metabase_api_key '
              'must all be set.')
        return

    ok = 0
    for endpoint, label in ENDPOINTS:
        url = f'{base_url}/api/database/{db_id}/{endpoint}'
        print(f'POST {url}')
        if post(url, api_key, args.timeout, label):
            ok += 1

    print(f'Metabase refresh: {ok}/{len(ENDPOINTS)} call(s) succeeded')


if __name__ == '__main__':
    main()