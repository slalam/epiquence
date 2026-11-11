"""Shared helper functions."""
import os
import re
import stat
import sys

DB_KEY_NAME = 'DB_URL'
API_KEY_NAME = 'API_KEY'

_LINE = re.compile(r'^\s*(?:export\s+)?([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$')


def read_settings_file(path):
    """Parse a ``KEY=value`` file."""
    values = {}
    try:
        st = os.stat(path)
    except OSError as exc:
        sys.exit(f'ERROR: settings file {path!r} is not readable ({exc})')
    if st.st_mode & (stat.S_IRWXG | stat.S_IRWXO):
        print(f'WARN: settings file {path} is accessible by group/other '
              f'(mode {stat.S_IMODE(st.st_mode):o}); run: chmod 600 {path}',
              file=sys.stderr)
    with open(path, encoding='utf-8') as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            m = _LINE.match(line)
            if not m:
                continue
            key, val = m.group(1), m.group(2)
            if len(val) >= 2 and val[0] == val[-1] and val[0] in ('"', "'"):
                val = val[1:-1]
            values[key] = val
    return values


def add_db_arguments(parser):
    """--settings-file (pipeline) or --db-url (manual runs)."""
    parser.add_argument('--settings-file', default=None, metavar='PATH',
                        help=f'KEY=value file holding {DB_KEY_NAME} (written by the pipeline)')
    parser.add_argument('--db-url', default=None,
                        help='SQLAlchemy URL, for running this script by hand')


def add_key_arguments(parser):
    """--settings-file (pipeline) or --api-key (manual runs)."""
    if not any(a.dest == 'settings_file' for a in parser._actions):
        parser.add_argument('--settings-file', default=None, metavar='PATH',
                            help=f'KEY=value file holding {API_KEY_NAME} (written by the pipeline)')
    parser.add_argument('--api-key', default='',
                        help='Metabase API key, for running this script by hand')


def _resolve(value, settings_file, key, label, required):
    if settings_file:
        got = read_settings_file(settings_file).get(key, '').strip()
        if got:
            return got
    if value:
        return value.strip()
    if required:
        sys.exit(f'ERROR: {label}: not found. Set surveillance_db_url in conf/env.config, '
                 f'or pass --db-url when running this script by hand.')
    return ''


def resolve_db(args):
    """Return the database URL (settings file first, then --db-url)."""
    return _resolve(args.db_url, args.settings_file, DB_KEY_NAME,
                    'database URL', required=True)


def resolve_key(args):
    """Return the Metabase API key, or '' when none is configured (optional)."""
    return _resolve(args.api_key, args.settings_file, API_KEY_NAME,
                    'Metabase API key', required=False)


def mask_url(url):
    """``postgresql+psycopg2://user:secret@host/db`` -> ``postgresql+psycopg2://user:***@host/db``."""
    return re.sub(r'(://[^/:@]+):[^@/]*@', r'\1:***@', url or '')
