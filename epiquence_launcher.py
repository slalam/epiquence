import marimo

__generated_with = "0.20.1"
app = marimo.App(width="full", app_title="EPIquence Workflow Launcher")


@app.cell
def _():
    import marimo as mo
    import re
    import os
    import csv
    import json
    import glob
    import subprocess
    from datetime import datetime
    return mo, re, os, csv, json, glob, subprocess, datetime


@app.cell
def _(os, json, re, datetime):
    CONFIG_PATH = os.environ.get("EPIQUENCE_CONFIG", os.path.join(os.path.dirname(os.path.abspath(__file__)), "epiquence_config.json"))
    def load_config():
        if os.path.exists(CONFIG_PATH):
            with open(CONFIG_PATH, 'r') as _f:
                return json.load(_f)
        return {"pipeline_dir": ".", "output_base": "/data/epiquence/results", "reads_base": "/data/epiquence/reads", "groups": {"GP001": "Healthcare-Associated Infections"}}
    def save_config(cfg):
        with open(CONFIG_PATH, 'w') as _f:
            json.dump(cfg, _f, indent=2)
    config = load_config()

    def pipeline_root(cfg=None):
        """Pipeline folder from the config; a relative path is resolved against this launcher's
        folder (same rule as the launch and resume code)."""
        _cfg = cfg or load_config()
        _pl = _cfg.get("pipeline_dir", ".")
        if not os.path.isabs(_pl):
            _pl = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), _pl))
        return _pl

    # Settings the launcher reads from conf/env.config to show the effective run setup and to
    # catch missing databases before a run is started.
    ENV_KEYS = ('sylph_db', 'sylph_taxonomy', 'gtdbtk_db', 'bakta_db',
                'max_cpus', 'max_memory')

    def read_env_settings():
        _vals = {}
        _env = os.path.join(pipeline_root(), "conf", "env.config")
        try:
            with open(_env, 'r') as _f:
                for _line in _f:
                    _s = _line.strip()
                    if not _s or _s.startswith('//') or _s.startswith('#'):
                        continue
                    _m = re.match(r'^(\w+)\s*=\s*(.+?)\s*(?://.*)?$', _s)
                    if _m and _m.group(1) in ENV_KEYS:
                        _v = _m.group(2).strip()
                        if len(_v) >= 2 and _v[0] == _v[-1] and _v[0] in '"\'':
                            _v = _v[1:-1]
                        _vals[_m.group(1)] = _v
        except OSError:
            _vals['_missing'] = _env
        return _vals

    # Launcher runs are always one of the two capacity profiles, and the profile fixes the species
    # engine and the Bakta database.
    PROFILES = ('high_capacity', 'low_capacity')
    PROFILE_ENGINE = {'high_capacity': 'gtdbtk', 'low_capacity': 'sylph'}

    def profile_engine(profile):
        """Species engine a launcher run uses for this profile."""
        return PROFILE_ENGINE.get(profile, 'gtdbtk')

    def preflight(profile):
        """Problems that would stop the pipeline at startup (errors) or are worth knowing
        before launching (warnings), for the given profile."""
        _errors, _warnings = [], []
        _env = read_env_settings()
        if '_missing' in _env:
            _errors.append(f"conf/env.config not found at {_env['_missing']}")
            return _errors, _warnings
        if profile not in PROFILES:
            _errors.append(f"Unknown compute profile '{profile}'")
            return _errors, _warnings
        if profile_engine(profile) == 'sylph':
            _db = _env.get('sylph_db', '')
            if not _db:
                _errors.append("Low capacity uses sylph: set sylph_db in conf/env.config")
            elif not os.path.exists(_db):
                _errors.append(f"sylph_db not found: {_db}")
            _tax = [t.strip() for t in _env.get('sylph_taxonomy', '').split(',') if t.strip()]
            if not _tax:
                _errors.append("Low capacity uses sylph: set sylph_taxonomy (GTDB taxonomy TSVs) in conf/env.config")
            for _t in _tax:
                if not os.path.exists(_t):
                    _errors.append(f"sylph_taxonomy file not found: {_t}")
        else:
            _g = _env.get('gtdbtk_db', '')
            if not _g or not os.path.exists(_g):
                _errors.append(f"High capacity uses GTDB-Tk: gtdbtk_db not found: {_g or '(not set)'}")
            _b = _env.get('bakta_db', '')
            if not _b or not os.path.exists(_b):
                _warnings.append(f"bakta_db not found ({_b or 'not set'}): core genome analysis and "
                                 "plasmid comparison maps will fail with the high capacity profile.")
        if not (_env.get('max_cpus') and _env.get('max_memory')):
            _warnings.append("conf/env.config has no max_cpus / max_memory, so the profile uses its built-in "
                             "host size. Add them (see conf/env.config.template) so the whole machine is used.")
        if profile == 'high_capacity':
            try:
                _c = int(_env.get('max_cpus', '0'))
                _g = float(re.sub(r'[^0-9.]', '', _env.get('max_memory', '0')) or 0)
                if (_c and _c < 64) or (_g and 'GB' in _env.get('max_memory', '').upper() and _g < 256):
                    _warnings.append(f"This host ({_c} threads / {_env.get('max_memory')}) is below the high "
                                     "capacity target (64+ threads, 256+ GB); consider low capacity.")
            except ValueError:
                pass
        return _errors, _warnings

    def workflow_now():
        """Current time in the pipeline's configured zone."""
        try:
            _env = os.path.join(pipeline_root(), "conf", "env.config")
            with open(_env, 'r') as _f:
                _m = re.search(r'^\s*report_timezone\s*=\s*["\']([^"\']+)["\']',
                               _f.read(), re.M)
            if _m:
                from zoneinfo import ZoneInfo
                return datetime.now(ZoneInfo(_m.group(1)))
        except Exception:
            pass
        return datetime.now()

    def list_dirs_sorted(base, cfg=None):
        """Immediate subdirectories of `base` as {label: full_path} for a dropdown."""
        cfg = load_config() if cfg is None else cfg
        _mode = str(cfg.get("reads_sort", "name")).lower()
        if _mode not in ("name", "modified"):
            _mode = "name"
        _desc = bool(cfg.get("reads_sort_desc", _mode != "name"))
        _show_date = bool(cfg.get("reads_show_date", _mode != "name"))

        _entries = []
        if os.path.isdir(base):
            for _d in os.listdir(base):
                _full = os.path.join(base, _d)
                if _d.startswith('.') or not os.path.isdir(_full):
                    continue
                try:
                    _ts = os.stat(_full).st_mtime
                except OSError:
                    _ts = 0.0
                _entries.append((_d, _full, _ts))

        if _mode == "name":
            _entries.sort(key=lambda _e: _e[0].lower(), reverse=_desc)
        else:
            # name is the tiebreaker so same-second directories stay stable
            _entries.sort(key=lambda _e: (_e[2], _e[0].lower()), reverse=_desc)

        _opts = {}
        for _d, _full, _ts in _entries:
            _label = _d
            if _show_date and _ts:
                _label = f"{_d}  -  {datetime.fromtimestamp(_ts).strftime('%m-%d-%y')}"
            _opts[_label] = _full
        return _opts

    return CONFIG_PATH, load_config, save_config, config, workflow_now, list_dirs_sorted, pipeline_root, read_env_settings, profile_engine, preflight


@app.cell
def _(mo):
    mo.Html("""<style>
    @import url('https://fonts.googleapis.com/css2?family=DM+Sans:opsz,wght@9..40,300;9..40,400;9..40,500;9..40,600;9..40,700&family=Space+Mono:wght@400;700&display=swap');
    :root{--eq-surface:#fff;--eq-border:#e2ddd6;--eq-text:#2d2926;--eq-muted:#737780;--eq-accent:#D85A3B;--eq-accent-light:#fdf0ec;--eq-accent-dark:#b84430;--eq-green:#3d6b5e;--eq-green-light:#e8f0ec;--eq-green-dark:#2a4d42;--eq-warm-light:#fdf4ec;--eq-danger:#c44d4d;--eq-r:8px;--eq-control-h:34px;}
    body,.marimo{background:#f7f5f2 !important;}
    .eq-wrap{font-family:'DM Sans',sans-serif;color:var(--eq-text);max-width:1100px;margin:0 auto;padding:0 12px;}
    .eq-header{display:flex;align-items:center;gap:14px;padding:16px 0 12px;border-bottom:1px solid var(--eq-border);margin-bottom:16px;}
    .eq-header svg{height:28px;width:auto;}
    .eq-header-text{font-size:11px;color:var(--eq-muted);letter-spacing:1.5px;text-transform:uppercase;font-family:'Space Mono',monospace;}
    .eq-row{display:flex;gap:16px;margin-bottom:16px;}.eq-row>*{flex:1;min-width:0;}
    .eq-card{background:var(--eq-surface);border:1px solid var(--eq-border);border-radius:var(--eq-r);padding:18px 20px;}
    /* Field rows: the control shrinks to the card, long option labels truncate
       with an ellipsis instead of widening the row and pushing buttons out.
       The inline label is kept on one line so the control box does not grow
       taller than the button; both are then centred on the row. */
    .eq-fieldrow{display:flex;gap:8px;align-items:center;min-width:0;}
    .eq-field{flex:1 1 auto;min-width:0;max-width:100%;}
    .eq-field-fixed{flex:0 0 auto;display:flex;align-items:center;}
    .eq-fieldrow>*,.eq-fieldrow>*>*{margin-top:0;margin-bottom:0;}
    .eq-field *{max-width:100%;box-sizing:border-box;}
    .eq-field label{white-space:nowrap;align-items:center;}
    .eq-field select,.eq-field [role="combobox"],.eq-field input,.eq-field-fixed button{height:var(--eq-control-h);line-height:normal;box-sizing:border-box;margin:0;vertical-align:middle;}
    .eq-field select,.eq-field [role="combobox"],.eq-field input{width:100%;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}
    .eq-field label,.eq-field>*{width:100%;min-width:0;}
    .eq-title{font-family:'Space Mono',monospace;font-size:11px;font-weight:700;color:var(--eq-accent-dark);text-transform:uppercase;letter-spacing:1.2px;margin:0 0 14px 0;padding-bottom:8px;border-bottom:2px solid var(--eq-accent-light);}
    .eq-tbl{width:100%;border-collapse:separate;border-spacing:0;font-size:12px;border:1px solid var(--eq-border);border-radius:6px;overflow:hidden;}
    .eq-tbl th{background:var(--eq-green-light);color:var(--eq-green-dark);font-family:'Space Mono',monospace;font-size:10px;font-weight:700;text-transform:uppercase;letter-spacing:.8px;padding:7px 10px;text-align:left;border-bottom:2px solid var(--eq-border);}
    .eq-tbl td{padding:6px 10px;border-bottom:1px solid #f0ede8;vertical-align:middle;}
    .eq-tbl tr:last-child td{border-bottom:none;}.eq-tbl tr:hover td{background:var(--eq-warm-light);}
    .eq-badge{display:inline-block;padding:1px 6px;border-radius:3px;font-size:10px;font-weight:600;font-family:'Space Mono',monospace;}
    .eq-b-l{background:#e3f2fd;color:#1565c0;}.eq-b-r1{background:#e8f5e9;color:#2e7d32;}.eq-b-r2{background:#fff3e0;color:#e65100;}
    .eq-stats{display:flex;gap:12px;padding:10px 14px;background:var(--eq-green-light);border-radius:6px;margin-bottom:12px;}
    .eq-stat{font-family:'Space Mono',monospace;font-size:11px;color:var(--eq-green-dark);}.eq-stat strong{font-size:16px;display:block;margin-bottom:1px;}
    .eq-ok{background:linear-gradient(135deg,var(--eq-green-dark),var(--eq-green));color:#fff;border-radius:var(--eq-r);padding:28px;text-align:center;}
    .eq-ok h3{font-family:'Space Mono',monospace;font-size:18px;margin:12px 0 6px;}.eq-ok p{opacity:.85;font-size:13px;margin:3px 0;}
    .eq-ok-id{display:inline-block;background:rgba(255,255,255,.15);border:1px solid rgba(255,255,255,.25);border-radius:5px;padding:6px 16px;margin-top:12px;font-family:'Space Mono',monospace;font-size:16px;font-weight:700;letter-spacing:1px;}
    .eq-ok-det{margin-top:14px;font-size:11px;opacity:.7;font-family:'Space Mono',monospace;}
    .eq-err{background:#fef2f2;border:1px solid #fecaca;color:var(--eq-danger);border-radius:6px;padding:8px 12px;font-size:12px;margin-bottom:10px;}
    .eq-foot{text-align:center;padding:14px 0;font-size:10px;color:var(--eq-muted);letter-spacing:1px;text-transform:uppercase;font-family:'Space Mono',monospace;}
    .eq-upload-hint{font-size:12px;color:var(--eq-muted);margin-bottom:10px;}
    .eq-upload-hint code{background:var(--eq-accent-light);padding:1px 4px;border-radius:3px;font-size:11px;}
    .eq-asm-grid{display:flex;gap:14px;flex-wrap:wrap;margin-top:2px;}
    .eq-asm-item{display:flex;align-items:center;gap:4px;font-size:12px;}
    .eq-asm-code{font-family:'Space Mono',monospace;font-size:10px;background:var(--eq-green-light);color:var(--eq-green-dark);padding:1px 5px;border-radius:3px;}
    .eq-prog{width:100%;height:4px;background:#eee;border-radius:2px;overflow:hidden;margin:8px 0 4px;}
    .eq-prog-bar{height:100%;background:var(--eq-accent);border-radius:2px;}
    .eq-prog-text{font-size:10px;color:var(--eq-muted);font-family:'Space Mono',monospace;}
    .eq-paths{margin-top:10px;font-size:11px;color:var(--eq-muted);line-height:1.6;}.eq-paths strong{color:var(--eq-text);}
    .eq-grp-msg{font-size:12px;color:var(--eq-green);font-family:'Space Mono',monospace;padding:4px 0;}
    .eq-hr{border:none;border-top:1px solid var(--eq-border);margin:12px 0;}
    .eq-scan-info{font-size:11px;color:var(--eq-green-dark);margin-top:8px;font-family:'Space Mono',monospace;}
    .eq-base-path{font-size:11px;color:var(--eq-muted);margin-bottom:8px;font-family:'Space Mono',monospace;}
    </style>""")
    return


@app.cell
def _(mo):
    mo.Html("""<div class="eq-wrap"><div class="eq-header">
    <svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 195 45'><g stroke='none' fill='#D85A3B' transform='matrix(0.65,0,0,0.65,0,-8)'><path d='M72 22H28C12.56 22 0 34.564 0 50c0 13.236 10.764 24 24 24h48c13.236 0 24-10.764 24-24v-4c0-13.236-10.764-24-24-24zM57.916 46.752a379.337 379.337 0 01-6.132 4.1C59.16 55.24 67.176 58.84 76 59.716V42c0-1.1.9-2 2-2s2 .9 2 2v18c1.1 0 2 .9 2 2s-.9 2-2 2c-12.132 0-22.592-4.92-32-10.736C38.592 59.072 28.132 64 16 64c-1.104 0-2-.9-2-2s.896-2 2-2c4.192 0 8.188-.7 12.04-1.813-.008-.063-.04-.115-.04-.187v-8c0-1.1.896-2 2-2s2 .9 2 2v6.9c4.272-1.601 8.328-3.729 12.22-6.045a567.82 567.82 0 01-6.136-4.096C31.34 42.168 25.336 38.124 20 36.64V50c0 1.1-.896 2-2 2s-2-.9-2-2V36c-1.104 0-2-.9-2-2s.896-2 2-2c7.504 0 15.304 5.304 24.332 11.444 2.476 1.684 5.04 3.412 7.668 5.083 2.628-1.679 5.191-3.407 7.668-5.083C64.695 37.304 72.496 32 80 32c1.1 0 2 .9 2 2s-.9 2-2 2c-4.721 0-10.04 2.824-16 6.68V50c0 1.1-.9 2-2 2s-2-.9-2-2v-4.656c-.687.468-1.376.928-2.084 1.408z'/></g><g fill='#D85A3B' transform='matrix(1.8,0,0,1.8,68,4)'><path d='M6.57 13.87C5.52 14.10 3.92 14.32 2.60 14.35L2.59 10.30L5.94 10.30C5.99 10.08 6.01 9.81 5.94 9.56L2.59 9.56L2.58 5.99C3.89 6.02 5.12 6.17 6.29 6.43C6.65 6.27 6.80 5.31 6.43 5.31L2.23 5.31C1.68 5.31 1.41 5.38 1.12 5.64C1.27 6.09 1.32 6.59 1.32 7.15L1.32 14.32C1.32 14.77 1.62 15.02 2.06 15.02L6.50 15.02C6.93 15.02 6.90 14.08 6.57 13.87ZM10.77 5.12C9.63 5.12 9.03 5.28 8.68 5.45C8.83 5.91 8.88 6.71 8.88 7.01L8.88 14.49C8.88 14.92 9.18 15.19 9.62 15.19C9.98 15.19 10.25 15.01 10.36 14.90C10.19 14.48 10.16 13.72 10.16 13.40L10.16 11.17L10.53 11.17C13.36 11.17 14.83 10.09 14.83 7.99C14.83 6.01 13.48 5.12 10.77 5.12ZM10.54 10.64L10.16 10.64L10.16 5.70C10.29 5.67 10.49 5.67 10.74 5.67C12.67 5.67 13.43 6.52 13.43 8.08C13.43 9.72 12.59 10.64 10.54 10.64ZM17.51 15.19C17.88 15.19 18.13 15.01 18.26 14.90C18.09 14.45 18.06 13.89 18.06 13.34L18.06 5.85C18.06 5.42 17.74 5.17 17.30 5.17C16.95 5.17 16.69 5.33 16.58 5.45C16.73 5.91 16.77 6.47 16.77 7.01L16.77 14.49C16.77 14.92 17.08 15.19 17.51 15.19ZM25.17 18.19C25.65 18.19 25.93 17.92 25.93 17.50L25.93 8.64C25.37 8.11 24.29 7.70 23.06 7.70C20.68 7.70 19.89 9.73 19.89 11.66C19.89 14.15 20.83 15.22 22.33 15.22C23.46 15.22 24.32 14.50 24.65 13.50L24.65 15.20C24.65 16.14 24.67 17.19 24.44 17.89C24.60 18.05 24.82 18.20 25.17 18.19ZM22.65 14.48C21.70 14.48 21.25 13.57 21.25 11.70C21.24 10.14 21.50 8.30 23.11 8.30C23.77 8.30 24.39 8.55 24.67 8.79L24.67 11.24C24.67 12.52 24.09 14.48 22.65 14.48ZM34.31 14C34.22 14.17 34.05 14.20 33.94 14.20C33.67 14.20 33.56 14 33.56 13.48L33.54 8.46C33.54 7.98 33.22 7.70 32.82 7.70C32.45 7.70 32.20 7.88 32.09 7.98C32.30 8.57 32.30 9.06 32.30 11.23C32.30 13.61 31.35 14.55 30.46 14.55C29.47 14.55 29.15 13.85 29.15 12.12L29.15 8.40C29.15 7.95 28.84 7.70 28.42 7.70C28.08 7.70 27.86 7.85 27.69 7.99C27.85 8.44 27.89 9.07 27.89 9.55L27.89 12.29C27.89 14.43 28.59 15.22 30.09 15.22C31.21 15.22 32.05 14.59 32.42 13.48L32.42 13.66C32.42 14.67 32.84 15.16 33.59 15.16C34.15 15.16 34.47 14.88 34.47 14.46C34.47 14.31 34.43 14.15 34.31 14ZM41.09 13.17C40.80 13.93 40.01 14.36 39.00 14.36C37.79 14.36 36.85 13.72 36.72 11.79C37.21 11.84 37.70 11.87 38.21 11.87C40.14 11.87 41.36 11.14 41.36 9.74C41.36 8.40 40.25 7.70 38.74 7.70C36.54 7.70 35.35 9.16 35.35 11.48C35.35 13.90 36.64 15.20 38.68 15.20C39.89 15.20 40.85 14.76 41.20 13.89C41.40 13.43 41.27 13.24 41.09 13.17ZM38.68 8.23C39.33 8.23 40.00 8.51 40.00 9.56C40.00 11.00 38.88 11.33 37.65 11.33C37.32 11.31 37.00 11.30 36.71 11.27C36.72 9.28 37.49 8.23 38.68 8.23ZM48.99 13.34L48.99 10.60C48.99 8.69 48.27 7.70 46.73 7.70C45.58 7.70 44.74 8.32 44.35 9.39L44.35 9.23C44.35 8.22 43.93 7.71 43.19 7.71C42.63 7.71 42.29 8.01 42.29 8.43C42.29 8.60 42.35 8.74 42.46 8.89C42.55 8.72 42.73 8.69 42.83 8.69C43.09 8.69 43.20 8.89 43.22 9.41L43.22 14.43C43.22 14.91 43.53 15.19 43.96 15.19C44.31 15.19 44.58 15.02 44.69 14.91C44.48 14.34 44.48 13.83 44.48 11.66C44.48 9.32 45.46 8.41 46.35 8.41C47.39 8.41 47.71 9.24 47.73 10.77L47.73 14.49C47.73 14.92 48.05 15.19 48.47 15.19C48.79 15.19 49.03 15.05 49.20 14.90C49.03 14.45 49 13.82 48.99 13.34ZM53.93 15.20C55.19 15.20 56.15 14.67 56.46 13.68C56.60 13.20 56.49 13.01 56.31 12.92C56.08 13.86 55.27 14.36 54.28 14.36C52.98 14.36 51.98 13.54 51.98 11.28C51.98 9.28 52.77 8.23 53.90 8.23C54.60 8.23 55.17 8.55 55.17 9.38C55.17 9.83 55.03 10.16 54.73 10.30C54.89 10.51 55.17 10.60 55.44 10.60C56.03 10.60 56.45 10.15 56.45 9.48C56.45 8.47 55.61 7.70 54.03 7.70C51.84 7.70 50.62 9.14 50.62 11.48C50.62 13.90 51.86 15.20 53.93 15.20ZM63.46 13.17C63.17 13.93 62.38 14.36 61.38 14.36C60.16 14.36 59.22 13.72 59.09 11.79C59.58 11.84 60.07 11.87 60.58 11.87C62.51 11.87 63.73 11.14 63.73 9.74C63.73 8.40 62.62 7.70 61.11 7.70C58.91 7.70 57.72 9.16 57.72 11.48C57.72 13.90 59.01 15.20 61.05 15.20C62.26 15.20 63.22 14.76 63.57 13.89C63.77 13.43 63.64 13.24 63.46 13.17ZM61.05 8.23C61.70 8.23 62.37 8.51 62.37 9.56C62.37 11.00 61.25 11.33 60.02 11.33C59.70 11.31 59.37 11.30 59.08 11.27C59.09 9.28 59.86 8.23 61.05 8.23Z'/></g></svg>
    <span class="eq-header-text">Workflow Launcher</span>
    </div></div>""")
    return


@app.cell
def _(mo, os, config, list_dirs_sorted):
    workflow_desc = mo.ui.text(label="Workflow Description *", full_width=True)
    _groups = config.get("groups", {})
    _gopts = {f"{_k} - {_v}": _k for _k, _v in sorted(_groups.items())}
    group_dropdown = mo.ui.dropdown(options=_gopts, label="Select Group", value=list(_gopts.keys())[0] if _gopts else None)
    new_group_input = mo.ui.text(placeholder="Create a new group", full_width=True)
    add_group_btn = mo.ui.run_button(label="Add", kind="neutral")
    asm_s = mo.ui.checkbox(label="Illumina", value=True)
    asm_l = mo.ui.checkbox(label="Oxford Nanopore", value=False)
    asm_hl = mo.ui.checkbox(label="Hybrid Long Read First", value=False)
    asm_hs = mo.ui.checkbox(label="Hybrid Short Read First", value=False)

    # Reads directory dropdown from reads_base, ordered per "reads_sort" in the config
    _rb = config.get("reads_base", "/data/epiquence/reads")
    _dir_opts = list_dirs_sorted(_rb, config)
    reads_dropdown = mo.ui.dropdown(options=_dir_opts, label="Reads Directory")
    assemblies_dropdown = mo.ui.dropdown(options=_dir_opts, label="Assemblies Directory")
    scan_btn = mo.ui.run_button(label="Scan", kind="neutral")

    outbreak_check = mo.ui.checkbox(label="Enable Outbreak Detection", value=True)
    core_genome_check = mo.ui.checkbox(label="Core Genome Analysis", value=False)
    mash_threshold = mo.ui.number(value=0.05, start=0.001, stop=0.1, step=0.001, label="Outbreak Mash Threshold")

    skip_rasusa = mo.ui.checkbox(label="Skip Downsampling (use full reads)", value=True)
    rasusa_coverage = mo.ui.number(value=100, start=10, stop=500, step=10, label="Target Coverage")
    rasusa_genome_size = mo.ui.number(value=6, start=1, stop=15, step=0.5, label="Genome Size (Mb)")

    hgt_check = mo.ui.checkbox(label="Enable Plasmid HGT Detection", value=True)
    hgt_threshold = mo.ui.number(value=0.005, start=0.001, stop=0.05, step=0.001, label="HGT Mash Threshold")
    hgt_sketch = mo.ui.number(value=5000, start=1000, stop=50000, step=1000, label="Sketch Size")
    hgt_maps_check = mo.ui.checkbox(label="Generate ORF Comparisons", value=True)

    surveillance_db_check = mo.ui.checkbox(label="Enable Surveillance", value=True)

    _profile_opts = {
        "High capacity | 64+ threads, 256+ GB RAM": "high_capacity",
        "Low capacity | 8 to 32 threads, 16 to 64 GB RAM": "low_capacity",
    }
    _def_profile = config.get("default_profile") or "high_capacity"
    capacity_profile = mo.ui.dropdown(
        options=_profile_opts, label="Compute Profile", full_width=True,
        value=next((k for k, v in _profile_opts.items() if v == _def_profile), list(_profile_opts)[0]))

    return capacity_profile, workflow_desc, group_dropdown, new_group_input, add_group_btn, asm_s, asm_l, asm_hl, asm_hs, reads_dropdown, assemblies_dropdown, scan_btn, outbreak_check, core_genome_check, mash_threshold, hgt_check, hgt_threshold, hgt_sketch, hgt_maps_check, rasusa_coverage, rasusa_genome_size, skip_rasusa, surveillance_db_check


@app.cell
def _(mo, os, reads_dropdown):
    # Metadata choices: CSV/TSV files in the selected Reads Directory (rebuilt when the folder changes).
    _opts = {"None": ""}
    _dir = reads_dropdown.value
    if _dir and os.path.isdir(_dir):
        for _f in sorted(os.listdir(_dir), key=str.lower):
            _p = os.path.join(_dir, _f)
            if not _f.startswith('.') and os.path.isfile(_p) and _f.lower().endswith(('.csv', '.tsv')):
                _opts[_f] = _p
    metadata_path = mo.ui.dropdown(options=_opts, value="None", label="Metadata File")
    return (metadata_path,)


@app.cell
def _(mo, config, capacity_profile, read_env_settings, preflight, workflow_desc, group_dropdown, new_group_input, add_group_btn, asm_s, asm_l, asm_hl, asm_hs, reads_dropdown, assemblies_dropdown, scan_btn, outbreak_check, core_genome_check, metadata_path, mash_threshold, hgt_check, hgt_threshold, hgt_sketch, hgt_maps_check, rasusa_coverage, rasusa_genome_size, skip_rasusa, surveillance_db_check):
    _rb = config.get("reads_base", "/data/epiquence/reads")

    _group_section = f"""<hr class="eq-hr">
    <div class="eq-fieldrow">
        <div class="eq-field">{mo.as_html(group_dropdown).text}</div>
    </div>
    <div class="eq-fieldrow" style="margin-top:6px;">
        <div class="eq-field">{mo.as_html(new_group_input).text}</div>
        <div class="eq-field-fixed">{mo.as_html(add_group_btn).text}</div>
    </div>
    <hr class="eq-hr">"""

    _left = mo.Html(f"""<div style="display:flex;flex-direction:column;flex:1;gap:16px;">
    <div class="eq-card"><div class="eq-title">\u2699 Workflow Configuration</div>
    {mo.as_html(workflow_desc).text}
    {_group_section}
    <div style="font-size:11px;font-weight:600;color:var(--eq-muted);text-transform:uppercase;letter-spacing:1px;margin:8px 0 2px;">Surveillance</div>
    <div style="display:flex;gap:12px;align-items:center;">
        {mo.as_html(surveillance_db_check).text}
    </div>
    <hr class="eq-hr">
    <div style="font-size:11px;font-weight:600;color:var(--eq-muted);text-transform:uppercase;letter-spacing:1px;margin:8px 0 2px;">Assembly Types</div>
    <div class="eq-asm-grid">
        <div class="eq-asm-item">{mo.as_html(asm_s).text} <span class="eq-asm-code">S</span></div>
        <div class="eq-asm-item">{mo.as_html(asm_l).text} <span class="eq-asm-code">L</span></div>
        <div class="eq-asm-item">{mo.as_html(asm_hl).text} <span class="eq-asm-code">HL</span></div>
        <div class="eq-asm-item">{mo.as_html(asm_hs).text} <span class="eq-asm-code">HS</span></div>
    </div>
    <hr class="eq-hr">
    <div style="font-size:11px;font-weight:600;color:var(--eq-muted);text-transform:uppercase;letter-spacing:1px;margin:0 0 2px;">Long Read Subsampling</div>
    <div style="margin-bottom:4px;">{mo.as_html(skip_rasusa).text}</div>
    <div style="display:flex;gap:12px;">
        <div style="flex:1;">{mo.as_html(rasusa_coverage).text}</div>
        <div style="flex:1;">{mo.as_html(rasusa_genome_size).text}</div>
    </div>
    <div class="eq-upload-hint">Subsample long reads to target coverage using Rasusa. Check above to skip downsampling and use full reads.</div>
    </div>
    <div class="eq-card" style="flex:1;"><div class="eq-title">\U0001F4C1 Sequencing Reads</div>
    <div class="eq-upload-hint">
        Select a directory containing <code>.fastq.gz</code> files:<br>
        <code>&lt;id&gt;_R1.fastq.gz</code> <code>&lt;id&gt;_R2.fastq.gz</code> (short) \u00b7 <code>&lt;id&gt;_L.fastq.gz</code> (long)
    </div>
    <div class="eq-fieldrow">
        <div class="eq-field">{mo.as_html(reads_dropdown).text}</div>
        <div class="eq-field-fixed">{mo.as_html(scan_btn).text}</div>
    </div>
    <div class="eq-upload-hint" style="margin-top:12px;">
        Optional: pre-assembled <code>.fasta</code> / <code>.fna</code> files (samples get <code>_A</code> suffix, skip QC/assembly)
    </div>
    <div class="eq-fieldrow">
        <div class="eq-field">{mo.as_html(assemblies_dropdown).text}</div>
    </div>
    </div>
    </div>""")

    # ── Compute profile: live summary of what Nextflow will do ──
    _prof = capacity_profile.value or "high_capacity"
    _env = read_env_settings()
    if _env.get('max_cpus') and _env.get('max_memory'):
        _host = f"{_env['max_cpus']} threads | {_env['max_memory']} (conf/env.config)"
    else:
        _host = "not set in conf/env.config (max_cpus / max_memory)"
    _errs, _warns = preflight(_prof)
    _note_style = {"err": ("rgba(216,90,59,.12)", "\u26D4"), "warn": ("rgba(230,170,40,.14)", "\u26A0")}
    _notes = "".join(
        '<div style="margin-top:6px;padding:6px 10px;border-radius:6px;font-size:12px;'
        f'background:{_note_style[_kind][0]};">{_note_style[_kind][1]} {_msg}</div>'
        for _kind, _list in (("err", _errs), ("warn", _warns)) for _msg in _list)
    _compute_card = f"""<div class="eq-card"><div class="eq-title">\u26A1 Compute Profile</div>
    <div>{mo.as_html(capacity_profile).text}</div>
    <div class="eq-upload-hint" style="margin-top:8px;">
        Host size: {_host}
    </div>
    {_notes}
    </div>"""

    _right = mo.Html(f"""<div style="display:flex;flex-direction:column;flex:1;gap:16px;">
    {_compute_card}
    <div class="eq-card" style="flex:1;"><div class="eq-title">\U0001F50D Outbreak Detection</div>
    <div style="margin-bottom:6px;">{mo.as_html(outbreak_check).text}</div>
    <div style="margin-bottom:4px;">{mo.as_html(mash_threshold).text}</div>
    <div class="eq-upload-hint" style="margin-bottom:10px;">Mash distance cutoff for grouping same-species samples into clusters. Lower values require closer relatedness.</div>
    <div style="margin-bottom:6px;">{mo.as_html(core_genome_check).text}</div>
    <div class="eq-upload-hint" style="margin-bottom:10px;">Run Bakta + Panaroo core genome alignment for clusters with 5+ samples. Produces additional SNP matrix and ML tree based on shared core genes.</div>
    <div>{mo.as_html(metadata_path).text}</div>
    <div class="eq-upload-hint" style="margin-top:8px;">
        Optional CSV/TSV from the selected Reads Directory, with a <code>sample_id</code> column, to annotate phylogenetic trees with epidemiological context.
    </div>
    </div>
    <div class="eq-card" style="flex:1;"><div class="eq-title">\U0001F517 Plasmid HGT Detection</div>
    <div style="margin-bottom:6px;">{mo.as_html(hgt_check).text}</div>
    <div style="margin-bottom:6px;">{mo.as_html(hgt_threshold).text}</div>
    <div style="margin-bottom:6px;">{mo.as_html(hgt_sketch).text}</div>
    <div style="margin-bottom:6px;">{mo.as_html(hgt_maps_check).text}</div>
    <div class="eq-upload-hint" style="margin-top:8px;">
        Detects potential horizontal gene transfer by clustering plasmids using Mash distance and generates plasmid comparison visualizations with gene annotations.
    </div>
    </div>
    </div>""")

    mo.Html(f'<div class="eq-wrap"><div class="eq-row" style="align-items:stretch;">{_left.text}{_right.text}</div></div>')
    return


@app.cell
def _(add_group_btn, new_group_input, save_config, load_config):
    add_group_msg = ""
    if add_group_btn and add_group_btn.value and new_group_input.value.strip():
        _cfg = load_config()
        _grps = _cfg.get("groups", {})
        _nums = []
        for _k in _grps:
            try:
                _nums.append(int(_k.replace("GP", "")))
            except ValueError:
                pass
        _nxt = max(_nums) + 1 if _nums else 1
        _nid = f"GP{_nxt:03d}"
        _grps[_nid] = new_group_input.value.strip()
        _cfg["groups"] = _grps
        save_config(_cfg)
        add_group_msg = f"\u2713 Created {_nid}: {new_group_input.value.strip()} | reload the page to see it in the dropdown"
    return (add_group_msg,)


@app.cell
def _(mo, add_group_msg):
    mo.Html(f'<div class="eq-wrap"><div class="eq-grp-msg">{add_group_msg}</div></div>') if add_group_msg else mo.Html("")
    return


@app.cell
def _(mo, re, os, glob, scan_btn, reads_dropdown, assemblies_dropdown):
    _sd = {}
    _pe = []
    scan_status = ""
    _total = 0
    _asm_samples = []

    has_reads_dir = reads_dropdown and reads_dropdown.value
    has_asm_dir = assemblies_dropdown and assemblies_dropdown.value

    if scan_btn and scan_btn.value and (has_reads_dir or has_asm_dir):
        if has_reads_dir:
            _dir = reads_dropdown.value
            if os.path.isdir(_dir):
                _files = sorted(glob.glob(os.path.join(_dir, "*.fastq.gz")) + glob.glob(os.path.join(_dir, "*.fq.gz")))
                _total = len(_files)
                for _fp in _files:
                    _fn = os.path.basename(_fp)
                    _b = _fn.replace('.fastq.gz', '').replace('.fq.gz', '')
                    _mr1 = re.match(r'^(.+)_R1$', _b)
                    _mr2 = re.match(r'^(.+)_R2$', _b)
                    _ml = re.match(r'^(.+)_L$', _b)
                    if _mr1:
                        _k = _mr1.group(1); _sd.setdefault(_k, {'name': _k, 'R1': None, 'R2': None, 'L': None, 'R1_path': None, 'R2_path': None, 'L_path': None})
                        _sd[_k]['R1'] = _fn; _sd[_k]['R1_path'] = _fp
                    elif _mr2:
                        _k = _mr2.group(1); _sd.setdefault(_k, {'name': _k, 'R1': None, 'R2': None, 'L': None, 'R1_path': None, 'R2_path': None, 'L_path': None})
                        _sd[_k]['R2'] = _fn; _sd[_k]['R2_path'] = _fp
                    elif _ml:
                        _k = _ml.group(1); _sd.setdefault(_k, {'name': _k, 'R1': None, 'R2': None, 'L': None, 'R1_path': None, 'R2_path': None, 'L_path': None})
                        _sd[_k]['L'] = _fn; _sd[_k]['L_path'] = _fp
                    else:
                        _pe.append(_fn)

        if has_asm_dir:
            _adir = assemblies_dropdown.value
            if os.path.isdir(_adir):
                _asm_files = sorted(
                    glob.glob(os.path.join(_adir, "*.fasta")) +
                    glob.glob(os.path.join(_adir, "*.fna")) +
                    glob.glob(os.path.join(_adir, "*.fa"))
                )
                for _af in _asm_files:
                    _abase = os.path.splitext(os.path.basename(_af))[0]
                    _asm_samples.append({'name': _abase, 'path': _af})

        parts = []
        if has_reads_dir and _total > 0:
            parts.append(f"{_total} read files ({len(_sd)} samples)")
        elif has_reads_dir:
            parts.append("no read files found")
        if has_asm_dir and _asm_samples:
            parts.append(f"{len(_asm_samples)} assemblies")
        elif has_asm_dir:
            parts.append("no assembly files found")
        scan_status = "\u2713 Scanned: " + ", ".join(parts) if parts else "No files found"

    elif scan_btn and scan_btn.value and not has_reads_dir and not has_asm_dir:
        scan_status = "Select a reads or assemblies directory first"

    sample_list = sorted(_sd.values(), key=lambda _x: _x['name'])
    assembly_sample_list = _asm_samples
    parse_errors = list(_pe)
    name_inputs = {_it['name']: mo.ui.text(value=_it['name'], full_width=True) for _it in sample_list}
    select_toggle = mo.ui.switch(label="Select All", value=True) if sample_list or assembly_sample_list else None
    return sample_list, assembly_sample_list, parse_errors, name_inputs, select_toggle, scan_status


@app.cell
def _(mo, sample_list, assembly_sample_list, select_toggle):
    _all_on = select_toggle.value if select_toggle else True
    sample_checks = {_it['name']: mo.ui.checkbox(value=_all_on) for _it in sample_list}
    asm_checks = {_it['name']: mo.ui.checkbox(value=_all_on) for _it in assembly_sample_list}
    return (sample_checks, asm_checks)


@app.cell
def _(mo, scan_status):
    def _show_scan():
        if not scan_status:
            return mo.Html("")
        _ok = "\u2713" in scan_status
        _c = "var(--eq-green-dark)" if _ok else "var(--eq-danger)"
        _w = "100%" if _ok else "0%"
        return mo.Html(f'<div class="eq-wrap"><div class="eq-prog"><div class="eq-prog-bar" style="width:{_w};background:{_c};"></div></div><div class="eq-prog-text" style="color:{_c};">{scan_status}</div></div>')
    _show_scan()
    return


@app.cell
def _(mo, os, sample_list, assembly_sample_list, parse_errors, name_inputs, sample_checks, asm_checks, select_toggle):
    def _tbl():
        if not sample_list and not assembly_sample_list:
            return mo.Html('<div class="eq-wrap"><div class="eq-card" style="text-align:center;padding:28px;color:var(--eq-muted);"><div style="font-size:28px;margin-bottom:8px;">\U0001F4C2</div><div style="font-size:13px;">No samples detected</div><div style="font-size:11px;margin-top:4px;">Select a reads and/or assemblies directory, then click Scan</div></div></div>')
        _nt = len(sample_list)
        _na = len(assembly_sample_list)
        _nl = sum(1 for _z in sample_list if _z['L'])
        _ns = sum(1 for _z in sample_list if _z['R1'] and _z['R2'])
        _nh = sum(1 for _z in sample_list if _z['L'] and _z['R1'] and _z['R2'])
        _al = ""
        for _e in parse_errors:
            _al += f'<div class="eq-err">\u26A0 Unrecognized: <code>{_e}</code></div>'
        for _z in sample_list:
            if _z['R1'] and not _z['R2']:
                _al += f'<div class="eq-err">\u26A0 {_z["name"]}: R1 found, R2 missing</div>'
            if _z['R2'] and not _z['R1']:
                _al += f'<div class="eq-err">\u26A0 {_z["name"]}: R2 found, R1 missing</div>'
        def _bg(_t, _v):
            if not _v:
                return '<span style="color:#ccc;">none</span>'
            _c = {'L': 'eq-b-l', 'R1': 'eq-b-r1', 'R2': 'eq-b-r2'}[_t]
            return f'<span class="eq-badge {_c}">{_t}</span> <span style="font-size:10px;color:var(--eq-muted);">{_v}</span>'
        _toggle_html = mo.as_html(select_toggle).text if select_toggle else ""
        _sel_bar = f'<div style="display:flex;gap:8px;align-items:center;margin-bottom:8px;"><div style="margin-left:auto;">{_toggle_html}</div></div>'
        _tr = ""
        for _i, _z in enumerate(sample_list):
            _w = mo.as_html(name_inputs[_z['name']]).text
            _chk = mo.as_html(sample_checks[_z['name']]).text
            _tr += f'<tr><td style="width:30px;text-align:center;">{_chk}</td><td style="width:30px;text-align:center;color:var(--eq-muted);font-family:Space Mono,monospace;font-size:10px;">{_i+1}</td><td style="width:180px;">{_w}</td><td>{_bg("L", _z["L"])}</td><td>{_bg("R1", _z["R1"])}</td><td>{_bg("R2", _z["R2"])}</td></tr>'

        _asm_html = ""
        if assembly_sample_list:
            _asm_rows = ""
            for _i, _z in enumerate(assembly_sample_list):
                _chk = mo.as_html(asm_checks[_z['name']]).text
                _fname = os.path.basename(_z['path'])
                _asm_rows += f'<tr><td style="width:30px;text-align:center;">{_chk}</td><td style="width:30px;text-align:center;color:var(--eq-muted);font-family:Space Mono,monospace;font-size:10px;">{_i+1}</td><td>{_z["name"]}_A</td><td><span class="eq-badge" style="background:#e8daef;color:#6c3483;">ASM</span> <span style="font-size:10px;color:var(--eq-muted);">{_fname}</span></td></tr>'
            _asm_html = f'<div style="margin-top:12px;font-size:12px;font-weight:600;color:#6c3483;">Assembly Inputs ({_na})</div><table class="eq-tbl" style="margin-top:4px;"><thead><tr><th></th><th>#</th><th>Sample ID</th><th>Assembly</th></tr></thead><tbody>{_asm_rows}</tbody></table>'

        _stats_html = f'<div class="eq-stats"><div class="eq-stat"><strong>{_nt}</strong>Read Samples</div><div class="eq-stat"><strong>{_nl}</strong>Long</div><div class="eq-stat"><strong>{_ns}</strong>Short Pairs</div><div class="eq-stat"><strong>{_nh}</strong>Hybrid Ready</div>'
        if _na > 0:
            _stats_html += f'<div class="eq-stat"><strong>{_na}</strong>Assemblies</div>'
        _stats_html += '</div>'

        _reads_table = ""
        if sample_list:
            _reads_table = f'<table class="eq-tbl"><thead><tr><th></th><th>#</th><th>Sample ID (editable)</th><th>Long</th><th>R1</th><th>R2</th></tr></thead><tbody>{_tr}</tbody></table>'

        return mo.Html(f'<div class="eq-wrap"><div class="eq-card"><div class="eq-title">\U0001F9EC Detected Samples</div>{_al}{_stats_html}{_sel_bar}{_reads_table}{_asm_html}</div></div>')
    _tbl()
    return


@app.cell
def _(mo, sample_list, assembly_sample_list):
    submit_btn = mo.ui.run_button(label="Launch Workflow", kind="danger") if (sample_list or assembly_sample_list) else None
    mo.Html(f'<div class="eq-wrap" style="text-align:center;margin-bottom:12px;">{mo.as_html(submit_btn).text}</div>') if submit_btn else mo.Html("")
    return (submit_btn,)


@app.cell
def _(mo, os, csv, subprocess, datetime, submit_btn, sample_list, assembly_sample_list, asm_checks, name_inputs, sample_checks, capacity_profile, profile_engine, preflight, workflow_desc, group_dropdown, asm_s, asm_l, asm_hl, asm_hs, load_config, outbreak_check, core_genome_check, metadata_path, mash_threshold, hgt_check, hgt_threshold, hgt_sketch, rasusa_coverage, rasusa_genome_size, skip_rasusa, assemblies_dropdown, surveillance_db_check, hgt_maps_check, workflow_now):
    def _go():
        if not submit_btn or not submit_btn.value:
            return mo.Html("")
        if not sample_list and not assembly_sample_list:
            return mo.Html('<div class="eq-wrap"><div class="eq-err">No samples detected.</div></div>')

        # Filter to only checked samples
        _selected = [_it for _it in sample_list if sample_checks.get(_it['name'], mo.ui.checkbox(value=True)).value]
        _selected_asm = [_it for _it in assembly_sample_list if asm_checks.get(_it['name'], mo.ui.checkbox(value=True)).value]

        if not _selected and not _selected_asm:
            return mo.Html('<div class="eq-wrap"><div class="eq-err">No samples selected. Use the checkboxes to select samples to include.</div></div>')
        if not workflow_desc.value or not workflow_desc.value.strip():
            return mo.Html('<div class="eq-wrap"><div class="eq-err">Workflow Description is required.</div></div>')

        # Same checks the pipeline makes at startup (missing sylph/GTDB-Tk databases), surfaced here
        # instead of in the log.
        _prof = capacity_profile.value or "high_capacity"
        _eng = profile_engine(_prof)
        _pf_errors, _ = preflight(_prof)
        if _pf_errors:
            return mo.Html('<div class="eq-wrap"><div class="eq-err"><strong>Cannot launch:</strong><br>'
                           + '<br>'.join(_pf_errors) + '</div></div>')

        # Assembly types only needed when there are read samples
        _types = []
        if _selected:
            if asm_s.value:
                _types.append("S")
            if asm_l.value:
                _types.append("L")
            if asm_hl.value:
                _types.append("HL")
            if asm_hs.value:
                _types.append("HS")
            if not _types:
                return mo.Html('<div class="eq-wrap"><div class="eq-err">Select at least one assembly type for read-based samples.</div></div>')

            # Validate assembly types against available reads in selected samples
            _has_long = any(_it.get('L_path') for _it in _selected)
            _has_short = any(_it.get('R1_path') and _it.get('R2_path') for _it in _selected)
            _need_long = any(t in ('L', 'HL') for t in _types)
            _need_short = any(t in ('S', 'HS') for t in _types)
            _need_hybrid = any(t in ('HL', 'HS') for t in _types)
            _mismatches = []
            if _need_long and not _has_long:
                _long_types = [t for t in _types if t in ('L', 'HL')]
                _mismatches.append(f"Assembly type(s) <strong>{', '.join(_long_types)}</strong> require long reads (ONT), but no <code>*_L.fastq.gz</code> files were found.")
            if _need_short and not _has_short:
                _short_types = [t for t in _types if t in ('S', 'HS')]
                _mismatches.append(f"Assembly type(s) <strong>{', '.join(_short_types)}</strong> require short reads (Illumina), but no <code>*_R1/*_R2.fastq.gz</code> pairs were found.")
            if _need_hybrid and not (_has_long and _has_short):
                _hyb_types = [t for t in _types if t in ('HL', 'HS')]
                if _has_long and not _has_short:
                    _mismatches.append(f"Hybrid assembly type(s) <strong>{', '.join(_hyb_types)}</strong> require both long and short reads, but no short read pairs were found.")
                elif _has_short and not _has_long:
                    _mismatches.append(f"Hybrid assembly type(s) <strong>{', '.join(_hyb_types)}</strong> require both long and short reads, but no long reads were found.")
            if _mismatches:
                _err_html = ''.join(f'<div class="eq-err">{m}</div>' for m in _mismatches)
                return mo.Html(f'<div class="eq-wrap">{_err_html}</div>')
        else:
            _types = ["L"]  # Dummy type for assembly-only runs

        _as = ",".join(_types)
        _cfg = load_config()
        _pl = _cfg.get("pipeline_dir", ".")
        # Resolve a relative pipeline_dir against this launcher's directory
        if not os.path.isabs(_pl):
            _pl = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), _pl))
        _ob = _cfg.get("output_base", "/data/epiquence/results")
        _gid = group_dropdown.value if group_dropdown.value else "GP001"
        _gdesc = _cfg.get("groups", {}).get(_gid, "")
        _wid = "W" + workflow_now().strftime("%m%d%y%H%M%S")
        _wd = os.path.join(_ob, _wid)
        _rd = os.path.join(_wd, "raw_reads")
        _od = os.path.join(_wd, "results")
        try:
            os.makedirs(_rd, exist_ok=True)
            os.makedirs(_od, exist_ok=True)
            _sr = []
            for _it in _selected:
                _o = _it['name']
                _n = name_inputs[_o].value.strip() or _o
                _lp = _r1p = _r2p = ""
                if _it.get('L_path'):
                    _dp = os.path.join(_rd, f"{_n}_L.fastq.gz")
                    os.symlink(_it['L_path'], _dp)
                    _lp = _dp
                if _it.get('R1_path'):
                    _dp = os.path.join(_rd, f"{_n}_R1.fastq.gz")
                    os.symlink(_it['R1_path'], _dp)
                    _r1p = _dp
                if _it.get('R2_path'):
                    _dp = os.path.join(_rd, f"{_n}_R2.fastq.gz")
                    os.symlink(_it['R2_path'], _dp)
                    _r2p = _dp
                _sr.append({'sample_id': _n, 'long_R': _lp, 'short_R1': _r1p, 'short_R2': _r2p})
            _sp = os.path.join(_wd, "sample_sheet.csv")
            with open(_sp, 'w', newline='') as _sf:
                _wr = csv.DictWriter(_sf, fieldnames=['sample_id', 'long_R', 'short_R1', 'short_R2'])
                _wr.writeheader()
                _wr.writerows(_sr)
            _cmd = ["nextflow", "run", os.path.join(_pl, "main.nf"), "--sample_sheet", _sp, "--outdir", _od, "--assembly_type", _as, "--workflow_id", _wid, "--workflow_description", workflow_desc.value or "", "--group_id", _gid, "--group_description", _gdesc, "-work-dir", os.path.join(_wd, "work"), "-bg"]
            if assemblies_dropdown.value:
                _cmd.extend(["--assemblies_dir", assemblies_dropdown.value])
            if skip_rasusa.value:
                _cmd.append("--skip_rasusa")
            else:
                _cmd.extend(["--rasusa_coverage", str(int(rasusa_coverage.value))])
                _cmd.extend(["--rasusa_genome_size", f"{rasusa_genome_size.value}mb"])
            if outbreak_check.value:
                _cmd.append("--outbreak")
                _cmd.extend(["--outbreak_threshold", str(mash_threshold.value)])
                if core_genome_check.value:
                    _cmd.append("--core_genome")
                if metadata_path.value and metadata_path.value.strip():
                    _cmd.extend(["--metadata", metadata_path.value.strip()])
            if hgt_check.value:
                _cmd.append("--plasmid_hgt")
                _cmd.extend(["--hgt_mash_threshold", str(hgt_threshold.value)])
                _cmd.extend(["--hgt_sketch_size", str(int(hgt_sketch.value))])
                if hgt_maps_check.value:
                    _cmd.append("--hgt_maps")
            if surveillance_db_check.value:
                _cmd.append("--surveillance")
            # Compute profile (Nextflow option, single dash) plus the engine it implies.
            _cmd.extend(["-profile", _prof, "--taxonomy_tool", _eng])
            _rs = os.path.join(_wd, "run.sh")
            with open(_rs, 'w') as _rf:
                _rf.write("#!/bin/bash\n")
                _rf.write(f"cd {_wd}\n")
                _rf.write(" ".join(f'"{_c}"' if ' ' in _c else _c for _c in _cmd) + "\n")
            os.chmod(_rs, 0o755)
            _logp = os.path.join(_wd, "nextflow_launch.log")
            with open(_logp, 'w') as _lf:
                subprocess.Popen(_cmd, cwd=_wd, stdout=_lf, stderr=subprocess.STDOUT, start_new_session=True)
            _det = [f"{len(_sr)} samples", f"{_as} assemblies", _gid,
                    f"profile: {_prof}", f"species: {_eng}"]
            if outbreak_check.value:
                _det.append("outbreak")
            if hgt_check.value:
                _det.append("plasmid HGT")
            if surveillance_db_check.value:
                _det.append("surveillance")
            _det_txt = " \u00b7 ".join(_det)
            return mo.Html(f'<div class="eq-wrap"><div class="eq-ok"><div style="font-size:40px;">\u2713</div><h3>Workflow Submitted Successfully</h3><p>Pipeline is running in the background. Nextflow will handle the rest.</p><div class="eq-ok-id">{_wid}</div><div class="eq-ok-det">{_det_txt}<br>{_wd}</div></div></div>')
        except Exception as _exc:
            return mo.Html(f'<div class="eq-wrap"><div class="eq-err"><strong>Error:</strong> {str(_exc)}</div></div>')
    _go()
    return


@app.cell
def _(mo, os, re, load_config):
    """Scan output_base for previous workflow directories."""
    _cfg = load_config()
    _ob = _cfg.get("output_base", "/data/epiquence/results")
    _wf_opts = {}
    if os.path.isdir(_ob):
        for _d in sorted(os.listdir(_ob), reverse=True):
            _full = os.path.join(_ob, _d)
            if os.path.isdir(_full) and _d.startswith('W'):
                _ss = os.path.join(_full, "sample_sheet.csv")
                _rs = os.path.join(_full, "run.sh")
                if os.path.exists(_ss) or os.path.exists(_rs):
                    # Try to get sample count and status
                    _info_parts = [_d]
                    if os.path.exists(_ss):
                        try:
                            with open(_ss) as _f:
                                _n_samples = sum(1 for _ in _f) - 1
                            _info_parts.append(f"{_n_samples} samples")
                        except Exception:
                            pass
                    # Check for run.sh to extract flags
                    if os.path.exists(_rs):
                        try:
                            with open(_rs) as _f:
                                _run_content = _f.read()
                            if '--outbreak' in _run_content:
                                _info_parts.append("outbreak")
                            if '--core_genome' in _run_content:
                                _info_parts.append("CGA")
                            if '--plasmid_hgt' in _run_content:
                                _info_parts.append("HGT")
                            if '--hgt_maps' in _run_content:
                                _info_parts.append("ORF Maps")
                            import re as _re
                            _pm = _re.search(r'(?<!\S)-profile\s+(\S+)', _run_content)
                            if _pm:
                                _info_parts.append(_pm.group(1))
                            _tm = _re.search(r'--taxonomy_tool\s+(\S+)', _run_content)
                            if _tm:
                                _info_parts.append(_tm.group(1))
                            _asm_match = _re.search(r'--assembly_type\s+(\S+)', _run_content)
                            if _asm_match:
                                _info_parts.append(_asm_match.group(1))
                            _desc_match = _re.search(r'--workflow_description\s+"([^"]*)"', _run_content)
                            if not _desc_match:
                                _desc_match = _re.search(r'--workflow_description\s+(\S+)', _run_content)
                            if _desc_match and _desc_match.group(1).strip():
                                _info_parts.append(_desc_match.group(1).strip())
                        except Exception:
                            pass
                    # Check if completed or still running
                    _report = os.path.join(_full, "results", "multiqc", "multiqc_report.html")
                    if os.path.exists(_report):
                        _info_parts.append("\u2705 complete")
                    elif os.path.exists(os.path.join(_full, ".nextflow.log")):
                        _info_parts.append("\u23F3 ran")
                    _label = " \u00b7 ".join(_info_parts)
                    _wf_opts[_label] = _full
    resume_dropdown = mo.ui.dropdown(options=_wf_opts, label="Select Workflow") if _wf_opts else None
    resume_btn = mo.ui.run_button(label="Resume", kind="warn") if _wf_opts else None
    resume_metadata = mo.ui.text(placeholder="Path to metadata CSV/TSV (optional, adds to -resume)", full_width=True) if _wf_opts else None
    return resume_dropdown, resume_btn, resume_metadata


@app.cell
def _(mo, os, re, resume_dropdown, config):
    """Parse selected workflow's run.sh to detect existing flags and thresholds."""
    _has_outbreak = False
    _has_hgt = False
    _has_hgt_maps = False
    _has_core_genome = False
    _has_surveillance = False
    _outbreak_thresh = 0.05
    _hgt_thresh = 0.005
    _profile = ""
    if resume_dropdown and resume_dropdown.value:
        _rs = os.path.join(resume_dropdown.value, "run.sh")
        if os.path.exists(_rs):
            try:
                with open(_rs) as _f:
                    _content = _f.read()
                _has_outbreak = '--outbreak' in _content
                _has_hgt = '--plasmid_hgt' in _content
                _has_hgt_maps = '--hgt_maps' in _content
                _has_core_genome = '--core_genome' in _content
                _has_surveillance = '--surveillance' in _content
                # Parse existing threshold values
                _m = re.search(r'--outbreak_threshold\s+([\d.]+)', _content)
                if _m:
                    _outbreak_thresh = float(_m.group(1))
                _m = re.search(r'--hgt_mash_threshold\s+([\d.]+)', _content)
                if _m:
                    _hgt_thresh = float(_m.group(1))
                _m = re.search(r'(?<!\S)-profile\s+(\S+)', _content)
                if _m:
                    _profile = _m.group(1)
            except Exception:
                pass
    # Unified toggle checkboxes: pre-checked based on current run state
    resume_toggle_outbreak = mo.ui.checkbox(label="Outbreak Detection", value=_has_outbreak) if resume_dropdown else None
    resume_toggle_core_genome = mo.ui.checkbox(label="Core Genome Analysis", value=_has_core_genome) if resume_dropdown else None
    resume_toggle_hgt = mo.ui.checkbox(label="Plasmid HGT Detection", value=_has_hgt) if resume_dropdown else None
    resume_toggle_hgt_maps = mo.ui.checkbox(label="ORF Comparisons", value=_has_hgt_maps) if resume_dropdown else None
    resume_toggle_surveillance = mo.ui.checkbox(label="Surveillance", value=_has_surveillance) if resume_dropdown else None
    # Threshold inputs: show when outbreak/hgt is active
    resume_outbreak_thresh = mo.ui.number(value=_outbreak_thresh, start=0.001, stop=0.1, step=0.001, label="Outbreak Mash Threshold") if resume_dropdown else None
    resume_hgt_thresh = mo.ui.number(value=_hgt_thresh, start=0.001, stop=0.05, step=0.001, label="HGT Mash Threshold") if resume_dropdown else None
    # Compute profile, preset from the run being resumed.
    _r_profile_opts = {
        "High capacity | 64+ threads, 256+ GB RAM": "high_capacity",
        "Low capacity | 8 to 32 threads, 16 to 64 GB RAM": "low_capacity",
    }
    _preset = _profile if _profile in _r_profile_opts.values() else (config.get("default_profile") or "high_capacity")
    resume_profile = mo.ui.dropdown(
        options=_r_profile_opts, label="Compute Profile", full_width=True,
        value=next((k for k, v in _r_profile_opts.items() if v == _preset), list(_r_profile_opts)[0])) if resume_dropdown is not None else None
    resume_had_profile = _profile
    resume_has_outbreak = _has_outbreak
    resume_has_hgt = _has_hgt
    resume_has_hgt_maps = _has_hgt_maps
    resume_has_core_genome = _has_core_genome
    resume_has_surveillance = _has_surveillance
    return resume_toggle_outbreak, resume_toggle_core_genome, resume_toggle_hgt, resume_toggle_hgt_maps, resume_toggle_surveillance, resume_has_outbreak, resume_has_hgt, resume_has_hgt_maps, resume_has_core_genome, resume_has_surveillance, resume_outbreak_thresh, resume_hgt_thresh, resume_profile, resume_had_profile


@app.cell
def _(mo, resume_dropdown, resume_btn, resume_metadata, resume_toggle_outbreak, resume_toggle_core_genome, resume_toggle_hgt, resume_toggle_hgt_maps, resume_toggle_surveillance, resume_has_outbreak, resume_has_hgt, resume_has_hgt_maps, resume_has_core_genome, resume_has_surveillance, resume_outbreak_thresh, resume_hgt_thresh, resume_profile, resume_had_profile, preflight):
    def _resume_card():
        if not resume_dropdown:
            return mo.Html("")
        _dd = mo.as_html(resume_dropdown).text
        _btn = mo.as_html(resume_btn).text
        _meta = mo.as_html(resume_metadata).text

        _toggle_items = []
        if resume_toggle_outbreak:
            _toggle_items.append(mo.as_html(resume_toggle_outbreak).text)
        if resume_toggle_core_genome:
            _toggle_items.append(mo.as_html(resume_toggle_core_genome).text)
        if resume_toggle_hgt:
            _toggle_items.append(mo.as_html(resume_toggle_hgt).text)
        if resume_toggle_hgt_maps:
            _toggle_items.append(mo.as_html(resume_toggle_hgt_maps).text)
        if resume_toggle_surveillance:
            _toggle_items.append(mo.as_html(resume_toggle_surveillance).text)
        _toggles_html = ""
        if _toggle_items:
            _toggles_html = f'<div style="display:flex;gap:16px;margin-top:8px;flex-wrap:wrap;">{" ".join(_toggle_items)}</div>'

        _thresh_html = ""
        _thresh_items = []
        if resume_outbreak_thresh:
            _thresh_items.append(mo.as_html(resume_outbreak_thresh).text)
        if resume_hgt_thresh:
            _thresh_items.append(mo.as_html(resume_hgt_thresh).text)
        if _thresh_items:
            _thresh_html = f'<div style="display:flex;gap:16px;margin-top:8px;">{" ".join(_thresh_items)}</div>'

        # Compute profile for the resumed run
        _compute_html = ""
        if resume_profile is not None and resume_dropdown.value:
            _p = resume_profile.value or "high_capacity"
            _msgs = []
            if _p != (resume_had_profile or ""):
                _msgs.append(("warn", "Changing the compute profile changes thread counts, the species engine and "
                                      "the Bakta image, so most steps re-run instead of coming from cache."))
            _errs, _warns = preflight(_p)
            _msgs = [("err", m) for m in _errs] + [("warn", m) for m in _warns] + _msgs
            _note_style = {"err": ("rgba(216,90,59,.12)", "\u26D4"), "warn": ("rgba(230,170,40,.14)", "\u26A0")}
            _notes = "".join(
                '<div style="margin-top:6px;padding:6px 10px;border-radius:6px;font-size:12px;'
                f'background:{_note_style[k][0]};">{_note_style[k][1]} {m}</div>' for k, m in _msgs)
            _compute_html = f"""<div style="margin-top:10px;">{mo.as_html(resume_profile).text}</div>
            {_notes}"""

        return mo.Html(f"""<div class="eq-wrap"><div class="eq-card">
        <div class="eq-title">\U0001F504 Resume Previous Workflow</div>
        <div class="eq-upload-hint" style="margin-bottom:8px;">Select a previous workflow to resume with <code>-resume</code>. Completed tasks will be cached.</div>
        {_dd}
        {_compute_html}
        {_toggles_html}
        {_thresh_html}
        <div style="margin-top:8px;">{_meta}</div>
        <div style="text-align:center;margin-top:12px;">{_btn}</div>
        </div></div>""")
    _resume_card()
    return


@app.cell
def _(mo, os, subprocess, resume_btn, resume_dropdown, resume_metadata, resume_toggle_outbreak, resume_toggle_core_genome, resume_toggle_hgt, resume_toggle_hgt_maps, resume_toggle_surveillance, resume_has_outbreak, resume_has_hgt, resume_has_hgt_maps, resume_has_core_genome, resume_has_surveillance, resume_outbreak_thresh, resume_hgt_thresh, resume_profile, profile_engine, preflight, load_config):
    def _do_resume():
        if not resume_btn or not resume_btn.value or not resume_dropdown or not resume_dropdown.value:
            return mo.Html("")
        _r_prof = (resume_profile.value or "high_capacity") if resume_profile is not None else "high_capacity"
        _r_eng = profile_engine(_r_prof)
        _pf_errors, _ = preflight(_r_prof)
        if _pf_errors:
            return mo.Html('<div class="eq-wrap"><div class="eq-err"><strong>Cannot resume:</strong><br>'
                           + '<br>'.join(_pf_errors) + '</div></div>')
        _wd = resume_dropdown.value
        _wid = os.path.basename(_wd)
        _rs = os.path.join(_wd, "run.sh")
        if not os.path.exists(_rs):
            return mo.Html(f'<div class="eq-wrap"><div class="eq-err">No run.sh found in {_wid}. Cannot resume.</div></div>')
        try:
            with open(_rs) as _f:
                _content = _f.read()
            # Extract the original nextflow command
            _lines = [l.strip() for l in _content.strip().split('\n') if l.strip() and not l.strip().startswith('#') and not l.strip().startswith('cd ')]
            if not _lines:
                return mo.Html(f'<div class="eq-wrap"><div class="eq-err">Could not parse run.sh in {_wid}.</div></div>')
            _orig_cmd = _lines[-1]
            # Remove -bg flag: subprocess.Popen already runs in background
            _orig_cmd = _orig_cmd.replace(' -bg', '')
            # Clean up deprecated params from old run.sh files
            import re as _re_clean
            _orig_cmd = _re_clean.sub(r"\s*--store_db_url\s+'[^']*'", '', _orig_cmd)
            _orig_cmd = _re_clean.sub(r'\s*--store_db_url\s+\S+', '', _orig_cmd)
            _orig_cmd = _re_clean.sub(r'\s*--store_db\b', '', _orig_cmd)
            _orig_cmd = _re_clean.sub(r"\s*--surveillance_db_url_input\s+'[^']*'", '', _orig_cmd)
            _orig_cmd = _re_clean.sub(r'\s*--surveillance_db_url_input\s+\S+', '', _orig_cmd)
            _orig_cmd = _re_clean.sub(r"\s*--surveillance_db_url\s+'[^']*'", '', _orig_cmd)
            _orig_cmd = _re_clean.sub(r'\s*--surveillance_db_url\s+\S+', '', _orig_cmd)
            # Update pipeline path to current version
            _cfg = load_config()
            _pl = _cfg.get("pipeline_dir", ".")
            if not os.path.isabs(_pl):
                _pl = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), _pl))
            _current_pl = os.path.join(_pl, "main.nf")
            import re as _re_pipe
            _orig_cmd = _re_pipe.sub(r'nextflow\s+run\s+\S+/main\.nf', f'nextflow run {_current_pl}', _orig_cmd)
            # Add -resume flag if not already present
            if '-resume' not in _orig_cmd:
                _orig_cmd += ' -resume'
            # Add metadata if provided
            _meta_val = resume_metadata.value.strip() if resume_metadata and resume_metadata.value else ""
            if _meta_val:
                if '--metadata' in _orig_cmd:
                    # Replace existing metadata path
                    import re as _re
                    _orig_cmd = _re.sub(r'--metadata\s+\S+', f'--metadata {_meta_val}', _orig_cmd)
                else:
                    _orig_cmd += f' --metadata {_meta_val}'

            # Apply module toggles: strip all module flags first, then add back what's checked
            import re as _re_mod
            # Strip all module flags
            _orig_cmd = _re_mod.sub(r'\s*--outbreak\b', '', _orig_cmd)
            _orig_cmd = _re_mod.sub(r'\s*--outbreak_threshold\s+\S+', '', _orig_cmd)
            _orig_cmd = _re_mod.sub(r'\s*--core_genome\b', '', _orig_cmd)
            _orig_cmd = _re_mod.sub(r'\s*--plasmid_hgt\b', '', _orig_cmd)
            _orig_cmd = _re_mod.sub(r'\s*--hgt_mash_threshold\s+\S+', '', _orig_cmd)
            _orig_cmd = _re_mod.sub(r'\s*--hgt_sketch_size\s+\S+', '', _orig_cmd)
            _orig_cmd = _re_mod.sub(r'\s*--hgt_maps\b', '', _orig_cmd)
            _orig_cmd = _re_mod.sub(r'\s*--surveillance\b', '', _orig_cmd)
            _orig_cmd = _re_mod.sub(r'\s*--surveillance_db_url\s+\S+', '', _orig_cmd)
            _orig_cmd = _re_mod.sub(r"\s*--surveillance_db_url\s+'[^']*'", '', _orig_cmd)
            # Compute profile and the engine it implies are re-applied from the control, replacing
            # whatever the original run used.
            _orig_cmd = _re_mod.sub(r'(?<!\S)-profile\s+\S+', '', _orig_cmd)
            _orig_cmd = _re_mod.sub(r'\s*--taxonomy_tool\s+\S+', '', _orig_cmd)
            _orig_cmd += f' -profile {_r_prof} --taxonomy_tool {_r_eng}'

            # Add back what's checked
            _want_outbreak = resume_toggle_outbreak and resume_toggle_outbreak.value
            _want_cga = resume_toggle_core_genome and resume_toggle_core_genome.value
            _want_hgt = resume_toggle_hgt and resume_toggle_hgt.value
            _want_maps = resume_toggle_hgt_maps and resume_toggle_hgt_maps.value
            _want_surv = resume_toggle_surveillance and resume_toggle_surveillance.value

            if _want_outbreak:
                _orig_cmd += ' --outbreak'
                if resume_outbreak_thresh:
                    _orig_cmd += f' --outbreak_threshold {resume_outbreak_thresh.value}'
                if _want_cga:
                    _orig_cmd += ' --core_genome'

            if _want_hgt:
                _orig_cmd += ' --plasmid_hgt'
                if resume_hgt_thresh:
                    _orig_cmd += f' --hgt_mash_threshold {resume_hgt_thresh.value}'
                if _want_maps:
                    _orig_cmd += ' --hgt_maps'

            if _want_surv:
                _orig_cmd += ' --surveillance'

            # Add metadata if provided
            _meta_val = resume_metadata.value.strip() if resume_metadata and resume_metadata.value else ""
            if _meta_val:
                if '--metadata' in _orig_cmd:
                    _orig_cmd = _re_mod.sub(r'--metadata\s+\S+', f'--metadata {_meta_val}', _orig_cmd)
                else:
                    _orig_cmd += f' --metadata {_meta_val}'
                # Ensure outbreak is enabled if metadata is added
                if '--outbreak' not in _orig_cmd:
                    _orig_cmd += ' --outbreak'

            # Add --surveillance if checkbox ticked and not already present
            if _want_surv and '--surveillance' not in _orig_cmd:
                _orig_cmd += ' --surveillance'

            # Thresholds are already included above via the toggle logic.
            import re as _re_thresh

            # Determine which processes need clearing for re-run
            _processes_to_clear = ['MULTIQC', 'GENERATE_MULTIQC_CONFIG']

            # Outbreak toggled on (was off)
            if _want_outbreak and not resume_has_outbreak:
                _processes_to_clear.extend(['RENDER_TREES', 'BUILD_OUTBREAK_REPORT'])
            # Outbreak toggled off (was on): clear outbreak caches
            if not _want_outbreak and resume_has_outbreak:
                _processes_to_clear.extend(['OUTBREAK_CLUSTER', 'SKA_ALIGNMENT', 'SNP_DISTS_CORE', 'IQTREE', 'ASSEMBLE_OUTBREAK_RESULTS', 'RENDER_TREES', 'BUILD_OUTBREAK_REPORT'])

            # Core genome toggled on (was off)
            if _want_cga and not resume_has_core_genome:
                _processes_to_clear.extend(['BAKTA_CGA', 'PANAROO', 'SNP_DISTS_CGA', 'IQTREE_CGA', 'ASSEMBLE_CGA_RESULTS', 'RENDER_CGA_TREES', 'BUILD_OUTBREAK_REPORT'])

            # HGT toggled on (was off)
            if _want_hgt and not resume_has_hgt:
                _processes_to_clear.extend(['COLLECT_PLASMID_FASTAS', 'MASH_DIST_PLASMIDS', 'PLASMID_ANI_AF', 'ASSEMBLE_HGT_RESULTS', 'BUILD_HGT_REPORT'])
            # HGT toggled off (was on): clear HGT caches
            if not _want_hgt and resume_has_hgt:
                _processes_to_clear.extend(['COLLECT_PLASMID_FASTAS', 'MASH_DIST_PLASMIDS', 'PLASMID_ANI_AF', 'ASSEMBLE_HGT_RESULTS', 'BUILD_HGT_REPORT', 'RENDER_HGT_TREES', 'UMAP_PLASMIDS'])

            # HGT maps toggled on (was off)
            if _want_maps and not resume_has_hgt_maps:
                _processes_to_clear.extend(['COMPARISON_MAP', 'BUILD_HGT_REPORT'])

            # Clear work dirs for processes that need re-running
            _work_dir = os.path.join(_wd, "work")
            if os.path.isdir(_work_dir) and len(_processes_to_clear) > 1:
                _nf_log = os.path.join(_wd, ".nextflow.log")
                if os.path.exists(_nf_log):
                    import re as _re2
                    _dirs_to_clear = set()
                    with open(_nf_log) as _lf:
                        for _line in _lf:
                            if any(p in _line for p in _processes_to_clear):
                                _wdir_match = _re2.search(r'work/([0-9a-f]{2}/[0-9a-f]+)', _line)
                                if _wdir_match:
                                    _dirs_to_clear.add(os.path.join(_work_dir, _wdir_match.group(1)))
                    import shutil as _shutil
                    for _d in _dirs_to_clear:
                        if os.path.isdir(_d):
                            _shutil.rmtree(_d)
            _resume_sh = os.path.join(_wd, "run_resume.sh")
            with open(_resume_sh, 'w') as _f:
                _f.write("#!/bin/bash\n")
                _f.write(f"cd {_wd}\n")
                _f.write(_orig_cmd + "\n")
            os.chmod(_resume_sh, 0o755)
            # Also update run.sh so flags are detected on next resume
            _rs_update = os.path.join(_wd, "run.sh")
            with open(_rs_update, 'w') as _f:
                _f.write("#!/bin/bash\n")
                _f.write(f"cd {_wd}\n")
                _f.write(_orig_cmd + "\n")
            os.chmod(_rs_update, 0o755)
            _logp = os.path.join(_wd, "nextflow_resume.log")
            import shlex
            _cmd_parts = shlex.split(_orig_cmd)
            with open(_logp, 'w') as _lf:
                subprocess.Popen(_cmd_parts, cwd=_wd, stdout=_lf, stderr=subprocess.STDOUT, start_new_session=True)
            _enabled = []
            _disabled = []
            if _want_outbreak:
                _enabled.append("Outbreak")
            elif resume_has_outbreak:
                _disabled.append("Outbreak")
            if _want_cga:
                _enabled.append("Core Genome")
            elif resume_has_core_genome:
                _disabled.append("Core Genome")
            if _want_hgt:
                _enabled.append("HGT")
            elif resume_has_hgt:
                _disabled.append("HGT")
            if _want_maps:
                _enabled.append("ORF Comparisons")
            elif resume_has_hgt_maps:
                _disabled.append("ORF Comparisons")
            if _want_surv:
                _enabled.append("Surveillance")
            elif resume_has_surveillance:
                _disabled.append("Surveillance")
            if _meta_val:
                _enabled.append(f"Metadata: {_meta_val}")
            _enabled.append(f"profile: {_r_prof}")
            _enabled.append(f"species: {_r_eng}")
            _detail_lines = ""
            if _enabled:
                _detail_lines += f"<br>Enabled: {', '.join(_enabled)}"
            if _disabled:
                _detail_lines += f"<br>Disabled: {', '.join(_disabled)}"
            return mo.Html(f'<div class="eq-wrap"><div class="eq-ok"><div style="font-size:40px;">\U0001F504</div><h3>Workflow Resumed</h3><p>Resuming with <code>-resume</code>. Cached tasks will be reused.{_detail_lines}</p><div class="eq-ok-id">{_wid}</div><div class="eq-ok-det">{_wd}</div></div></div>')
        except Exception as _exc:
            return mo.Html(f'<div class="eq-wrap"><div class="eq-err"><strong>Error:</strong> {str(_exc)}</div></div>')
    _do_resume()
    return


@app.cell
def _(mo):
    mo.Html('<div class="eq-wrap"><div class="eq-foot">EPIquence \u00b7 Sequence to Surveillance</div></div>')
    return


if __name__ == "__main__":
    app.run()