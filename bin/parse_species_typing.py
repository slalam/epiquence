#!/usr/bin/env python3
"""Parse species-based typing results (ECTyper, Kleborate, emmtyper)."""
import csv
import os
import json
from collections import defaultdict


def parse_ectyper_dirs(ectyper_dirs_path):
    """Parse ECTyper output directories."""
    results = {}
    if not os.path.exists(ectyper_dirs_path):
        return results

    for d in sorted(os.listdir(ectyper_dirs_path)):
        ec_dir = os.path.join(ectyper_dirs_path, d)
        output_file = os.path.join(ec_dir, 'output.tsv')
        if not os.path.isdir(ec_dir) or not os.path.exists(output_file):
            continue
        sample_name = d.replace('_ectyper', '')
        with open(output_file) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                results[sample_name] = dict(row)
                break  # One row per sample
    return results


def parse_kleborate_dirs(kleborate_dirs_path):
    """Parse Kleborate output directories."""
    results = {}
    if not os.path.exists(kleborate_dirs_path):
        return results

    for d in sorted(os.listdir(kleborate_dirs_path)):
        kleb_dir = os.path.join(kleborate_dirs_path, d)
        if not os.path.isdir(kleb_dir):
            continue
        sample_name = d.replace('_kleborate', '')

        # Kleborate v3 outputs: klebsiella_pneumo_complex_output.txt or
        # klebsiella_oxytoca_complex_output.txt
        for fname in ['klebsiella_pneumo_complex_output.txt', 'klebsiella_oxytoca_complex_output.txt']:
            output_file = os.path.join(kleb_dir, fname)
            if os.path.exists(output_file):
                preset = 'kpsc' if 'pneumo' in fname else 'kosc'
                with open(output_file) as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        row_dict = dict(row)
                        row_dict['_preset'] = preset
                        results[sample_name] = row_dict
                        break
    return results


def parse_emmtyper_results(emmtyper_dir_path):
    """Parse emmtyper TSV results (no headers in output)."""
    results = {}
    if not os.path.exists(emmtyper_dir_path):
        return results

    for fname in sorted(os.listdir(emmtyper_dir_path)):
        if not fname.endswith('_emmtyper.tsv'):
            continue
        filepath = os.path.join(emmtyper_dir_path, fname)
        sample_name = fname.replace('_emmtyper.tsv', '')
        with open(filepath) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    results[sample_name] = {
                        'filename': parts[0],
                        'num_hits': parts[1],
                        'num_clusters': parts[2],
                        'emm_type': parts[3],
                        'position_1': parts[4],
                        'possible_emm_2': parts[5],
                        'position_2': parts[6],
                        'emm_cluster': parts[7],
                    }
                elif len(parts) >= 4:
                    results[sample_name] = {
                        'filename': parts[0],
                        'num_hits': parts[1] if len(parts) > 1 else '',
                        'num_clusters': parts[2] if len(parts) > 2 else '',
                        'emm_type': parts[3] if len(parts) > 3 else '',
                        'position_1': parts[4] if len(parts) > 4 else '',
                        'possible_emm_2': parts[5] if len(parts) > 5 else '',
                        'position_2': parts[6] if len(parts) > 6 else '',
                        'emm_cluster': parts[7] if len(parts) > 7 else '',
                    }
                break  # One line per sample
    return results

def _first_row(path, delimiter):
    """Return the first data row of a delimited file with a header, as a dict."""
    with open(path) as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            return dict(row)
    return {}


def parse_sistr_dirs(sistr_dirs_path):
    """Parse SISTR output directories (Salmonella)."""
    results = {}
    if not os.path.exists(sistr_dirs_path):
        return results
    for d in sorted(os.listdir(sistr_dirs_path)):
        sistr_dir = os.path.join(sistr_dirs_path, d)
        if not os.path.isdir(sistr_dir):
            continue
        sample_name = d.replace('_sistr', '')
        out_csv = os.path.join(sistr_dir, 'sistr-output.csv')
        if os.path.exists(out_csv):
            row = _first_row(out_csv, ',')
            if row:
                results[sample_name] = row
    return results


def parse_seroba_dirs(seroba_dirs_path):
    """Parse SeroBA output directories (Streptococcus pneumoniae)."""
    results = {}
    if not os.path.exists(seroba_dirs_path):
        return results
    for d in sorted(os.listdir(seroba_dirs_path)):
        seroba_dir = os.path.join(seroba_dirs_path, d)
        if not os.path.isdir(seroba_dir):
            continue
        sample_name = d.replace('_seroba', '')
        pred = None
        for fname in ('pred.csv', 'pred.tsv'):
            p = os.path.join(seroba_dir, fname)
            if os.path.exists(p):
                pred = p
                break
        if not pred:
            continue
        # Sniff delimiter from the header line
        with open(pred) as fh:
            header = fh.readline()
        delim = '\t' if header.count('\t') >= header.count(',') else ','
        row = _first_row(pred, delim)
        if row:
            results[sample_name] = row
    return results


def parse_meningotype_results(meningotype_dir_path):
    """Parse meningotype --finetype TSV results (Neisseria meningitidis)."""
    results = {}
    if not os.path.exists(meningotype_dir_path):
        return results
    for fname in sorted(os.listdir(meningotype_dir_path)):
        if not fname.endswith('_meningotype.tsv'):
            continue
        filepath = os.path.join(meningotype_dir_path, fname)
        sample_name = fname.replace('_meningotype.tsv', '')
        row = _first_row(filepath, '\t')
        if row:
            results[sample_name] = row
    return results


def parse_elgato_dirs(elgato_dirs_path):
    """Parse el_gato output directories (Legionella pneumophila)."""
    results = {}
    if not os.path.exists(elgato_dirs_path):
        return results
    for d in sorted(os.listdir(elgato_dirs_path)):
        elgato_dir = os.path.join(elgato_dirs_path, d)
        if not os.path.isdir(elgato_dir):
            continue
        sample_name = d.replace('_elgato', '')
        out_txt = os.path.join(elgato_dir, 'possible_mlsts.txt')
        if os.path.exists(out_txt):
            row = _first_row(out_txt, '\t')
            if row:
                results[sample_name] = row
    return results