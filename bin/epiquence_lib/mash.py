"""Mash distance-matrix (``mash dist -t``) parsing."""

import csv
import os


def read_matrix(path, basename=False, strip_fasta=False):
    """Read a Mash ``-t`` distance matrix."""
    def clean(name):
        if basename:
            name = os.path.basename(name)
        if strip_fasta:
            name = name.replace(".fasta", "")
        return name.strip()

    row_names, col_names, rows = [], [], []
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        col_names = [clean(h) for h in header[1:]]
        for row in reader:
            row_names.append(clean(row[0]))
            parsed = []
            for x in row[1:]:
                try:
                    parsed.append(float(x))
                except ValueError:
                    parsed.append(1.0)
            rows.append(parsed)
    return row_names, col_names, rows


def pairwise_dict(row_names, col_names, rows, symmetric=True):
    """Build ``{(a, b): distance}`` from a parsed matrix."""
    dist = {}
    for i, a in enumerate(row_names):
        for j, b in enumerate(col_names):
            if j >= len(rows[i]):
                continue
            d = rows[i][j]
            dist[(a, b)] = d
            if symmetric:
                dist[(b, a)] = d
    return dist