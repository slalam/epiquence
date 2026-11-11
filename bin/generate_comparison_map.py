#!/usr/bin/env python3
"""Generates JSON from Bakta annotations and BLAST ORF comparisons, then renders SVG using the actual CGView.js library (server-side via Node.js)."""

import argparse
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path



def parse_fasta(path):
    seqs = {}
    circularity = {}  # contig_name -> bool
    hdr = None
    buf = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if hdr:
                seqs[hdr] = "".join(buf)
            full_header = line[1:]
            hdr = full_header.split()[0]
            buf = []
            # Check for circularity markers in full header
            import re
            is_circ = bool(re.search(r'circular[= _]true|topology[= _]circular|\[circularised\]|circular_contig|rotated[= _][Tt]rue|\bcircular\b', full_header, re.IGNORECASE))
            circularity[hdr] = is_circ
        else:
            buf.append(line)
    if hdr:
        seqs[hdr] = "".join(buf)
    return seqs, circularity

def parse_gff3(path):
    features = []
    keep = {"CDS","tRNA","rRNA","tmRNA","ncRNA","regulatory_region",
            "oriT","oriC","ncRNA_region","repeat_region"}
    for line in open(path):
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 9 or p[2] not in keep: continue
        attrs = {}
        for a in p[8].split(";"):
            if "=" in a:
                k, v = a.split("=", 1)
                attrs[k] = v
        name = attrs.get("gene", attrs.get("Name", attrs.get("ID", "unknown")))
        features.append({
            "contig": p[0],  # GFF3 column 1 = seqid/contig name
            "type": p[2], "start": int(p[3]), "end": int(p[4]),
            "strand": 1 if p[6] == "+" else -1,
            "name": name, "product": attrs.get("product", ""),
            "locus_tag": attrs.get("locus_tag", ""), "id": attrs.get("ID", ""),
        })
    return features

def parse_blast6(path):
    hits = []
    for line in open(path):
        p = line.strip().split("\t")
        if len(p) < 12: continue
        hits.append({
            "qseqid": p[0], "sseqid": p[1], "pident": float(p[2]),
            "length": int(p[3]), "mismatch": int(p[4]), "gapopen": int(p[5]),
            "qstart": int(p[6]), "qend": int(p[7]),
            "sstart": int(p[8]), "send": int(p[9]),
            "evalue": float(p[10]), "bitscore": float(p[11]),
        })
    return hits



def gc_content_plot(seq, window=500, step=250):
    data = []
    s = seq.upper()
    for i in range(0, len(s), step):
        chunk = s[i:i+window]
        if not chunk: continue
        gc = (chunk.count("G") + chunk.count("C")) / len(chunk)
        data.append({"start": i+1, "stop": min(i+window, len(s)), "score": gc})
    return data

def gc_skew_plot(seq, window=500, step=250):
    data = []
    s = seq.upper()
    for i in range(0, len(s), step):
        chunk = s[i:i+window]
        g, c = chunk.count("G"), chunk.count("C")
        skew = (g-c)/(g+c) if (g+c) else 0
        data.append({"start": i+1, "stop": min(i+window, len(s)), "score": skew})
    return data


# BAKTA / BLAST RUNNERS

def find_tool(name, custom_path=None):
    """Find an executable, checking custom path first, then PATH, then common conda locations."""
    import shutil
    if custom_path:
        if os.path.isfile(custom_path) and os.access(custom_path, os.X_OK):
            return custom_path
        # Maybe it's a directory containing the tool
        candidate = os.path.join(custom_path, name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
        raise FileNotFoundError(
            f"Tool '{name}' not found at specified path: {custom_path}"
        )
    found = shutil.which(name)
    if found:
        return found
    # Check common conda/mamba locations
    home = os.path.expanduser("~")
    for base in [
        os.path.join(home, "miniforge3"),
        os.path.join(home, "miniconda3"),
        os.path.join(home, "anaconda3"),
        os.path.join(home, "mambaforge"),
    ]:
        for env in ["envs/bakta/bin", "envs/blast/bin", "bin"]:
            candidate = os.path.join(base, env, name)
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                return candidate
    return None

def run_bakta(fasta, outdir, db, prefix=None, threads=4, bakta_path=None):
    os.makedirs(outdir, exist_ok=True)
    bakta_bin = find_tool("bakta", bakta_path)
    if not bakta_bin:
        raise FileNotFoundError(
            "bakta not found in PATH.\n"
            "  Options:\n"
            "    1. Activate your bakta conda env:  conda activate bakta\n"
            "    2. Specify path:  --bakta-path /path/to/bakta\n"
            "    3. Skip annotation:  --skip-bakta --ref-dir <bakta_output_dir>"
        )
    cmd = [bakta_bin,"--db",db,"--output",outdir,"--threads",str(threads),"--force"]
    if prefix: cmd += ["--prefix", prefix]
    cmd.append(str(fasta))
    print(f"[BAKTA] {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        print(f"[BAKTA] STDOUT: {r.stdout[-500:]}" if r.stdout else "")
        print(f"[BAKTA] STDERR: {r.stderr[-500:]}" if r.stderr else "")
        raise RuntimeError(f"Bakta failed for {fasta}")
    print(f"[BAKTA] Done: {fasta}")
    return outdir

def find_bakta(d):
    d = Path(d)
    gff = list(d.glob("*.gff3"))
    faa = list(d.glob("*.faa"))
    fna = list(d.glob("*.fna"))
    if not gff: raise FileNotFoundError(f"No .gff3 in {d}")
    if not faa: raise FileNotFoundError(f"No .faa in {d}")
    return str(gff[0]), str(faa[0]), str(fna[0]) if fna else None

def run_blast(ref_orfs_faa, comp_orfs_faa, out, evalue=1e-4, threads=4, blast_path=None):
    """Run BLAST matching ORF-vs-ORF approach: - makeblastdb on COMPARISON ORFs (subject/db) -
    blastp with REFERENCE ORFs as query - -seg yes for low-complexity filtering -
    -num_alignments 10000 (not max_target_seqs 1)"""
    makeblastdb_bin = find_tool("makeblastdb", blast_path)
    blastp_bin = find_tool("blastp", blast_path)
    if not makeblastdb_bin or not blastp_bin:
        raise FileNotFoundError(
            "BLAST+ (makeblastdb/blastp) not found in PATH.\n"
            "  Options:\n"
            "    1. Install:  conda install -c bioconda blast\n"
            "    2. Specify path:  --blast-path /path/to/blast/bin/\n"
            "    3. Skip BLAST:  --skip-blast --blast-files <file1> <file2>"
        )
    db = out + "_db"
    # DB is the COMPARISON (subject): standard
    r = subprocess.run([makeblastdb_bin, "-in", comp_orfs_faa, "-dbtype", "prot", "-out", db],
                       capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"makeblastdb failed: {r.stderr[:500]}")
    # Query is the REFERENCE ORFs: standard
    r = subprocess.run([blastp_bin,
                        "-query", ref_orfs_faa,
                        "-db", db,
                        "-outfmt", "6",
                        "-evalue", str(evalue),
                        "-num_threads", str(threads),
                        "-num_alignments", "10000",
                        "-seg", "yes",
                        "-out", out],
                       capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"blastp failed: {r.stderr[:500]}")
    print(f"[BLAST] Done: {Path(ref_orfs_faa).stem} vs {Path(comp_orfs_faa).stem}")
    return out


def extract_orfs(fasta, outdir, min_len=30, genetic_code=11):
    """Extract ORFs from all 6 reading frames, matching ORF-vs-ORF approach."""
    os.makedirs(outdir, exist_ok=True)
    prefix = Path(fasta).stem
    faa_path = os.path.join(outdir, f"{prefix}_orfs.faa")

    seqs, _ = parse_fasta(fasta)

    codon_table = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L',
        'CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S',
        'TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A',
        'GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N',
        'AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R',
        'CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }

    # Bacterial genetic code 11 start codons
    start_codons = {'ATG', 'GTG', 'TTG'}

    def revcomp(s):
        rc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
        return ''.join(rc.get(c, 'N') for c in reversed(s.upper()))

    def find_orfs_in_frame(dna, strand, frame_offset, seq_len):
        """Find all ORFs (start...stop) in a single reading frame."""
        orfs = []
        i = frame_offset
        while i < len(dna) - 2:
            codon = dna[i:i+3].upper()
            if codon in start_codons:
                # First aa is always M regardless of actual start codon (bacterial convention)
                protein = ['M']
                j = i + 3
                while j < len(dna) - 2:
                    c = dna[j:j+3].upper()
                    aa = codon_table.get(c, 'X')
                    if aa == '*':
                        j += 3
                        break
                    protein.append(aa)
                    j += 3
                prot_seq = ''.join(protein)
                if len(prot_seq) >= min_len:
                    if strand == 1:
                        nt_start = i + 1
                        nt_end = min(j, len(dna))
                    else:
                        nt_start = seq_len - j + 1
                        nt_end = seq_len - i
                    rf = frame_offset + 1
                    orfs.append((nt_start, nt_end, prot_seq, rf))
                i = j  # skip past this ORF
            else:
                i += 3
        return orfs

    all_orfs = []
    for contig_name, seq in seqs.items():
        dna = seq.upper()
        sl = len(dna)
        rc_dna = revcomp(dna)
        for frame in range(3):
            for start, end, prot, rf in find_orfs_in_frame(dna, 1, frame, sl):
                all_orfs.append((contig_name, start, end, 1, rf, prot))
            for start, end, prot, rf in find_orfs_in_frame(rc_dna, -1, frame, sl):
                all_orfs.append((contig_name, start, end, -1, rf, prot))

    # Filter enclosed ORFs: remove smaller ORFs fully contained within larger ones.
    all_orfs.sort(key=lambda x: (x[0], -(x[2] - x[1])))
    filtered = []
    for orf in all_orfs:
        contig, start, end, strand, rf, prot = orf
        enclosed = False
        for fc, fs, fe, fstrand, _, _ in filtered:
            if fc == contig and fs <= start and fe >= end:
                enclosed = True
                break
        if not enclosed:
            filtered.append(orf)

    removed = len(all_orfs) - len(filtered)
    all_orfs = filtered

    with open(faa_path, 'w') as f:
        for idx, (contig, start, end, strand, rf, prot) in enumerate(all_orfs, 1):
            header = f"{contig};orf_{idx};start={start};end={end};strand={strand};rf={rf}"
            f.write(f">{header}\n{prot}\n")

    print(f"[ORF] {prefix}: {len(all_orfs)} ORFs (6-frame, min {min_len} aa, "
          f"starts=ATG/GTG/TTG, {removed} enclosed removed)")
    return faa_path





AMR_PREFIXES = {"bla","aac","aph","tet","mcr","qnr","sul","dfr",
                "erm","mef","van","oxa","ndm","vim","imp","kpc","ctx",
                "shv","tem","cmx","cat","fos","str","arr"}
# ant(2''), ant(3'') etc.: require '(' to avoid matching 'antirestriction' etc.
AMR_PATTERNS_WITH_PAREN = {"ant"}

def is_amr(name):
    nl = name.lower()
    if any(nl.startswith(p) for p in AMR_PREFIXES):
        return True
    if any(nl.startswith(p + "(") or nl.startswith(p + "'") for p in AMR_PATTERNS_WITH_PAREN):
        return True
    return False


# IDENTITY → COLOR (Dark-blue → cyan → green)

def identity_color(pident):
    if pident >= 100: return "rgba(10,10,58,1)"
    if pident >= 99:  return "rgba(13,27,107,1)"
    if pident >= 95:  return "rgba(26,58,138,1)"
    if pident >= 90:  return "rgba(32,96,170,1)"
    if pident >= 80:  return "rgba(40,136,187,1)"
    if pident >= 70:  return "rgba(48,170,204,1)"
    if pident >= 60:  return "rgba(69,204,204,1)"
    if pident >= 50:  return "rgba(102,221,187,1)"
    if pident >= 40:  return "rgba(136,238,170,1)"
    return "rgba(187,244,204,1)"



def build_cgview_json(seq_len, ref_contigs, features, blast_comps,
                      gc_data=None, skew_data=None, display_name=None,
                      identity_cutoffs=None, circularity=None, amr_gene_data=None,
                      fasta_to_gff=None, ref_fasta=None, color_palette='warm'):
    """Build a CGView.js JSON in CGView.js JSON format."""
    if circularity is None:
        circularity = {}
    if amr_gene_data is None:
        amr_gene_data = {}
    if fasta_to_gff is None:
        fasta_to_gff = {}
    if ref_fasta is None:
        ref_fasta = ""
    if display_name is None:
        display_name = "Plasmid Map"

    if identity_cutoffs is None:
        identity_cutoffs = [0, 85, 90, 95, 100]
    identity_cutoffs = sorted(set(identity_cutoffs))
    # Ensure 0 and 100 are always present
    if identity_cutoffs[0] != 0:
        identity_cutoffs.insert(0, 0)
    if identity_cutoffs[-1] != 100:
        identity_cutoffs.append(100)

    # ── Generate BLAST identity color gradient ──
    palettes = {
        'warm': [
            (0.0,  (230, 218, 200)),  # warm sand
            (0.25, (210, 180, 130)),  # warm gold
            (0.50, (120, 185, 170)),  # sage teal
            (0.75, (65, 145, 155)),   # ocean teal
            (1.0,  (35, 80, 110)),    # deep steel blue
        ],
        'cool': [
            (0.0,  (210, 230, 240)),  # ice blue
            (0.25, (160, 200, 220)),  # light steel
            (0.50, (80, 150, 190)),   # medium blue
            (0.75, (40, 100, 160)),   # royal blue
            (1.0,  (15, 50, 100)),    # deep navy
        ],
        'earth': [
            (0.0,  (235, 225, 210)),  # cream
            (0.25, (195, 180, 155)),  # tan
            (0.50, (140, 155, 130)),  # sage green
            (0.75, (100, 120, 95)),   # olive
            (1.0,  (55, 70, 55)),     # dark forest
        ],
    }
    palette_name = color_palette if color_palette in palettes else 'warm'
    stops = palettes[palette_name]

    def blast_color(t):
        """Interpolate identity colormap."""
        for i in range(len(stops) - 1):
            t0, c0 = stops[i]
            t1, c1 = stops[i + 1]
            if t0 <= t <= t1:
                f = (t - t0) / (t1 - t0) if t1 > t0 else 0
                r = int(c0[0] + f * (c1[0] - c0[0]))
                g = int(c0[1] + f * (c1[1] - c0[1]))
                b = int(c0[2] + f * (c1[2] - c0[2]))
                return f"rgba({r},{g},{b},1)"
        return f"rgba({stops[-1][1][0]},{stops[-1][1][1]},{stops[-1][1][2]},1)"

    n_cutoffs = len(identity_cutoffs)
    blast_legend_items = []
    for i, cutoff in enumerate(identity_cutoffs):
        # Map cutoff index to color: 0th = lightest, last = darkest
        t = i / max(1, n_cutoffs - 1)
        color = blast_color(t)
        if cutoff == 100:
            name = "BLAST = 100% Identical"
        elif cutoff == 0:
            next_cutoff = identity_cutoffs[1] if len(identity_cutoffs) > 1 else 100
            name = f"BLAST < {next_cutoff}% Identical"
        else:
            name = f"BLAST >= {cutoff}% Identical"
        blast_legend_items.append({
            "name": name,
            "swatchColor": color,
            "decoration": "arc",
            "meta": {"blast_meta": {"blast_formatter": {"cutoff": cutoff}}},
        })

    ref_lookup = {}
    for f in features:
        for k in [f["locus_tag"], f["id"], f["name"]]:
            if k: ref_lookup[k] = f

    # ── Detect AMR genes: no longer colored separately in Bakta track ──

    # ── Legend items ──
    legend_items = [
        {"name": "CDS", "swatchColor": "rgba(176,93,55,1)", "decoration": "arrow"},
        {"name": "oriT", "swatchColor": "rgba(72,191,184,1)", "decoration": "arrow"},
        {"name": "oriC", "swatchColor": "rgba(219,80,146,1)", "decoration": "arrow"},
    ]
    legend_items.extend(blast_legend_items)

    # ── Bakta features ──
    bakta_source = "bakta-1.1"
    cgv_features = []

    for f in features:
        ftype = f["type"]
        name = f["name"]
        product = f["product"]

        if name and name != "unknown":
            label = name
        elif product:
            label = product[:50]
        else:
            label = "hypothetical protein"

        # Determine type and legend
        if ftype == "CDS":
            cgv_type = "cds"
            legend_name = "CDS"
        elif ftype == "oriT" or (ftype == "CDS" and "origin of transfer" in product.lower()):
            cgv_type = "oriT"
            legend_name = "oriT"
            label = "origin of transfer"
        elif ftype == "oriC":
            cgv_type = "oriC"
            legend_name = "oriC"
            label = "origin of replication"
        elif ftype in ("ncRNA", "ncRNA_region", "ncRNA-region", "tRNA", "rRNA", "tmRNA"):
            cgv_type = ftype.lower()
            legend_name = "CDS"
        else:
            cgv_type = "cds"
            legend_name = "CDS"

        cgv_features.append({
            "name": label,
            "type": cgv_type,
            "start": f["start"],
            "stop": f["end"],
            "strand": f["strand"],
            "source": bakta_source,
            "legend": legend_name,
            "contig": f.get("contig", list(ref_contigs.keys())[0]),
            "tags": [],
            "favorite": True,
        })

    # ── BLAST features ──
    # With the standard BLAST approach.
    import re

    def parse_orf_coords(orf_id):
        """Parse coordinates from ORF ID like 'contig;orf_1;start=100;end=500;strand=1;rf=1'
        Returns (contig_name, start, end, strand)"""
        m_start = re.search(r'start=(\d+)', orf_id)
        m_end = re.search(r'end=(\d+)', orf_id)
        m_strand = re.search(r'strand=(-?\d+)', orf_id)
        # Contig name is everything before the first ';'
        contig = orf_id.split(';')[0] if ';' in orf_id else None
        if m_start and m_end:
            return contig, int(m_start.group(1)), int(m_end.group(1)), int(m_strand.group(1)) if m_strand else 1
        return None, None, None, None

    for comp_idx, (comp_name, hits) in enumerate(blast_comps):
        source_name = f"blast-{comp_idx + 1}.1"

        # Keep best hit per QUERY ORF (reference ORF), by bitscore
        best = {}
        for h in hits:
            qid = h["qseqid"]
            if qid not in best or h["bitscore"] > best[qid]["bitscore"]:
                best[qid] = h

        for qid, hit in best.items():
            # Parse reference ORF coordinates from query ID
            ref_contig, ref_start, ref_end, ref_strand = parse_orf_coords(qid)
            if ref_start is None:
                # Fallback: try Bakta lookup (for skip-blast with old-format files)
                ref_orf = ref_lookup.get(qid)
                if not ref_orf:
                    for k, v in ref_lookup.items():
                        if qid in k or k in qid:
                            ref_orf = v
                            break
                if not ref_orf:
                    continue
                ref_start = ref_orf["start"]
                ref_end = ref_orf["end"]
                ref_strand = ref_orf["strand"]
                ref_contig = ref_orf.get("contig", list(ref_contigs.keys())[0])

            pident = hit["pident"]

            # Legend assignment: find the highest cutoff this hit meets
            lg = blast_legend_items[0]["name"]
            for bli in blast_legend_items:
                cutoff = bli["meta"]["blast_meta"]["blast_formatter"]["cutoff"]
                if cutoff == 100:
                    if pident >= 100:
                        lg = bli["name"]
                elif pident >= cutoff:
                    lg = bli["name"]

            feat = {
                "type": "blast",
                "start": ref_start,
                "stop": ref_end,
                "strand": ref_strand,
                "source": source_name,
                "legend": lg,
                "contig": ref_contig or list(ref_contigs.keys())[0],
                "tags": [],
                "meta": {
                    "subject": hit["sseqid"],
                    "percent_identity": round(pident, 3),
                    "mismatches": hit.get("mismatch", 0),
                    "alignment_length": hit.get("length", 0),
                    "gap_openings": hit.get("gapopen", 0),
                    "subject_start": hit.get("sstart", 1),
                    "subject_stop": hit.get("send", 0),
                    "evalue": hit.get("evalue", 0),
                    "bit_score": round(hit.get("bitscore", 0)),
                },
            }
            if pident < 100:
                feat["score"] = round(pident / 100, 3)

            cgv_features.append(feat)

    # ── Filter enclosed BLAST features ──
    # Remove BLAST features that are fully enclosed within larger features in the same source.
    ann_features = [f for f in cgv_features if f.get("type") != "blast"]
    blast_features = [f for f in cgv_features if f.get("type") == "blast"]

    filtered_blast = []
    for source in set(f["source"] for f in blast_features):
        src_feats = [f for f in blast_features if f["source"] == source]
        # Sort by span size descending (largest first)
        src_feats.sort(key=lambda f: -(f["stop"] - f["start"]))
        kept = []
        for feat in src_feats:
            enclosed = False
            feat_contig = feat.get("contig", "")
            for kf in kept:
                # Only compare features on the SAME contig
                if kf.get("contig", "") == feat_contig and kf["start"] <= feat["start"] and kf["stop"] >= feat["stop"]:
                    enclosed = True
                    break
            if not enclosed:
                kept.append(feat)
        filtered_blast.extend(kept)

    removed = len(blast_features) - len(filtered_blast)
    if removed > 0:
        print(f"  [FILTER] Removed {removed} enclosed BLAST features")
    cgv_features = ann_features + filtered_blast

    # ── Tracks ──
    # Layout: AMR (outermost, outside) → BLAST tracks (outside) → BACKBONE → Bakta (inside) → GC
    # Content → GC Skew (innermost)
    tracks = []

    # Calculate average identity for each comparison to sort by similarity
    comp_avg_identity = []
    for comp_idx, (comp_name, hits) in enumerate(blast_comps):
        if hits:
            avg_id = sum(h["pident"] for h in hits) / len(hits)
        else:
            avg_id = 0
        comp_avg_identity.append((comp_idx, comp_name, avg_id))

    # Sort: LOWEST identity first (outermost ring)
    comp_avg_identity.sort(key=lambda x: x[2])

    # Separate reference self-BLAST (index 0) from comparisons
    ref_entry = None
    comp_entries = []
    for entry in comp_avg_identity:
        if entry[0] == 0:
            ref_entry = entry
        else:
            comp_entries.append(entry)

    # ── Track layout ──
    # Order determines proximity to backbone (first = closest) around tracks straddle the
    # backbone; inside tracks go inside; outside goes outermost

    tracks.append({
        "name": "Bakta Annotation",
        "dataType": "feature",
        "dataMethod": "source",
        "dataKeys": bakta_source,
        "separateFeaturesBy": "strand",
        "position": "both",
        "thicknessRatio": 1,
        "favorite": True,
    })

    if ref_entry:
        comp_idx, comp_name, avg_id = ref_entry
        source_name = f"blast-{comp_idx + 1}.1"
        track_name = f"{comp_name} (reference)"
        tracks.append({
            "name": track_name,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": source_name,
            "separateFeaturesBy": "none",
            "position": "outside",
            "thicknessRatio": 1,
            "drawOrder": "score",
        })

    for comp_idx, comp_name, avg_id in comp_entries:
        source_name = f"blast-{comp_idx + 1}.1"
        track_name = comp_name
        if not track_name.endswith(".fasta"):
            track_name = f"{comp_name}.fasta"
        tracks.append({
            "name": track_name,
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": source_name,
            "separateFeaturesBy": "none",
            "position": "outside",
            "thicknessRatio": 1,
            "drawOrder": "score",
        })

    # ── AMR gene track: outermost outside ring ──
    amr_source = "amr-genes"
    amr_track_obj = None
    amr_gene_names = set()

    if amr_gene_data:
        # Build reverse mapping: original contig name -> cgview contig name
        orig_to_cgview = {}
        if fasta_to_gff:
            for orig, mapped in fasta_to_gff.items():
                orig_to_cgview[orig] = mapped
        try:
            with open(ref_fasta) as _rf_amr:
                for _line in _rf_amr:
                    if _line.startswith('>'):
                        _hdr = _line[1:].strip()
                        _cname = _hdr.split()[0]
                        _orig_match = re.search(r'original_name=(\S+)', _hdr)
                        if _orig_match:
                            orig_to_cgview[_orig_match.group(1)] = _cname
        except Exception:
            pass

        contig_keys = list(ref_contigs.keys())

        for gene_key, rf in amr_gene_data.items():
            gene = rf['gene']
            pident = rf['pident']
            coverage = rf['coverage']
            rf_contig = rf.get('contig', '')
            rf_start = rf.get('start', 0)
            rf_end = rf.get('end', 0)

            # Map AMRFinderPlus contig name to CGView contig name
            cgview_contig = orig_to_cgview.get(rf_contig, '')
            if not cgview_contig:
                if rf_contig in ref_contigs:
                    cgview_contig = rf_contig
                else:
                    # Contig not in reference (e.g., chromosomal AMR): skip
                    continue

            strand = 1 if rf_end >= rf_start else -1
            feat_start = min(rf_start, rf_end)
            feat_end = max(rf_start, rf_end)

            label = gene
            amr_gene_names.add(gene)

            cgv_features.append({
                "name": label,
                "type": "amr",
                "start": feat_start,
                "stop": feat_end,
                "strand": strand,
                "source": amr_source,
                "legend": gene,
                "contig": cgview_contig,
                "tags": {
                    "identity": f"{pident:.1f}%",
                    "coverage": f"{coverage:.1f}%",
                },
                "favorite": True,
            })

    if amr_gene_names:
        # Add legend entry per unique AMR gene
        amr_colors = [
            "rgba(220,60,60,1)", "rgba(235,87,87,1)", "rgba(200,40,40,1)",
            "rgba(239,108,77,1)", "rgba(178,34,34,1)", "rgba(205,92,92,1)",
            "rgba(233,72,72,1)", "rgba(210,50,50,1)", "rgba(195,65,65,1)",
            "rgba(224,78,78,1)",
        ]
        for i, gene in enumerate(sorted(amr_gene_names)):
            legend_items.append({
                "name": gene,
                "swatchColor": amr_colors[i % len(amr_colors)],
                "decoration": "arc",
            })

        # AMR track: outermost outside ring
        amr_track_obj = {
            "name": "AMR Genes",
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": amr_source,
            "separateFeaturesBy": "none",
            "position": "outside",
            "thicknessRatio": 1,
        }

    if amr_track_obj:
        tracks.append(amr_track_obj)

    tracks.append({
        "name": "GC Content",
        "separateFeaturesBy": "strand",
        "position": "inside",
        "thicknessRatio": 4,
        "dataType": "plot",
        "dataMethod": "sequence",
        "dataKeys": "gc-content",
        "dataOptions": {},
    })

    tracks.append({
        "name": "GC Skew",
        "separateFeaturesBy": "strand",
        "position": "inside",
        "thicknessRatio": 4,
        "dataType": "plot",
        "dataMethod": "sequence",
        "dataKeys": "gc-skew",
        "dataOptions": {},
    })

    # GC Content and GC Skew legend items
    legend_items.extend([
        {"name": "GC Content", "swatchColor": "rgba(0,0,0,1)", "decoration": "arc"},
        {"name": "GC Skew+", "swatchColor": "rgba(0,153,0,1)", "decoration": "arc"},
        {"name": "GC Skew-", "swatchColor": "rgba(153,0,153,1)", "decoration": "arc"},
    ])

    # ── Build sequence with actual sequence data (multi-contig) ──
    contigs_list = []
    for contig_name, contig_seq in ref_contigs.items():
        is_circ = circularity.get(contig_name, False)
        contigs_list.append({
            "name": contig_name,
            "orientation": "+",
            "length": len(contig_seq),
            "seq": contig_seq,
        })

    # ── Contig boundary track (draft multi-contig references) ──
    # When the reference is a fragmented MOB-suite bin rather than a closed replicon, the contigs
    # are drawn end-to-end.
    if len(contigs_list) > 1:
        boundary_source = "contig-boundaries"
        shades = [
            "rgba(90,110,130,0.85)", "rgba(140,120,100,0.85)",
            "rgba(110,140,120,0.85)", "rgba(140,105,125,0.85)",
        ]
        for idx, cinfo in enumerate(contigs_list, 1):
            label = f"contig {idx} ({cinfo['length']:,} bp)"
            legend_items.append({
                "name": label,
                "swatchColor": shades[(idx - 1) % len(shades)],
                "decoration": "arc",
            })
            cgv_features.append({
                "name": label,
                "type": "contig",
                "start": 1,
                "stop": cinfo["length"],
                "strand": 1,
                "source": boundary_source,
                "legend": label,
                "contig": cinfo["name"],
                "tags": [],
                "meta": {
                    "contig_name": cinfo["name"],
                    "length": cinfo["length"],
                    "circular": bool(circularity.get(cinfo["name"], False)),
                },
            })
        tracks.append({
            "name": "Reference Contigs",
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": boundary_source,
            "separateFeaturesBy": "none",
            "position": "inside",
            "thicknessRatio": 0.6,
        })

    # Update display name with circularity info
    circ_names = [k for k, v in circularity.items() if v]
    if circ_names:
        display_name += f" [circular: {', '.join(circ_names)}]"

    seq_obj = {
        "font": "sans-serif,plain,14",
        "color": "rgb(0,0,0)",
        "contigs": contigs_list,
    }

    # ── Assemble final JSON ──
    import datetime
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    cgview_json = {
        "cgview": {
            "version": "1.8.1",
            "created": now,
            "updated": now,
            "name": display_name,
            "settings": {
                "format": "circular",
                "geneticCode": 11,
                "backgroundColor": "rgba(255,255,255,1)",
                "showShading": True,
                "arrowHeadLength": 0.3,
                "initialMapThicknessProportion": 0.1,
                "maxMapThicknessProportion": 0.5,
            },
            "backbone": {
                "color": "rgba(128,128,128,1)",
                "colorAlternate": "rgba(200,200,200,1)",
                "thickness": 5,
                "decoration": "arc",
            },
            "ruler": {
                "font": "sans-serif,plain,10",
                "color": "rgba(100,100,100,1)",
                "minTickWidth": 9999,
            },
            "annotation": {
                "font": "monospace,plain,12",
                "onlyDrawFavorites": True,
                "labelPlacement": "default",
                "visible": True,
            },
            "dividers": {
                "slot": {
                    "visible": False,
                },
            },
            "highlighter": {"visible": False},
            "captions": [],
            "legend": {
                "position": "top-right",
                "textAlignment": "left",
                "defaultFont": "sans-serif,plain,13",
                "defaultFontColor": "rgba(0,0,0,1)",
                "defaultMinArcLength": 1,
                "backgroundColor": "rgba(255,255,255,0.7)",
                "items": legend_items,
            },
            "sequence": seq_obj,
            "features": cgv_features,
            "plots": [],
            "bookmarks": [],
            "tracks": tracks,
            "meta": {
                "blast_meta": {
                    "blast_formatter": {
                        "colorScheme": "custom",
                    },
                },
            },
        }
    }

    return cgview_json


# GENERATE INTERACTIVE HTML (CGView.js via CDN)

def _find_template(name="cgview_template.html"):
    """Locate a bundled HTML template."""
    here = os.path.dirname(os.path.realpath(__file__))
    candidates = [
        os.path.join(here, "templates", name),
        os.path.join(here, name),
    ]
    for cand in candidates:
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(
        f"HTML template not found: {name} (searched: {', '.join(candidates)})"
    )


def generate_html(json_path, html_path, display_name="Plasmid Map", amr_genes=None, caption=""):
    """Generate HTML from template with embedded JSON."""
    import os
    json_data = open(json_path).read()
    tmpl_path = _find_template("cgview_template.html")
    tmpl = open(tmpl_path).read()
    amr_str = ', '.join(sorted(amr_genes)) if amr_genes else ''
    # Caption is injected into a JS string literal: escape quotes/newlines.
    caption_str = (caption or '').replace('\\', '\\\\').replace('"', '\\"').replace('\n', ' ')
    html = (tmpl
        .replace("DISPLAY_NAME_PLACEHOLDER", display_name)
        .replace("JSON_DATA_PLACEHOLDER", json_data)
        .replace("AMR_GENES_PLACEHOLDER", amr_str)
        .replace("MAP_CAPTION_PLACEHOLDER", caption_str))
    with open(html_path, "w") as f:
        f.write(html)
    print("[HTML] Written to: " + html_path)
    return html_path




def main():
    parser = argparse.ArgumentParser(
        description="Plasmid CGView.js Pipeline: Bakta + BLAST + CGView.js HTML",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Workflow (standard):
  - Reference plasmid: Full Bakta annotation (gene names, AMR, oriT/oriC)
  - ORF extraction: 6-frame translation of all sequences (not just annotated CDSs)
  - BLASTp: ref ORFs (query) vs comparison ORFs (subject/db)
  - CGView.js: interactive HTML with CDN rendering

Examples:
  # Full pipeline:
  python pipeline.py --ref ref.fasta --comp comp1.fasta comp2.fasta \\
      --bakta-db /path/to/db -o results/

  # Custom identity cutoffs:
  python pipeline.py --ref ref.fasta --comp comp1.fasta \\
      --bakta-db /path/to/db --identity-cutoffs 0 90 95 99 100 -o results/

  # Skip Bakta (use pre-computed reference annotation):
  python pipeline.py --skip-bakta --ref-dir bakta_ref/ \\
      --ref-fasta ref.fasta --comp comp1.fasta comp2.fasta -o results/
        """)
    parser.add_argument("--ref", help="Reference plasmid FASTA")
    parser.add_argument("--comp", nargs="+", help="Comparison plasmid FASTAs")
    parser.add_argument("--bakta-db", help="Bakta database path")
    parser.add_argument("--bakta-path", help="Path to bakta executable (if not in PATH)")
    parser.add_argument("--blast-path", help="Path to BLAST+ bin directory (if not in PATH)")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-bakta", action="store_true",
                        help="Skip Bakta for reference; use pre-computed output via --ref-dir")
    parser.add_argument("--ref-dir", help="Pre-computed Bakta output dir for reference")
    parser.add_argument("--ref-fasta", help="Reference FASTA (with --skip-bakta)")
    parser.add_argument("--skip-blast", action="store_true",
                        help="Skip BLAST; use pre-computed BLAST files via --blast-files")
    parser.add_argument("--blast-files", nargs="+")
    parser.add_argument("--comp-names", nargs="+")
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--identity-cutoffs", nargs="+", type=float, default=[0, 85, 90, 95, 100],
                        help="BLAST identity cutoffs for legend color bins (default: 0 85 90 95 100). "
                             "Example: --identity-cutoffs 0 70 85 95 99 100")
    parser.add_argument("--color-palette", default="warm",
                        choices=["warm", "cool", "earth"],
                        help="BLAST identity color palette (default: warm)")
    parser.add_argument("--ref-display-name", default=None,
                        help="Override reference display name in the map title")
    parser.add_argument("--amrfinder-tsv", default=None,
                        help="AMRFinderPlus TSV for AMR gene identity/coverage metrics")
    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Annotate REFERENCE with Bakta ──
    print("=" * 60)
    print("STEP 1: Annotate reference plasmid (Bakta)")
    print("=" * 60)

    if args.skip_bakta:
        ref_bakta = args.ref_dir
        ref_fasta = args.ref_fasta or args.ref
    else:
        if not args.ref or not args.bakta_db:
            parser.error("--ref and --bakta-db required")
        ref_fasta = args.ref
        nm = Path(args.ref).stem
        ref_bakta = str(outdir / f"bakta_{nm}")
        run_bakta(args.ref, ref_bakta, args.bakta_db, prefix=nm,
                  threads=args.threads, bakta_path=args.bakta_path)

    gff, ref_faa, bakta_fna = find_bakta(ref_bakta)
    feats = parse_gff3(gff)
    ref_display = args.ref_display_name if args.ref_display_name else Path(ref_fasta).stem
    # Safe version for filenames (no spaces, parens, or special chars)
    ref_display_safe = re.sub(r'[^a-zA-Z0-9_\-.]', '_', ref_display)

    # ── Contig name resolution ──
    # GFF contig names are the single source of truth (e.g., contig_1, contig_2).

    orig_seqs, orig_circularity = parse_fasta(ref_fasta)
    fasta_contig_names = list(orig_seqs.keys())

    # Get unique GFF contig names in numeric order
    gff_contig_names = sorted(set(f["contig"] for f in feats),
                              key=lambda x: int(x.split('_')[-1]) if x.split('_')[-1].isdigit() else 0)

    fasta_to_gff = {}
    if set(gff_contig_names) != set(fasta_contig_names):
        for i, fasta_name in enumerate(fasta_contig_names):
            if i < len(gff_contig_names):
                fasta_to_gff[fasta_name] = gff_contig_names[i]
            else:
                fasta_to_gff[fasta_name] = fasta_name
        seqs = {fasta_to_gff.get(k, k): v for k, v in orig_seqs.items()}
        circularity = {fasta_to_gff.get(k, k): v for k, v in orig_circularity.items()}
        print(f"  Contig mapping: {fasta_to_gff}")
    else:
        seqs = dict(orig_seqs)
        circularity = dict(orig_circularity)
        fasta_to_gff = {}

    # Build concatenated sequence and total length across all contigs
    seq_len = sum(len(s) for s in seqs.values())
    cds_n = sum(1 for f in feats if f["type"] == "CDS")
    circ_contigs = [k for k, v in circularity.items() if v]
    circ_info = f", circular: {', '.join(circ_contigs)}" if circ_contigs else ""
    print(f"  {ref_display}: {seq_len:,} bp ({len(seqs)} contig(s)){circ_info}, {cds_n} CDS, {len(feats)} features total")

    # ── Step 2: ORF-vs-ORF BLAST Comparisons (standard) ──
    # Standard ORF-vs-ORF approach.
    print("\n" + "=" * 60)
    print("STEP 2: ORF-vs-ORF BLAST Comparisons")
    print("=" * 60)

    blast_comps = []
    if args.skip_blast:
        names = args.comp_names or [Path(f).stem for f in (args.blast_files or [])]
        for bf, cn in zip(args.blast_files or [], names):
            hits = parse_blast6(bf)
            blast_comps.append((cn, hits))
            print(f"  {cn}: {len(hits)} hits (pre-computed)")
    else:
        # Extract reference ORFs (6-frame): used as QUERY for all comparisons
        orf_dir = str(outdir / "orfs")
        os.makedirs(orf_dir, exist_ok=True)

        # Use a FASTA with GFF-matching contig names for ORF extraction
        if fasta_to_gff:
            renamed_fasta = os.path.join(orf_dir, f"{ref_display_safe}_renamed.fasta")
            with open(renamed_fasta, 'w') as fout:
                for orig_name, seq_str in orig_seqs.items():
                    new_name = fasta_to_gff.get(orig_name, orig_name)
                    fout.write(f">{new_name}\n{seq_str}\n")
            ref_orfs_faa = extract_orfs(renamed_fasta, orf_dir)
        else:
            ref_orfs_faa = extract_orfs(ref_fasta, orf_dir)

        # Self-BLAST: reference vs itself (innermost comparison track for visual confirmation)
        print(f"  Self-comparison (reference vs reference)...")
        self_bout = str(outdir / f"blast_{ref_display_safe}_vs_self.txt")
        run_blast(ref_orfs_faa, ref_orfs_faa, self_bout, evalue=args.evalue,
                  threads=args.threads, blast_path=args.blast_path)
        self_hits = parse_blast6(self_bout)
        blast_comps.append((ref_display, self_hits))
        print(f"  {ref_display} (self): {len(self_hits)} hits")

        # Comparison plasmids (no self-BLAST - reference is already shown as backbone + CDS tracks)
        comp_fastas = args.comp or []
        if not comp_fastas:
            print("  (no comparison plasmids provided)")
        for cf in comp_fastas:
            cn = Path(cf).stem
            # Extract ORFs from comparison (6-frame)
            comp_orfs_faa = extract_orfs(cf, orf_dir)
            # BLAST: reference ORFs (query) vs comparison ORFs (db)
            bout = str(outdir / f"blast_{ref_display_safe}_vs_{cn}.txt")
            run_blast(ref_orfs_faa, comp_orfs_faa, bout, evalue=args.evalue,
                      threads=args.threads, blast_path=args.blast_path)
            hits = parse_blast6(bout)
            blast_comps.append((cn, hits))
            print(f"  {cn}: {len(hits)} hits")

    if not blast_comps:
        print("  (no comparisons)")

    # ── Step 3: Build CGView JSON ──
    print("\n" + "=" * 60)
    print("STEP 3: Build CGView.js JSON (standard)")
    print("=" * 60)

    # Build display name: the --ref-display-name is already clean from select_comparison_ref.py
    display_name = f"{ref_display} - {seq_len:,} bp"
    # Parse AMRFinderPlus output for AMR gene identity/coverage
    amr_data = {}  # gene_name -> {pident, coverage, contig, start, end}
    if args.amrfinder_tsv and os.path.exists(args.amrfinder_tsv):
        import csv as _csv_amr
        with open(args.amrfinder_tsv) as _af:
            for row in _csv_amr.DictReader(_af, delimiter='\t'):
                # Acquired core genes only: POINT rows are mutations with no standalone feature to draw, and
                # Scope 'plus' (efflux, stress response) is excluded from the report throughout.
                if (row.get('Subtype') or '').strip().upper() != 'AMR':
                    continue
                if (row.get('Scope') or '').strip().lower() == 'plus':
                    continue
                gene = (row.get('Element symbol') or '').strip()
                if not gene:
                    continue
                pident = float(row.get('% Identity to reference') or 0)
                coverage = float(row.get('% Coverage of reference') or 0)
                contig = (row.get('Contig id') or '').strip()
                contig = contig.split()[0] if contig else ''
                start = int(float(row.get('Start') or 0))
                end = int(float(row.get('Stop') or 0))
                amr_data[gene.lower()] = {
                    'gene': gene, 'pident': pident, 'coverage': coverage,
                    'contig': contig, 'start': start, 'end': end,
                }
        print(f"  AMRFinderPlus: {len(amr_data)} AMR genes with identity/coverage")

    cgview = build_cgview_json(
        seq_len=seq_len,
        ref_contigs=seqs,
        features=feats,
        blast_comps=blast_comps,
        display_name=display_name,
        identity_cutoffs=args.identity_cutoffs,
        circularity=circularity,
        amr_gene_data=amr_data,
        fasta_to_gff=fasta_to_gff,
        ref_fasta=ref_fasta,
        color_palette=getattr(args, 'color_palette', 'warm'),
    )
    json_path = str(outdir / f"{ref_display_safe}_cgview.json")
    with open(json_path, "w") as f:
        json.dump(cgview, f, indent=2)
    n_feats = len(cgview["cgview"]["features"])
    n_tracks = len(cgview["cgview"]["tracks"])
    print(f"  JSON: {json_path}")
    print(f"  {n_feats} features, {n_tracks} tracks")

    # ── Step 4: Generate interactive HTML ──
    print("\n" + "=" * 60)
    print("STEP 4: Generate interactive HTML (CGView.js CDN)")
    print("=" * 60)

    html_path = str(outdir / f"{ref_display_safe}_map.html")
    # Collect AMR gene names from the CGView features actually placed on the map
    amr_gene_list = sorted(set(
        f['legend'] for f in cgview['cgview']['features']
        if f.get('source') == 'amr-genes' and f.get('legend')
    ))
    # Caption: only when the reference is a fragmented multi-contig bin, where contig adjacency on
    # the circle is an artifact of ordering, not evidence.
    caption = ""
    if len(seqs) > 1:
        caption = (
            f"Draft reference: {len(seqs)} contigs from a MOB-suite plasmid bin, drawn end-to-end "
            "in descending length/score order. Contig boundaries are marked on the 'Reference Contigs' "
            "track: adjacency between contigs is an ordering artifact and is NOT evidence of physical "
            "linkage. AMR genes absent from this reference cannot appear on the map; see "
            "amr_gene_summary.tsv for the full cluster gene set."
        )
    generate_html(json_path, html_path, display_name=display_name, amr_genes=amr_gene_list,
                  caption=caption)

    print(f"\n{'=' * 60}")
    print("DONE")
    print(f"{'=' * 60}")
    print(f"  JSON: {json_path}")
    print(f"  HTML: {html_path}")
    print(f"\n  Open the HTML file in a browser for interactive viewing.")
    print(f"  The JSON is compatible with CGView.js viewers")


if __name__ == "__main__":
    main()