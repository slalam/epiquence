#!/usr/bin/env python3
"""Parse MOB-recon output into MultiQC-compatible TSV."""
import argparse, csv, os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--mob-dir', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    num_plasmid_contigs = num_chromosomes = 0
    plasmid_sizes = []
    rep_types = set()
    relaxase_types_list = []
    mash_nearest = set()
    primary_cluster_ids = set()
    secondary_cluster_ids = set()

    contig_report = f"{args.mob_dir}/contig_report.txt"
    if os.path.exists(contig_report):
        with open(contig_report, 'r') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                mol_type = row.get('molecule_type', '')
                if mol_type == 'plasmid':
                    num_plasmid_contigs += 1
                    try: plasmid_sizes.append(int(row.get('size', row.get('contig_length', '0'))))
                    except ValueError: pass
                    rep = row.get('rep_type', row.get('rep_type(s)', '-'))
                    if rep and rep != '-':
                        rep_types.update(r.strip() for r in rep.split(','))
                    pc = row.get('primary_cluster_id', '-')
                    if pc and pc != '-': primary_cluster_ids.add(pc)
                elif mol_type == 'chromosome':
                    num_chromosomes += 1

    mobtyper_results = f"{args.mob_dir}/mobtyper_results.txt"
    if os.path.exists(mobtyper_results):
        with open(mobtyper_results, 'r') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                rep = row.get('rep_type(s)', row.get('rep_type', '-'))
                if rep and rep != '-': rep_types.update(r.strip() for r in rep.split(','))
                relaxase = row.get('relaxase_type(s)', row.get('relaxase_type', '-'))
                if relaxase and relaxase != '-': relaxase_types_list.extend(r.strip() for r in relaxase.split(','))
                mash = row.get('mash_nearest_neighbor', row.get('mash_neighbor_identification', '-'))
                if mash and mash != '-': mash_nearest.add(mash)
                pc = row.get('primary_cluster_id', '-')
                if pc and pc != '-': primary_cluster_ids.add(pc)
                sc = row.get('secondary_cluster_id', '-')
                if sc and sc != '-': secondary_cluster_ids.add(sc)

    join_or_dash = lambda s: '; '.join(sorted(s)) if s else '-'
    total_plasmid_size = sum(plasmid_sizes)
    avg_plasmid_size = int(total_plasmid_size / len(plasmid_sizes)) if plasmid_sizes else 0

    with open(args.output, 'w') as out:
        out.write("Sample\tNum_Plasmids\tNum_Plasmid_Contigs\tNum_Chromosomes\tTotal_Plasmid_Size_bp\tAvg_Plasmid_Size_bp\tReplicon_Types\tRelaxase_Types\tMash_Nearest\tPrimary_Cluster_IDs\tSecondary_Cluster_IDs\n")
        out.write(f"{args.sample}\t{len(primary_cluster_ids)}\t{num_plasmid_contigs}\t{num_chromosomes}\t{total_plasmid_size}\t{avg_plasmid_size}\t{join_or_dash(rep_types)}\t{'; '.join(relaxase_types_list) if relaxase_types_list else '-'}\t{join_or_dash(mash_nearest)}\t{join_or_dash(primary_cluster_ids)}\t{join_or_dash(secondary_cluster_ids)}\n")

if __name__ == '__main__':
    main()