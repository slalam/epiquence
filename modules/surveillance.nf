/*
 * Cumulative surveillance: SKA2 + ReporTree on chromosome bins
 * Pages per species in surveillance_results_dir: all/, <group_id>/, ungrouped/.
 */

include { runSettingsArg } from './utils'

// Root for the published surveillance views.
def surveillanceResultsRoot() {
    return params.surveillance_results_dir
        ? params.surveillance_results_dir
        : "${params.surveillance_sketch_dir}/surveillance-results"
}

// '--group <id>' for a group scope, nothing for the species-wide one.
def groupArg(group) {
    if (group == '__ALL__') return ''
    return "--group '" + group.toString().replace("'", "'\"'\"'") + "'"
}

process EXTRACT_CHROMOSOME {
    label 'process_low'
    tag "${sample_id}_${asm_type}"

    input:
    tuple val(sample_id), val(asm_type), path(mob_dir)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}.fasta"), emit: chromosome, optional: true

    script:
    def min_bp = params.surveillance_min_chromosome_bp ?: 500000
    """
    CHROMO="${mob_dir}/chromosome.fasta"

    # No chromosome bin at all: MOB-recon assigned everything to plasmid/unassigned.
    if [ ! -s "\$CHROMO" ]; then
        echo "WARN: no chromosome.fasta in ${mob_dir}" >&2
        echo "WARN: ${sample_id}_${asm_type} excluded from surveillance" >&2
        exit 0
    fi

    TOTAL_BP=\$(grep -v '^>' "\$CHROMO" | tr -d '\\n\\r ' | wc -c)

    # Implausibly small bin: usually a fragmented short-read assembly where chromosomal contigs
    # were misbinned as plasmid.
    if [ "\$TOTAL_BP" -lt ${min_bp} ]; then
        echo "WARN: chromosome bin for ${sample_id}_${asm_type} is \$TOTAL_BP bp (< ${min_bp} bp)" >&2
        echo "WARN: ${sample_id}_${asm_type} excluded from surveillance" >&2
        exit 0
    fi

    # Rename to {sample_id}_{asm_type} so the sketch carries the sample identity.
    cp "\$CHROMO" "${sample_id}_${asm_type}.fasta"
    echo "Chromosome for ${sample_id}_${asm_type}: \$TOTAL_BP bp"
    """
}

process SKA_BUILD {
    label 'process_low'
    tag "${species_name}_${sample_id}"
    publishDir "${params.surveillance_sketch_dir}/${species_safe}/sketches", mode: 'copy'

    input:
    tuple val(species_name), val(sample_id), path(fasta)

    output:
    tuple val(species_name), val(sample_id), path("${sample_id}.skf"), emit: sketch

    script:
    species_safe = species_name.replaceAll(/\s+/, '_')
    """
    # Archive any previous version of this isolate before the new sketch is published over it.
    SKDIR='${params.surveillance_sketch_dir}/${species_safe}/sketches'
    CHDIR='${params.surveillance_sketch_dir}/${species_safe}/chromosomes'
    mkdir -p "\$SKDIR/superseded" "\$CHDIR/superseded" "\$CHDIR"

    if [ -f "\$SKDIR/${sample_id}.skf" ]; then
        N=1
        while [ -e "\$SKDIR/superseded/${sample_id}_O\${N}.skf" ]; do N=\$((N+1)); done
        mv "\$SKDIR/${sample_id}.skf" "\$SKDIR/superseded/${sample_id}_O\${N}.skf"
        if [ -f "\$CHDIR/${sample_id}.fasta" ]; then
            mv "\$CHDIR/${sample_id}.fasta" "\$CHDIR/superseded/${sample_id}_O\${N}.fasta"
        fi
        echo "${sample_id}: previous version archived as ${sample_id}_O\${N}" >&2
    fi

    cp -f ${fasta} "\$CHDIR/${sample_id}.fasta"

    ska build -o ${sample_id} -k 31 --threads ${task.cpus} ${fasta}
    """
}

process SKA_MERGE {
    label 'process_medium'
    tag "${species_name}"
    publishDir "${params.surveillance_sketch_dir}/${species_safe}/merged", mode: 'copy'

    input:
    tuple val(species_name), path(sketches), path(prev_merged)

    output:
    tuple val(species_name), path("${species_safe}_merged.skf"), emit: merged

    script:
    species_safe = species_name.replaceAll(/\s+/, '_')
    def prev_arg = prev_merged.name != 'NO_PREV_MERGED' ? "${prev_merged}" : ""
    """
    # Rebuild from the per-isolate sketches rather than merging into the previous merged sketch.
    mkdir -p sketch_inputs

    # Historical sketches first...
    HIST='${params.surveillance_sketch_dir}/${species_safe}/sketches'
    if [ -d "\$HIST" ]; then
        cp \$HIST/*.skf sketch_inputs/ 2>/dev/null || true
    fi

    # ...then this run's on top, so a re-sequenced isolate replaces its older sketch by filename.
    for skf in ${sketches}; do
        if [ -f "\$skf" ]; then
            cp -f "\$skf" sketch_inputs/
        fi
    done

    # Superseded versions are archived out of sketches/ and are excluded by default.
    if [ "${params.surveillance_include_superseded}" = "true" ]; then
        SUP='${params.surveillance_sketch_dir}/${species_safe}/chromosomes/superseded'
        if [ -d "\$SUP" ]; then
            for fa in \$SUP/*.fasta; do
                [ -e "\$fa" ] || continue
                name=\$(basename "\$fa" .fasta)
                ska build -o "sketch_inputs/\$name" -k 31 --threads ${task.cpus} "\$fa"
                echo "Included superseded version \$name" >&2
            done
        fi
    fi

    echo "Merging \$(ls sketch_inputs/*.skf 2>/dev/null | wc -l) sketch(es) for ${species_name}" >&2

    N_SKETCHES=\$(ls sketch_inputs/*.skf 2>/dev/null | wc -l)

    if [ "\$N_SKETCHES" -ge 2 ]; then
        ska merge -o ${species_safe}_merged sketch_inputs/*.skf
    else
        # Single isolate: just copy the sketch as the "merged" file
        cp sketch_inputs/*.skf ${species_safe}_merged.skf
    fi
    """
}

process SKA_DISTANCE {
    label 'process_medium'
    tag "${species_name}"

    input:
    tuple val(species_name), path(merged_sketch)

    output:
    tuple val(species_name), path("${species_safe}_distances.tsv"), emit: distances

    script:
    species_safe = species_name.replaceAll(/\s+/, '_')
    """
    # Check how many samples are in the merged sketch
    ska nk ${merged_sketch} >&2 || true

    # ska distance computes pairwise distances between all samples within a merged .skf
    if ska distance ${merged_sketch} > ${species_safe}_distances.tsv; then
        echo "SKA distance completed. Lines: \$(wc -l < ${species_safe}_distances.tsv)" >&2
    else
        echo "SKA distance failed or single sample; creating empty distances" >&2
        echo -e "Sample1\\tSample2\\tDistance\\tMismatches (proportion)\\tMatch count\\tMismatch count" > ${species_safe}_distances.tsv
    fi
    """
}

process REPORTREE {
    label 'process_medium'
    tag "${species_name}"
    publishDir "${params.surveillance_sketch_dir}/${species_safe}/merged/reportree", mode: 'copy'

    input:
    tuple val(species_name), path(distances), path(prev_partitions)

    output:
    tuple val(species_name), path("reportree_out/partitions.tsv"), emit: partitions
    tuple val(species_name), path("reportree_out/nomenclature_changes.tsv"), emit: nomenclature_changes, optional: true
    tuple val(species_name), path("reportree_out"), emit: output_dir

    script:
    species_safe = species_name.replaceAll(/\s+/, '_')
    def prev_arg = prev_partitions.name != 'NO_PREV_PARTITIONS' ? "--nomenclature-file ${prev_partitions}" : ""
    """
    mkdir -p reportree_out

    # Count non-empty data lines (skip header, ignore blank lines)
    N_PAIRS=\$(tail -n +2 ${distances} | grep -c '[^[:space:]]' || true)
    echo "Distance file has \$N_PAIRS data rows" >&2

    if [ "\$N_PAIRS" -gt 0 ]; then
        # Convert SKA2 pairwise list to symmetric distance matrix for ReporTree
        python3 << 'PYEOF'
import csv, sys

pairs = {}
samples = set()
with open('${distances}') as f:
    reader = csv.DictReader(f, delimiter='\\t')
    print(f"Columns: {reader.fieldnames}", file=sys.stderr)
    for row in reader:
        rl = {k.strip().lower().replace(' ', ''): v for k, v in row.items()}
        s1 = rl.get('sample1', rl.get('sample 1', ''))
        s2 = rl.get('sample2', rl.get('sample 2', ''))
        dist = rl.get('distance', rl.get('snpdistance', rl.get('snps', '0')))
        if s1 and s2:
            try:
                d = int(float(dist))
            except (ValueError, TypeError):
                d = 0
            samples.add(s1)
            samples.add(s2)
            pairs[(s1, s2)] = d
            pairs[(s2, s1)] = d

samples = sorted(samples)
print(f"Parsed {len(samples)} unique samples from {len(pairs)//2} pairs", file=sys.stderr)

if len(samples) >= 2:
    with open('distance_matrix.tsv', 'w') as out:
        out.write('\\t' + '\\t'.join(samples) + '\\n')
        for s1 in samples:
            row = [s1]
            for s2 in samples:
                if s1 == s2:
                    row.append('0')
                else:
                    row.append(str(pairs.get((s1, s2), 0)))
            out.write('\\t'.join(row) + '\\n')

    # Create minimal metadata file (required for --nomenclature-code-levels)
    with open('metadata.tsv', 'w') as mf:
        mf.write('sequence\\n')
        for s in samples:
            mf.write(s + '\\n')
else:
    open('distance_matrix.tsv', 'w').close()
    print("No valid sample pairs found", file=sys.stderr)
PYEOF

        # Only run ReporTree if matrix has content
        if [ -s distance_matrix.tsv ]; then
            reportree.py \
                -d_mx distance_matrix.tsv \
                -m metadata.tsv \
                -o ReporTree \
                --HC-threshold single-10,single-20,single-50,single-100 \
                --nomenclature-code-levels 100,50,20,10 \
                ${prev_arg}

            echo "All output files:" >&2
            ls -la ReporTree* *.tsv *.nwk 2>/dev/null >&2

            # Collect ReporTree output into reportree_out/
            for f in ReporTree_partitions.tsv ReporTree_clusterComposition.tsv ReporTree_nomenclature_changes.tsv ReporTree_single_HC.nwk ReporTree.log ReporTree_metadata_w_partitions.tsv; do
                if [ -f "\$f" ]; then cp "\$f" reportree_out/; fi
            done

            # Copy partitions to expected name
            if [ -f reportree_out/ReporTree_partitions.tsv ]; then
                cp reportree_out/ReporTree_partitions.tsv reportree_out/partitions.tsv
            fi
            # Copy nomenclature changes if present
            if [ -f reportree_out/ReporTree_nomenclature_changes.tsv ]; then
                cp reportree_out/ReporTree_nomenclature_changes.tsv reportree_out/nomenclature_changes.tsv
            fi
        else
            echo "Empty distance matrix; creating minimal partitions" >&2
            echo -e "sequence\\tsingle-100x1.0\\tsingle-50x1.0\\tsingle-20x1.0\\tsingle-10x1.0" > reportree_out/partitions.tsv
        fi
    else
        # Single isolate or no pairs: create minimal partitions file
        echo -e "sequence\\tsingle-100x1.0\\tsingle-50x1.0\\tsingle-20x1.0\\tsingle-10x1.0" > reportree_out/partitions.tsv
    fi

    # Ensure output files exist
    if [ ! -f reportree_out/partitions.tsv ]; then
        echo "ERROR: partitions.tsv not generated" >&2
        echo "Files in reportree_out:" >&2
        ls -la reportree_out/ >&2
        FOUND=\$(find reportree_out -name "*partitions*" -print -quit)
        if [ -n "\$FOUND" ]; then
            echo "Found partitions file: \$FOUND, renaming" >&2
            cp "\$FOUND" reportree_out/partitions.tsv
        else
            exit 1
        fi
    fi

    # Log partition column count
    echo "Partitions columns:" >&2
    head -1 reportree_out/partitions.tsv | tr '\\t' '\\n' | wc -l >&2
    """
}

process PREPARE_ISOLATES_TSV {
    label 'process_low'
    tag "${species_name}"

    input:
    tuple val(species_name), val(sample_ids)
    path(general_stats)

    output:
    tuple val(species_name), path("${species_safe}_isolates.tsv"), emit: isolates_tsv

    script:
    species_safe = species_name.replaceAll(/\s+/, '_')
    def sample_list = sample_ids.join(',')
    """
    python3 -c "
import csv, sys

samples = '${sample_list}'.split(',')
gs_data = {}
if '${general_stats}' != 'NO_FILE':
    with open('${general_stats}') as f:
        for row in csv.DictReader(f, delimiter='\\t'):
            gs_data[row.get('sample', '')] = row

with open('${species_safe}_isolates.tsv', 'w') as out:
    w = csv.writer(out, delimiter='\\t')
    w.writerow(['isolate_id', 'sample', 'assembly_type', 'mlst_scheme', 'mlst_st', 'sketch_path', 'fasta_path'])
    for s in samples:
        gs = gs_data.get(s, {})
        # s is sample with suffix (e.g., SAMPLE-01_L).
        parts = s.rsplit('_', 1)
        if len(parts) == 2 and parts[1] in ('L','S','HL','HS'):
            iso_id = parts[0]
            asm = parts[1]
        else:
            iso_id = s
            asm = ''
        w.writerow([iso_id, s, asm, gs.get('mlst_scheme',''), gs.get('mlst_st',''), '', ''])
"
    """
}

process LOAD_SURVEILLANCE_DB {
    label 'process_low'
    tag "${species_name}"
    cache false
    maxForks 1

    input:
    tuple val(species_name), path(isolates_tsv), path(distances), path(partitions), path(nomenclature_changes), val(merged_sketch_path)

    output:
    tuple val(species_name), path("surveillance_load.log"), emit: log
    tuple val(species_name), path("partitions_for_next_run.tsv"), emit: saved_partitions, optional: true

    script:
    species_safe = species_name.replaceAll(/\s+/, '_')
    def run_id = "${params.workflow_id}_${species_safe}_${new Date().format('yyyyMMddHHmmss')}"
    def changes_arg = nomenclature_changes.name != 'NO_CHANGES' ? "--nomenclature-changes ${nomenclature_changes}" : ""
    def distances_arg = distances.name != 'NO_DISTANCES' ? "--distances-tsv ${distances}" : ""
    def partitions_arg = partitions.name != 'NO_PARTITIONS' ? "--partitions-tsv ${partitions}" : ""
    def partitions_dir = "${params.surveillance_sketch_dir}/${species_safe}/merged"
    def db_arg = runSettingsArg(workflow.launchDir)
    """
    load_surveillance_db.py \
        ${db_arg} \
        --species "${species_name}" \
        --workflow-id "${params.workflow_id ?: ''}" \
        --workflow-description "${params.workflow_description ?: ''}" \
        --group-id "${params.group_id ?: ''}" \
        --group-description "${params.group_description ?: ''}" \
        --run-id "${run_id}" \
        --isolates-tsv ${isolates_tsv} \
        ${distances_arg} \
        ${partitions_arg} \
        ${changes_arg} \
        --merged-sketch-path "${merged_sketch_path}" \
        --partitions-output partitions_for_next_run.tsv \
        2>&1 | tee surveillance_load.log
    """
}
process UPDATE_CUMULATIVE_ISOLATES {
    label 'process_low'
    tag "${species_name}"
    // Published under a species-qualified name so the next run's file().exists() lookup finds it.
    publishDir "${params.surveillance_sketch_dir}/${species_safe}", mode: 'copy', overwrite: true,
               saveAs: { "${species_safe}_isolates_cumulative.tsv" }

    input:
    tuple val(species_name), path(current_isolates), path(prev_cumulative)
    path(chromosomal_amr)

    output:
    tuple val(species_name), path("cumulative_isolates.tsv"), emit: cumulative

    script:
    species_safe = species_name.replaceAll(/\s+/, '_')
    def prev_arg = prev_cumulative.name != 'NO_PREV_ISOLATES' ? "--prev ${prev_cumulative}" : ''
    def amr_arg  = chromosomal_amr.name  != 'NO_AMR'          ? "--chromosomal-amr ${chromosomal_amr}" : ''
    // Shell-quoted: descriptions are free text and may contain apostrophes.
    def grp_desc = "'" + (params.group_description ?: '').toString().replace("'", "'\"'\"'") + "'"
    """
    update_cumulative_isolates.py \\
        --current ${current_isolates} \\
        ${prev_arg} \\
        ${amr_arg} \\
        --species '${species_name}' \\
        --workflow-id '${params.workflow_id ?: ''}' \\
        --group-id '${params.group_id ?: ''}' \\
        --group-description ${grp_desc} \\
        --timezone '${params.report_timezone ?: ''}' \\
        --output cumulative_isolates.tsv
    """
}

process SURVEILLANCE_TREE {
    label 'process_low'
    tag "${species_name}:${scope_dir}"
    // One HTML per species per scope, replaced in place each run.
    publishDir "${surveillanceResultsRoot()}/${scope_dir}/species-trees", mode: 'copy', overwrite: true

    input:
    tuple val(species_name), val(scope_dir), val(group), path(distances), path(cumulative_isolates), path(partitions)

    output:
    tuple val(species_name), path("${species_safe}_surveillance_tree.html"), emit: tree, optional: true

    script:
    species_safe = species_name.replaceAll(/\s+/, '_')
    def part_arg = partitions.name != 'NO_PARTITIONS' ? "--partitions ${partitions}" : ''
    def grp_arg  = groupArg(group)
    """
    surveillance_tree.py \\
        --distances ${distances} \\
        --isolates ${cumulative_isolates} \\
        ${part_arg} \\
        ${grp_arg} \\
        --species '${species_name}' \\
        --snp-cap ${params.surveillance_tree_snp_cap} \\
        --max-pairs ${params.surveillance_tree_max_pairs} \\
        --timezone '${params.report_timezone ?: ''}' \\
        --output ${species_safe}_surveillance_tree.html
    """
}

process SURVEILLANCE_DISTANCE_MATRIX {
    label 'process_low'
    tag "${species_name}:${scope_dir}"
    // One page (+ TSV) per species per scope, replaced in place each run.
    publishDir "${surveillanceResultsRoot()}/${scope_dir}/distance-matrix", mode: 'copy', overwrite: true

    input:
    tuple val(species_name), val(scope_dir), val(group), path(distances), path(cumulative_isolates)

    output:
    tuple val(species_name), path("${species_safe}_distance_matrix.html"), emit: matrix, optional: true
    // Square matrix covering every isolate, in the same clustered order as the rendered one.
    tuple val(species_name), path("${species_safe}_distance_matrix.tsv"), emit: matrix_tsv, optional: true

    script:
    species_safe = species_name.replaceAll(/\s+/, '_')
    def grp_arg  = groupArg(group)
    """
    surveillance_distance_matrix.py \\
        --distances ${distances} \\
        --isolates ${cumulative_isolates} \\
        ${grp_arg} \\
        --species '${species_name}' \\
        --color-break ${params.surveillance_distance_color_break} \\
        --vector-max ${params.surveillance_distance_vector_max} \\
        --max-isolates ${params.surveillance_distance_max_isolates} \\
        --timezone '${params.report_timezone ?: ''}' \\
        --matrix-out ${species_safe}_distance_matrix.tsv \\
        --output ${species_safe}_distance_matrix.html
    """
}