/* Plasmid HGT detection and comparison maps */

process COLLECT_PLASMID_FASTAS {
    label 'process_low'
    tag "${asm_type}"

    input:
    val(asm_type)
    path(mob_dirs)
    path(amrfinder_files)

    output:
    tuple val(asm_type), path("plasmids_${asm_type}"), emit: plasmid_dir

    script:
    // Short-read assemblies keep the MOB bin as the unit: one plasmid is genuinely spread over
    // many contigs there, so splitting would shatter a single plasmid into many points.
    def contig_level = params.plasmid_hgt_contig_types.split(',')*.trim().contains(asm_type) ? '1' : '0'
    """
    mkdir -p plasmids_${asm_type}
    CONTIG_LEVEL=${contig_level}
    echo "Plasmid unit for ${asm_type}: \$([ "\$CONTIG_LEVEL" = "1" ] && echo contig || echo 'MOB bin')" >&2

    for d in ${mob_dirs}; do
        # Only process dirs matching this assembly type
        if echo "\$d" | grep -qE "_${asm_type}_mob_recon\$"; then
            sample=\$(echo "\$d" | sed 's/_${asm_type}_mob_recon\$//')
            # MOB-recon outputs plasmid FASTAs as plasmid_*.fasta
            for pf in \$d/plasmid_*.fasta; do
                if [ -f "\$pf" ]; then
                    base=\$(basename "\$pf" .fasta)
                    PREFIX="\${sample}_${asm_type}_\${base}"
                    if [ "\$CONTIG_LEVEL" = "1" ]; then
                        # One FASTA per contig.
                        split_plasmid_contigs.py \\
                            --input "\$pf" \\
                            --prefix "\$PREFIX" \\
                            --outdir "plasmids_${asm_type}" \\
                            --min-length ${params.plasmid_min_length}
                    else
                        cp "\$pf" "plasmids_${asm_type}/\${PREFIX}.fasta"
                    fi
                fi
            done
        fi
    done

    # Co-locate the AMRFinderPlus hits so the reference-selection and map scripts can resolve them
    # by sample prefix.
    if [ "${amrfinder_files}" != "NO_AMRFINDER_FILES" ]; then
        for f in ${amrfinder_files}; do
            if echo "\$f" | grep -qE "_${asm_type}_amrfinder\\.tsv\$"; then
                cp "\$f" "plasmids_${asm_type}/\$(basename \$f)"
            fi
        done
    fi

    echo "Collected plasmids for ${asm_type}:"
    ls -la plasmids_${asm_type}/ | head -30
    echo "Total: \$(ls plasmids_${asm_type}/*.fasta 2>/dev/null | wc -l) plasmid sequence FASTAs"
    echo "AMRFinder files: \$(ls plasmids_${asm_type}/*_amrfinder.tsv 2>/dev/null | wc -l)"
    """
}

process MASH_DIST_PLASMIDS {
    label 'process_medium'
    tag "${asm_type}"
    publishDir "${params.outdir}/plasmid_hgt/${asm_type}/mash", mode: 'copy'

    input:
    tuple val(asm_type), path(plasmid_dir)

    output:
    tuple val(asm_type), path("mash_dist_matrix.tsv"), emit: mash_matrix
    tuple val(asm_type), path(plasmid_dir), emit: plasmid_dir_out

    script:
    """
    N_PLASMIDS=\$(ls ${plasmid_dir}/*.fasta 2>/dev/null | wc -l)
    echo "Found \$N_PLASMIDS plasmid FASTAs"

    if [ "\$N_PLASMIDS" -lt 2 ]; then
        echo "Too few plasmids for clustering"
        echo -e "plasmid" > mash_dist_matrix.tsv
        exit 0
    fi

    mash sketch -s ${params.hgt_sketch_size} -o all_plasmids ${plasmid_dir}/*.fasta
    mash dist all_plasmids.msh all_plasmids.msh -t > raw_mash_dist.tsv
    sed 's|[^\t]*/||g; s|\\.fasta||g' raw_mash_dist.tsv > mash_dist_matrix.tsv
    """
}

process CLUSTER_PLASMIDS {
    label 'process_low'
    tag "${asm_type}"
    publishDir "${params.outdir}/plasmid_hgt/${asm_type}/mash", mode: 'copy'

    input:
    tuple val(asm_type), path(mash_matrix)

    output:
    tuple val(asm_type), path("plasmid_clusters.tsv"), emit: clusters

    script:
    """
    plasmid_hgt_cluster.py \
        --mash-matrix ${mash_matrix} \
        --threshold ${params.hgt_mash_threshold} \
        --output plasmid_clusters.tsv
    """
}

process UMAP_PLASMIDS {
    label 'process_low'
    tag "${asm_type}"
    cache false
    publishDir "${params.outdir}/plasmid_hgt/${asm_type}/umap", mode: 'copy'

    input:
    tuple val(asm_type), path(mash_matrix), path(clusters), path(metadata)

    output:
    tuple val(asm_type), path("${asm_type}_umap_embedding.tsv"), emit: embedding, optional: true
    tuple val(asm_type), path("${asm_type}_umap_embed.html"), emit: umap_html, optional: true
    tuple val(asm_type), path("${asm_type}_umap_plot.html"), emit: umap_standalone, optional: true

    script:
    def meta_arg = metadata.name != 'NO_HGT_METADATA' ? "--metadata ${metadata}" : ""
    """
    plasmid_umap.py \
        --mash-matrix ${mash_matrix} \
        --clusters ${clusters} \
        --asm-type ${asm_type} \
        ${meta_arg} \
        --output-prefix ${asm_type}

    # Create empty files if UMAP skipped (too few plasmids)
    touch -a ${asm_type}_umap_embedding.tsv
    touch -a ${asm_type}_umap_embed.html
    touch -a ${asm_type}_umap_plot.html
    """
}

process SELECT_AMR_PLASMIDS {
    label 'process_low'
    tag "${asm_type}"
    // Python step is kept separate from mashtree: the staphb/mashtree image has no python3.
    container 'docker://quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(asm_type), path(plasmid_dir)

    output:
    tuple val(asm_type), path("amr_plasmids"), path("fallback.nwk"), emit: selected
    path("${asm_type}_amr_plasmids.tsv"), emit: summary, optional: true

    script:
    def min_genomes = params.amr_tree_min_genomes ?: 3
    """
    # fallback.nwk stays empty unless too few plasmids qualify for a real tree.
    touch fallback.nwk
    mkdir -p amr_plasmids

    set +e
    select_amr_plasmids.py \
        --plasmid-dir ${plasmid_dir} \
        --outdir amr_plasmids \
        --summary ${asm_type}_amr_plasmids.tsv \
        --min-genomes ${min_genomes} \
        --fallback-newick fallback.nwk
    STATUS=\$?
    set -e

    # 0 = enough for mashtree, 3 = fallback written; anything else is a real failure
    if [ "\$STATUS" -ne 0 ] && [ "\$STATUS" -ne 3 ]; then
        echo "ERROR: select_amr_plasmids.py failed with status \$STATUS" >&2
        exit \$STATUS
    fi
    """
}

process MASHTREE_AMR_PLASMIDS {
    label 'process_low'
    tag "${asm_type}"
    container 'docker://staphb/mashtree:1.4.6'
    publishDir "${params.outdir}/plasmid_hgt/${asm_type}/trees", mode: 'copy'

    input:
    tuple val(asm_type), path(amr_dir), path(fallback_nwk)

    output:
    tuple val(asm_type), path("${asm_type}_amr_plasmids.nwk"), emit: tree

    script:
    """
    if [ -s ${fallback_nwk} ]; then
        # Too few AMR plasmids for a real tree; use the placeholder verbatim.
        cp ${fallback_nwk} ${asm_type}_amr_plasmids.nwk
    else
        mashtree --numcpus ${task.cpus} ${amr_dir}/*.fasta > ${asm_type}_amr_plasmids.nwk
    fi

    if [ ! -s ${asm_type}_amr_plasmids.nwk ]; then
        echo "ERROR: no Newick produced for ${asm_type}" >&2
        exit 1
    fi
    """
}

process AMR_PLASMID_MATRIX {
    label 'process_low'
    tag "${asm_type}"
    cache false
    // Same numpy/scipy image the outbreak SNP heatmap uses: no new dependency.
    container 'docker://quay.io/biocontainers/mulled-v2-0fd299cadb7a80e2cc704b5d903ccc54893c512d:377fcafe3b6e7ef703094f8f47a64d081622ee09-2'
    publishDir "${params.outdir}/plasmid_hgt/${asm_type}/mash", mode: 'copy'

    input:
    tuple val(asm_type), path(mash_matrix), path(amr_dir), path(metadata)

    output:
    tuple val(asm_type), path("${asm_type}_amr_mash_matrix.html"), emit: matrix_html
    tuple val(asm_type), path("${asm_type}_amr_mash_matrix.tsv"),  emit: matrix_tsv
    tuple val(asm_type), path("${asm_type}_amr_mash_order.txt"),   emit: leaf_order

    script:
    def meta_arg = metadata.name != 'NO_HGT_METADATA' ? "--metadata ${metadata}" : ""
    def max_arg  = params.hgt_matrix_max_plasmids ?: 120
    """
    plasmid_amr_matrix.py \\
        --mash-matrix ${mash_matrix} \\
        --amr-dir ${amr_dir} \\
        --asm-type ${asm_type} \\
        ${meta_arg} \\
        --threshold ${params.hgt_mash_threshold} \\
        --max-plasmids ${max_arg} \\
        --output-prefix ${asm_type}

    # Keep the channel shape stable when the matrix was skipped (< 2 AMR plasmids).
    touch -a ${asm_type}_amr_mash_matrix.html \\
             ${asm_type}_amr_mash_matrix.tsv \\
             ${asm_type}_amr_mash_order.txt
    """
}


process SELECT_PLASMID_REF {
    label 'process_low'
    tag "${asm_type}_${cluster_id}"
    publishDir "${params.outdir}/plasmid_hgt/${asm_type}/comparison_maps/${cluster_id}", mode: 'copy',
               pattern: "{amr_presence_matrix.tsv,amr_gene_summary.tsv,ref_info.txt,groups_manifest.json}"

    input:
    tuple val(asm_type), val(cluster_id), path(plasmid_fastas), path(amrfinder_files)

    output:
    tuple val(asm_type), val(cluster_id), path("ref_output"), path(plasmid_fastas), path(amrfinder_files), emit: selected
    path("amr_presence_matrix.tsv"), optional: true, emit: amr_matrix
    path("amr_gene_summary.tsv"),    optional: true, emit: amr_summary
    path("ref_info.txt"),            optional: true, emit: ref_info
    path("groups_manifest.json"),    optional: true, emit: manifest

    script:
    // Fall back if the deployment's nextflow.config predates this param, otherwise.
    def group_mode = params.plasmid_map_group_mode ?: 'auto'
    """
    PRIORITY_ARG=""
    if [ -n "${params.amr_priority_genes}" ]; then
        PRIORITY_ARG="--priority-genes ${params.amr_priority_genes}"
    fi

    select_comparison_ref.py \
        --fastas *.fasta \
        --output-dir ref_output \
        --min-length ${params.plasmid_min_length} \
        --circular-min ${params.plasmid_circular_min} \
        --fragment-min ${params.plasmid_fragment_min} \
        --containment-cov ${params.plasmid_containment_cov} \
        --containment-id ${params.plasmid_containment_id} \
        --group-mode ${group_mode} \
        --asm-type ${asm_type} \
        \$PRIORITY_ARG

    # Surface companion/provenance files as task-level outputs for publishing
    for f in amr_presence_matrix.tsv amr_gene_summary.tsv ref_info.txt groups_manifest.json; do
        [ -f "ref_output/\$f" ] && cp "ref_output/\$f" "\$f"
    done
    true
    """
}

process BAKTA_PLASMID_REF {
    label 'process_high'
    tag "${asm_type}_${cluster_id}_${group_suffix}"
    maxForks 1

    input:
    tuple val(asm_type), val(cluster_id), val(group_suffix), path(ref_fasta), val(display_name), val(amrfinder_tsv), path(comp_fastas), path(amrfinder_files)

    output:
    tuple val(asm_type), val(cluster_id), val(group_suffix), path(ref_fasta), val(display_name), val(amrfinder_tsv), path("bakta_ref"), path(comp_fastas), path(amrfinder_files), emit: annotated

    script:
    """
    export MPLCONFIGDIR=\$(mktemp -d)
    mkdir -p bakta_ref

    bakta \
        --db ${params.bakta_db} \
        --output bakta_ref \
        --threads ${task.cpus} \
        --force \
        --prefix plasmid_ref \
        --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-crispr --skip-plot \
        ${ref_fasta}
    """
}

process COMPARISON_MAP {
    label 'process_high'
    tag "${asm_type}_${cluster_id}${group_suffix}"
    publishDir "${params.outdir}/plasmid_hgt/${asm_type}/comparison_maps/${cluster_id}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2
    maxForks 1

    input:
    tuple val(asm_type), val(cluster_id), val(group_suffix), path(ref_fasta), val(display_name), val(amrfinder_tsv), path(bakta_dir), path(comp_fastas), path(amrfinder_files)

    output:
    tuple val(asm_type), val(cluster_id), path("${cluster_id}${group_suffix}_map.html"), emit: map_html
    tuple val(asm_type), path("${cluster_id}${group_suffix}_map.html"), emit: map_html_flat

    script:
    def ref_sample_prefix = display_name.split('\\s')[0].split('\\(')[0].trim()
    def palette = params.comparison_map_palette ?: 'cool'
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    # Build AMRFinderPlus argument if available
    AMRFINDER_ARG=""
    if [ -n "${amrfinder_tsv}" ] && [ -f "${amrfinder_tsv}" ]; then
        AMRFINDER_ARG="--amrfinder-tsv ${amrfinder_tsv}"
        echo "AMRFinderPlus: ${amrfinder_tsv}"
    fi

    # Build comp list: exclude ref FASTA, ref sample's original FASTA, and non-FASTA files
    COMPS=""
    REF_SAMPLE_BASE="${ref_sample_prefix}"
    for f in *.fasta; do
        # Skip the group ref FASTA
        [ "\$f" = "${ref_fasta}" ] && continue
        # Skip other group ref FASTAs (ref_a.fasta, ref_b.fasta, etc.)
        echo "\$f" | grep -qE '^ref(_[a-z])?\\.fasta\$' && continue
        # Skip the reference sample's original plasmid FASTA
        STEM=\$(basename "\$f" .fasta)
        [ "\$STEM" = "\$REF_SAMPLE_BASE" ] && continue
        COMPS="\$COMPS \$f"
    done

    COMPS=\$(echo "\$COMPS" | sed 's/^ *//')

    echo "Reference: ${display_name}"
    echo "Comparisons: \$COMPS"

    if [ -z "\$COMPS" ]; then
        echo "No comparison samples; generating self-only map"
    fi

    # Run with pre-computed Bakta
    generate_comparison_map.py \
        --skip-bakta \
        --ref ${ref_fasta} \
        --ref-dir ${bakta_dir} \
        --comp \$COMPS \
        --threads ${task.cpus} \
        --identity-cutoffs 0 85 90 95 100 \
        --ref-display-name "${display_name}" \
        --color-palette ${palette} \
        \$AMRFINDER_ARG \
        --output results/

    cp results/*_map.html ${cluster_id}${group_suffix}_map.html
    """
}

process PLASMID_ANI_AF {
    label 'process_medium'
    tag "${asm_type}_${cluster_id}"
    publishDir "${params.outdir}/plasmid_hgt/${asm_type}/ani/${cluster_id}", mode: 'copy'

    input:
    tuple val(asm_type), val(cluster_id), path(plasmid_fastas)

    output:
    tuple val(asm_type), val(cluster_id), path("${cluster_id}_ani_af.tsv"), emit: ani_af

    script:
    """
    mkdir -p results

    FASTAS=(\$(ls *.fasta))
    N=\${#FASTAS[@]}

    echo -e "query\\treference\\tANI\\tAF_query\\tAF_ref\\tTotalSNPs\\tTotalIndels\\tAlignedBases\\tquery_size\\tref_size" > ${cluster_id}_ani_af.tsv

    # Pairwise ANI/AF using nucmer --maxmatch + dnadiff
    for i in \$(seq 0 \$((\$N - 1))); do
        for j in \$(seq \$((\$i + 1)) \$((\$N - 1))); do
            QRY=\${FASTAS[\$i]}
            REF_F=\${FASTAS[\$j]}
            QRY_NAME=\$(basename \$QRY .fasta)
            REF_F_NAME=\$(basename \$REF_F .fasta)

            PREFIX="results/\${QRY_NAME}_vs_\${REF_F_NAME}"
            nucmer --maxmatch -p \$PREFIX \$REF_F \$QRY 2>/dev/null || true
            if [ -f "\${PREFIX}.delta" ]; then
                dnadiff -d \${PREFIX}.delta -p \$PREFIX 2>/dev/null || true
            fi

            if [ -f "\${PREFIX}.report" ]; then
                ANI=\$(grep "AvgIdentity" \${PREFIX}.report | head -1 | awk '{print \$2}')
                AF_R_RAW=\$(grep "AlignedBases" \${PREFIX}.report | head -1 | awk '{print \$2}' | sed 's/(.*//')
                AF_Q_RAW=\$(grep "AlignedBases" \${PREFIX}.report | head -1 | awk '{print \$3}' | sed 's/(.*//')
                SNPS=\$(grep "TotalSNPs" \${PREFIX}.report | head -1 | awk '{print \$2}')
                INDELS=\$(grep "TotalIndels" \${PREFIX}.report | head -1 | awk '{print \$2}')
                ALIGNED=\$(grep "AlignedBases" \${PREFIX}.report | head -1 | awk '{print \$2}')

                # If ANI is 0 or empty but we have aligned bases and 0 SNPs, sequences are identical
                if [ "\${ANI}" = "0.00" ] || [ "\${ANI}" = "0" ] || [ -z "\${ANI}" ]; then
                    if [ "\${SNPS:-0}" = "0" ] && [ "\${AF_R_RAW:-0}" != "0" ]; then
                        ANI="100.00"
                    fi
                fi

                R_TOTAL=\$(grep "TotalBases" \${PREFIX}.report | head -1 | awk '{print \$2}')
                Q_TOTAL=\$(grep "TotalBases" \${PREFIX}.report | head -1 | awk '{print \$3}')
                AF_Q_PCT=\$(awk -v aligned="\$AF_Q_RAW" -v total="\$Q_TOTAL" 'BEGIN{if(total>0) printf "%.2f", aligned/total*100; else print 0}')
                AF_R_PCT=\$(awk -v aligned="\$AF_R_RAW" -v total="\$R_TOTAL" 'BEGIN{if(total>0) printf "%.2f", aligned/total*100; else print 0}')

                echo -e "\${QRY_NAME}\\t\${REF_F_NAME}\\t\${ANI}\\t\${AF_Q_PCT}\\t\${AF_R_PCT}\\t\${SNPS:-0}\\t\${INDELS:-0}\\t\${ALIGNED}\\t\${Q_TOTAL}\\t\${R_TOTAL}" >> ${cluster_id}_ani_af.tsv
                echo -e "\${REF_F_NAME}\\t\${QRY_NAME}\\t\${ANI}\\t\${AF_R_PCT}\\t\${AF_Q_PCT}\\t\${SNPS:-0}\\t\${INDELS:-0}\\t\${ALIGNED}\\t\${R_TOTAL}\\t\${Q_TOTAL}" >> ${cluster_id}_ani_af.tsv
            else
                echo -e "\${QRY_NAME}\\t\${REF_F_NAME}\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0" >> ${cluster_id}_ani_af.tsv
                echo -e "\${REF_F_NAME}\\t\${QRY_NAME}\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0" >> ${cluster_id}_ani_af.tsv
            fi
        done
    done
    """
}

process RENDER_HGT_TREES {
    label 'process_low'
    tag "${asm_type}"
    cache false
    publishDir "${params.outdir}/plasmid_hgt/${asm_type}/trees", mode: 'copy'

    input:
    tuple val(asm_type), path(tree_files), path(metadata)
    path('gtdbtk_files/*')

    output:
    tuple val(asm_type), path("hgt_tree_svgs_${asm_type}"), emit: svg_dir

    script:
    def meta_arg = metadata.name != 'NO_HGT_METADATA' ? "--metadata ${metadata}" : ""
    """
    plasmid_hgt_render_trees.py \
        --tree-files ${tree_files} \
        ${meta_arg} \
        --gtdbtk-dir gtdbtk_files \
        --outdir hgt_tree_svgs_${asm_type}
    """
}

process BUILD_HGT_REPORT {
    label 'process_low'
    tag "${asm_type}"
    cache false
    publishDir "${params.outdir}/plasmid_hgt/${asm_type}", mode: 'copy'

    input:
    tuple val(asm_type), path(cluster_results_dir), path(metadata), path(tree_svgs), path(map_htmls), path(umap_html), path(amr_matrix_html)

    output:
    path "plasmid_hgt_${asm_type}_table_mqc.yaml", emit: table_yaml
    path "plasmid_hgt_${asm_type}_plots_mqc.yaml", emit: plots_yaml
    path "plasmid_hgt_${asm_type}_summary_mqc.yaml", emit: summary_yaml
    path "plasmid_hgt_${asm_type}_combined.tsv", emit: combined_tsv, optional: true

    script:
    """
    # Stage map HTMLs into a directory for the report script
    mkdir -p comparison_maps
    for f in ${map_htmls}; do
        if [ -f "\$f" ]; then
            cp "\$f" comparison_maps/
        fi
    done

    plasmid_hgt_report.py \
        --asm-type ${asm_type} \
        --results-dir ${cluster_results_dir} \
        --metadata ${metadata} \
        --tree-svgs ${tree_svgs} \
        --comparison-maps comparison_maps \
        --umap-html ${umap_html} \
        --amr-matrix-html ${amr_matrix_html} \
        --table-yaml plasmid_hgt_${asm_type}_table_mqc.yaml \
        --plots-yaml plasmid_hgt_${asm_type}_plots_mqc.yaml \
        --combined-tsv plasmid_hgt_${asm_type}_combined.tsv \
        --workflow-id ${params.workflow_id ?: 'NA'} \
        --group-id ${params.group_id ?: 'NA'}
    """
}