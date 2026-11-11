/* Outbreak detection: Mash clustering, SKA2 SNPs, IQ-TREE, report */

process COLLECT_CHROMOSOMES {
    label 'process_low'
    tag "collect_chromosomes"

    input:
    path(mob_dirs)

    output:
    path "chromosomes", emit: chromo_dir

    script:
    """
    mkdir -p chromosomes
    for d in ${mob_dirs}; do
        if [ -d "\$d" ]; then
            sample=\$(basename "\$d" | sed 's/_mob_recon\$//')
            chromo="\$d/chromosome.fasta"
            if [ -f "\$chromo" ]; then
                cp "\$chromo" "chromosomes/\${sample}.fasta"
            fi
        fi
    done
    count=\$(ls chromosomes/*.fasta 2>/dev/null | wc -l)
    echo "Collected \$count chromosome assemblies"
    """
}

process MASHTREE {
    label 'process_medium'
    tag "mashtree"
    publishDir "${params.outdir}/outbreak/overview", mode: 'copy'

    input:
    path(chromo_dir)

    output:
    path "overview_tree.nwk", emit: tree, optional: true
    path "mash_distance_matrix.tsv", emit: mash_matrix

    script:
    """
    NUM_SAMPLES=\$(ls ${chromo_dir}/*.fasta 2>/dev/null | wc -l)

    if [ "\$NUM_SAMPLES" -ge 2 ]; then
        mashtree --numcpus ${task.cpus} ${chromo_dir}/*.fasta > overview_tree.nwk
    fi

    # Mash distance matrix (works with 1+ samples)
    mash sketch -o all_samples -s 10000 ${chromo_dir}/*.fasta
    mash dist all_samples.msh all_samples.msh -t > raw_mash_matrix.tsv

    # Clean sample names: remove path prefix and .fasta suffix
    sed 's|[^\t]*/||g; s|\\.fasta||g' raw_mash_matrix.tsv > mash_distance_matrix.tsv
    """
}

process OUTBREAK_CLUSTER {
    label 'process_low'
    tag "clustering"
    publishDir "${params.outdir}/outbreak", mode: 'copy'

    input:
    path(mash_matrix)
    path(chromo_dir)
    path(gtdbtk_files)
    path(metadata)

    output:
    path "outbreak_summary.json", emit: summary
    path "cluster_list.csv", emit: cluster_list
    path "clusters", emit: clusters_dir, optional: true

    script:
    def meta_arg = metadata.name != 'NO_METADATA' ? "--metadata ${metadata}" : ""
    """
    mkdir -p gtdbtk_dir
    for f in ${gtdbtk_files}; do
        if [ -f "\$f" ]; then cp "\$f" gtdbtk_dir/; fi
    done

    outbreak_cluster.py \
        --gtdbtk-dir gtdbtk_dir \
        --assemblies-dir ${chromo_dir} \
        --mash-matrix ${mash_matrix} \
        ${meta_arg} \
        --threshold ${params.outbreak_threshold} \
        --outdir .
    """
}

process SKA_ALIGNMENT {
    label 'process_high'
    tag "${cluster_name}"
    publishDir "${params.outdir}/outbreak/ska/${cluster_name}", mode: 'copy'

    input:
    tuple val(cluster_name), path(clusters_dir)

    output:
    tuple val(cluster_name), path("${cluster_name}_ska.aln"), emit: core_aln

    script:
    """
    mkdir -p input_seqs
    cp ${clusters_dir}/${cluster_name}/*.fasta input_seqs/

    echo "SKA2 input files for ${cluster_name}:"
    ls -la input_seqs/

    # Build ska sketch files for each sample
    for f in input_seqs/*.fasta; do
        SAMPLE=\$(basename "\$f" .fasta)
        ska build -o "\$SAMPLE" -k 31 --threads ${task.cpus} "\$f"
    done

    ska merge -o merged *.skf

    # Generate alignment (ska outputs to prefix without extension)
    ska align -o ska_out merged.skf --min-freq 0.8 --filter no-filter
    mv ska_out ${cluster_name}_ska.aln
    """
}


process SNP_DISTS_CORE {
    label 'process_low'
    tag "${cluster_name}"
    publishDir "${params.outdir}/outbreak/snp_dists/${cluster_name}", mode: 'copy'

    input:
    tuple val(cluster_name), path(core_aln)

    output:
    tuple val(cluster_name), path("${cluster_name}_snp_matrix.tsv"), emit: snp_matrix

    script:
    """
    snp-dists ${core_aln} > ${cluster_name}_snp_matrix.tsv
    # Clean .fasta from sample names in the matrix
    sed -i 's/\\.fasta//g' ${cluster_name}_snp_matrix.tsv
    """
}

process IQTREE {
    label 'process_high'
    tag "${cluster_name}"
    publishDir "${params.outdir}/outbreak/iqtree/${cluster_name}", mode: 'copy'

    input:
    tuple val(cluster_name), path(core_aln)

    output:
    tuple val(cluster_name), path("${cluster_name}_iqtree.treefile"), emit: tree
    tuple val(cluster_name), path("${cluster_name}_iqtree.*"), emit: all_files

    script:
    """
    N_SEQS=\$(grep -c "^>" ${core_aln})
    THREADS=${task.cpus}
    if [ "\$THREADS" -gt 8 ]; then THREADS=8; fi

    echo "IQ-TREE: \$N_SEQS sequences, \$THREADS threads"

    if [ "\$N_SEQS" -lt 3 ]; then
        echo "Skipping IQ-TREE: need at least 3 sequences, got \$N_SEQS"
        touch ${cluster_name}_iqtree.treefile
    else
        # Bootstrap requires 4+ sequences
        if [ "\$N_SEQS" -ge 4 ]; then
            BOOT_FLAG="-bb 1000"
        else
            BOOT_FLAG=""
        fi
        iqtree3 -s ${core_aln} -st DNA -m GTR+G \$BOOT_FLAG -T \$THREADS -pre ${cluster_name}_iqtree
    fi
    """
}



process RENDER_TREES {
    label 'process_low'
    tag "render_trees"
    publishDir "${params.outdir}/outbreak/trees", mode: 'copy'

    input:
    path(overview_tree)
    path(tree_files)
    path(metadata)
    path(outbreak_summary)
    path('mlst_files/*')
    path(chromosomal_amr)
    path('gtdbtk_files/*')

    output:
    path "tree_svgs", emit: svg_dir

    script:
    def meta_arg = metadata.name != 'NO_METADATA' ? "--metadata ${metadata}" : ""
    def tree_args = tree_files.name != 'NO_TREE_FILES' ? "--cluster-trees ${tree_files}" : ""
    def amr_arg = chromosomal_amr.name != 'NO_AMR' ? "--chromosomal-amr ${chromosomal_amr}" : ""
    """
    outbreak_render_trees.py \
        --overview-tree ${overview_tree} \
        ${tree_args} \
        ${meta_arg} \
        --outbreak-summary ${outbreak_summary} \
        --mlst-dir mlst_files \
        --gtdbtk-dir gtdbtk_files \
        ${amr_arg} \
        --outdir tree_svgs
    """
}

process BUILD_OUTBREAK_REPORT {
    label 'process_low'
    tag "build_report"
    cache false
    publishDir "${params.outdir}/outbreak", mode: 'copy'

    input:
    path(tree_svgs)
    path(cluster_results)
    path(outbreak_summary)
    path('ref_files/*')
    path('mlst_files/*')
    path(chromosomal_amr)
    path(metadata)
    path(cga_tree_svgs)
    path(cga_cluster_results)

    output:
    path "outbreak_section_mqc.yaml", emit: report_yaml
    path "outbreak_cluster_summary_mqc.yaml", emit: summary_table_yaml
    path "ska2_pairwise_distances.tsv", emit: pairwise_tsv
    path "cga_pairwise_distances.tsv", emit: cga_pairwise_tsv, optional: true

    script:
    def amr_arg = chromosomal_amr.name != 'NO_AMR' ? "--chromosomal-amr ${chromosomal_amr}" : ""
    def meta_arg = metadata.name != 'NO_METADATA' ? "--metadata ${metadata}" : ""
    def cga_tree_arg = cga_tree_svgs.name != 'NO_CGA_TREES' ? "--cga-tree-svgs ${cga_tree_svgs}" : ""
    def cga_results_arg = cga_cluster_results.name != 'NO_CGA_RESULTS' ? "--cga-cluster-results ${cga_cluster_results}" : ""
    """
    outbreak_build_report.py \
        --tree-svgs ${tree_svgs} \
        --cluster-results ${cluster_results} \
        --outbreak-summary ${outbreak_summary} \
        --ref-dir ref_files \
        --mlst-dir mlst_files \
        ${amr_arg} \
        ${meta_arg} \
        ${cga_tree_arg} \
        ${cga_results_arg} \
        --workflow-id "${params.workflow_id ?: ''}" \
        --workflow-timestamp "\$(date '+%Y-%m-%d %H:%M:%S')" \
        --group-id "${params.group_id ?: ''}" \
        --output outbreak_section_mqc.yaml
    """
}