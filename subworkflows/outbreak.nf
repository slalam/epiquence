/* Outbreak detection */

include { COLLECT_CHROMOSOMES; MASHTREE; OUTBREAK_CLUSTER; SKA_ALIGNMENT; SNP_DISTS_CORE; IQTREE; RENDER_TREES } from '../modules/outbreak'

workflow OUTBREAK_DETECTION {
    take:
    mob_recon_dirs
    gtdbtk_parsed
    metadata_file
    mlst_parsed
    chromosomal_amr

    main:
    mob_dirs_collected = mob_recon_dirs
        .map { id, asm_type, dir -> dir }
        .collect()

    COLLECT_CHROMOSOMES(mob_dirs_collected)
    MASHTREE(COLLECT_CHROMOSOMES.out.chromo_dir)

    gtdbtk_collected = gtdbtk_parsed.collect()

    OUTBREAK_CLUSTER(
        MASHTREE.out.mash_matrix,
        COLLECT_CHROMOSOMES.out.chromo_dir,
        gtdbtk_collected,
        metadata_file
    )

    // Parse cluster list into channels
    all_clusters = OUTBREAK_CLUSTER.out.cluster_list
        .splitCsv(header: true, strip: true)

    // ── All clusters: SKA2 -> SNP-dists + IQ-TREE ──
    ska_ch = all_clusters
        .filter { it.run_ska == 'True' }
        .map { row -> row.cluster_name }
        .combine(OUTBREAK_CLUSTER.out.clusters_dir)
        .map { name, cdir -> [name, cdir] }

    SKA_ALIGNMENT(ska_ch)
    SNP_DISTS_CORE(SKA_ALIGNMENT.out.core_aln)
    IQTREE(SKA_ALIGNMENT.out.core_aln)

    all_snp = SNP_DISTS_CORE.out.snp_matrix
        .map { name, f -> f }
        .collect()
        .ifEmpty(file('NO_SNP_MATRICES'))

    tree_collected = IQTREE.out.tree
        .map { name, f -> f }
        .collect()
        .ifEmpty(file('NO_TREE_FILES'))

    ASSEMBLE_OUTBREAK_RESULTS(all_snp, tree_collected)

    iqtree_files = IQTREE.out.tree
        .map { name, f -> f }
        .collect()
        .ifEmpty(file('NO_TREE_FILES'))

    mlst_collected = mlst_parsed
        .collect()
        .ifEmpty(file('NO_MLST'))

    RENDER_TREES(
        MASHTREE.out.tree,
        iqtree_files,
        metadata_file,
        OUTBREAK_CLUSTER.out.summary,
        mlst_collected,
        chromosomal_amr,
        gtdbtk_collected
    )

    emit:
    overview_tree   = MASHTREE.out.tree
    summary         = OUTBREAK_CLUSTER.out.summary
    cluster_list    = OUTBREAK_CLUSTER.out.cluster_list
    clusters_dir    = OUTBREAK_CLUSTER.out.clusters_dir
    tree_svgs       = RENDER_TREES.out.svg_dir
    ska_results     = ASSEMBLE_OUTBREAK_RESULTS.out.results_dir
    mlst_collected  = mlst_collected
}

process ASSEMBLE_OUTBREAK_RESULTS {
    label 'process_low'
    tag "assemble_results"

    input:
    path(snp_files)
    path(tree_files)

    output:
    path "cluster_results", emit: results_dir

    script:
    def has_trees = tree_files.name != 'NO_TREE_FILES'
    def has_snps = snp_files.name != 'NO_SNP_MATRICES'
    """
    mkdir -p cluster_results

    if [ "${has_snps}" = "true" ]; then
        for f in ${snp_files}; do
            if [ -f "\$f" ]; then
                cluster=\$(basename "\$f" _snp_matrix.tsv)
                mkdir -p "cluster_results/\$cluster"
                cp "\$f" "cluster_results/\$cluster/snp_matrix.tsv"
            fi
        done
    fi

    if [ "${has_trees}" = "true" ]; then
        for f in ${tree_files}; do
            if [ -f "\$f" ]; then
                cluster=\$(basename "\$f" _iqtree.treefile)
                mkdir -p "cluster_results/\$cluster"
                cp "\$f" "cluster_results/\$cluster/iqtree.treefile"
            fi
        done
    fi

    echo "Assembled cluster results:"
    find cluster_results -type f | sort || echo "No cluster results"
    """
}