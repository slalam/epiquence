/* Core genome analysis */

include { BAKTA_CGA; PANAROO; SNP_DISTS_CGA; IQTREE_CGA; RENDER_CGA_TREES; ASSEMBLE_CGA_RESULTS } from '../modules/core_genome'

workflow CORE_GENOME_ANALYSIS {
    take:
    cluster_list        // from OUTBREAK_CLUSTER.out.cluster_list
    clusters_dir        // from OUTBREAK_CLUSTER.out.clusters_dir
    metadata_file
    outbreak_summary    // from OUTBREAK_CLUSTER.out.summary
    mlst_parsed
    chromosomal_amr

    main:
    // Parse cluster list and filter for >=5 samples
    all_clusters = cluster_list
        .splitCsv(header: true, strip: true)

    cga_ch = all_clusters
        .filter { it.run_cga == 'True' }
        .map { row -> row.cluster_name }
        .combine(clusters_dir)
        .map { name, cdir -> [name, cdir] }

    // Expand to per-sample Bakta inputs
    bakta_in = cga_ch
        .flatMap { cluster_name, clusters_dir ->
            def cluster_path = new File("${clusters_dir}/${cluster_name}")
            def fastas = cluster_path.listFiles().findAll { it.name.endsWith('.fasta') }
            fastas.collect { f ->
                def sample = f.name.replace('.fasta', '')
                [cluster_name, sample, f.toPath()]
            }
        }

    BAKTA_CGA(bakta_in)

    panaroo_in = BAKTA_CGA.out.gff
        .map { cluster_name, sample_name, gff -> [cluster_name, gff] }
        .groupTuple()

    PANAROO(panaroo_in)
    SNP_DISTS_CGA(PANAROO.out.core_aln)
    IQTREE_CGA(PANAROO.out.core_aln)

    cga_snp_collected = SNP_DISTS_CGA.out.snp_matrix
        .map { name, f -> f }
        .collect()
        .ifEmpty(file('NO_CGA_SNP'))

    cga_tree_collected = IQTREE_CGA.out.tree
        .map { name, f -> f }
        .collect()
        .ifEmpty(file('NO_CGA_TREES'))

    ASSEMBLE_CGA_RESULTS(cga_snp_collected, cga_tree_collected)

    mlst_collected = mlst_parsed
        .collect()
        .ifEmpty(file('NO_MLST'))

    RENDER_CGA_TREES(
        cga_tree_collected,
        metadata_file,
        outbreak_summary,
        mlst_collected,
        chromosomal_amr
    )

    emit:
    cga_tree_svgs   = RENDER_CGA_TREES.out.svg_dir
    cga_results     = ASSEMBLE_CGA_RESULTS.out.results_dir
}