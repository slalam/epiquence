/* Core genome analysis: Bakta, Panaroo, SNP distances, IQ-TREE */

process BAKTA_CGA {
    label 'process_high'
    tag "${cluster_name}_${sample_name}"
    cache 'lenient'
    publishDir "${params.outdir}/outbreak/bakta/${cluster_name}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(cluster_name), val(sample_name), path(fasta)

    output:
    tuple val(cluster_name), val(sample_name), path("${sample_name}.gff3"), emit: gff

    script:
    """
    export DIAMOND_BLOCK_SIZE=2
    export MPLCONFIGDIR=\$(mktemp -d)
    bakta --db ${params.bakta_db} \
        --prefix ${sample_name} \
        --output . \
        --force \
        --threads ${task.cpus} \
        --skip-plot \
        --compliant \
        ${fasta}
    """
}

process PANAROO {
    label 'process_high'
    tag "${cluster_name}"
    publishDir "${params.outdir}/outbreak/panaroo/${cluster_name}", mode: 'copy'

    input:
    tuple val(cluster_name), path(gff_files)

    output:
    tuple val(cluster_name), path("panaroo_out/core_gene_alignment_filtered.aln"), emit: core_aln

    script:
    """
    # Resolve .fna files from the same directory as each .gff3 symlink target.
    mkdir -p fna_dir
    for gff in ${gff_files}; do
        REAL_GFF=\$(readlink -f "\$gff")
        GFF_DIR=\$(dirname "\$REAL_GFF")
        PREFIX=\$(basename "\$gff" .gff3)
        FNA="\${GFF_DIR}/\${PREFIX}.fna"
        if [ -f "\$FNA" ]; then
            cp "\$FNA" fna_dir/
        fi
    done

    # Convert Bakta GFF3 + FNA to Prokka-compatible combined format
    bakta_to_prokka_gff.py ${gff_files} --fna-dir fna_dir --outdir converted_gffs

    panaroo -i converted_gffs/*.gff3 \
        -o panaroo_out \
        --clean-mode strict \
        -a core \
        --remove-invalid-genes \
        --core_threshold 0.99 \
        --threads ${task.cpus}
    """
}

process SNP_DISTS_CGA {
    label 'process_low'
    tag "${cluster_name}"
    publishDir "${params.outdir}/outbreak/cga_snp_dists/${cluster_name}", mode: 'copy'

    input:
    tuple val(cluster_name), path(core_aln)

    output:
    tuple val(cluster_name), path("${cluster_name}_cga_snp_matrix.tsv"), emit: snp_matrix

    script:
    """
    snp-dists ${core_aln} > ${cluster_name}_cga_snp_matrix.tsv
    """
}

process IQTREE_CGA {
    label 'process_high'
    tag "${cluster_name}"
    publishDir "${params.outdir}/outbreak/cga_iqtree/${cluster_name}", mode: 'copy'

    input:
    tuple val(cluster_name), path(core_aln)

    output:
    tuple val(cluster_name), path("${cluster_name}_cga_iqtree.treefile"), emit: tree
    tuple val(cluster_name), path("${cluster_name}_cga_iqtree.*"), emit: all_files

    script:
    """
    N_SEQS=\$(grep -c "^>" ${core_aln})
    THREADS=${task.cpus}
    if [ "\$THREADS" -gt 8 ]; then THREADS=8; fi

    echo "IQ-TREE CGA: \$N_SEQS sequences, \$THREADS threads"

    if [ "\$N_SEQS" -lt 3 ]; then
        echo "Skipping IQ-TREE: need at least 3 sequences, got \$N_SEQS"
        touch ${cluster_name}_cga_iqtree.treefile
    else
        if [ "\$N_SEQS" -ge 4 ]; then
            BOOT_FLAG="-bb 1000"
        else
            BOOT_FLAG=""
        fi
        iqtree3 -s ${core_aln} -st DNA -m GTR+G \$BOOT_FLAG -T \$THREADS -pre ${cluster_name}_cga_iqtree
    fi
    """
}

process RENDER_CGA_TREES {
    label 'process_low'
    tag "render_cga_trees"
    publishDir "${params.outdir}/outbreak/cga_trees", mode: 'copy'

    input:
    path(tree_files)
    path(metadata)
    path(outbreak_summary)
    path('mlst_files/*')
    path(chromosomal_amr)

    output:
    path "cga_tree_svgs", emit: svg_dir

    script:
    def meta_arg = metadata.name != 'NO_METADATA' ? "--metadata ${metadata}" : ""
    def tree_args = tree_files.name != 'NO_CGA_TREES' ? "--cluster-trees ${tree_files}" : ""
    def amr_arg = chromosomal_amr.name != 'NO_AMR' ? "--chromosomal-amr ${chromosomal_amr}" : ""
    """
    outbreak_render_trees.py \
        ${tree_args} \
        ${meta_arg} \
        --outbreak-summary ${outbreak_summary} \
        --mlst-dir mlst_files \
        ${amr_arg} \
        --outdir cga_tree_svgs \
        --prefix cga_
    """
}

process ASSEMBLE_CGA_RESULTS {
    label 'process_low'
    tag "assemble_cga_results"

    input:
    path(snp_files)
    path(tree_files)

    output:
    path "cga_cluster_results", emit: results_dir

    script:
    def has_trees = tree_files.name != 'NO_CGA_TREES'
    def has_snps = snp_files.name != 'NO_CGA_SNP'
    """
    mkdir -p cga_cluster_results

    if [ "${has_snps}" = "true" ]; then
        for f in ${snp_files}; do
            if [ -f "\$f" ]; then
                cluster=\$(basename "\$f" _cga_snp_matrix.tsv)
                mkdir -p "cga_cluster_results/\$cluster"
                cp "\$f" "cga_cluster_results/\$cluster/cga_snp_matrix.tsv"
            fi
        done
    fi

    if [ "${has_trees}" = "true" ]; then
        for f in ${tree_files}; do
            if [ -f "\$f" ]; then
                cluster=\$(basename "\$f" _cga_iqtree.treefile)
                mkdir -p "cga_cluster_results/\$cluster"
                cp "\$f" "cga_cluster_results/\$cluster/cga_iqtree.treefile"
            fi
        done
    fi

    echo "Assembled CGA cluster results:"
    find cga_cluster_results -type f | sort || echo "No CGA cluster results"
    """
}