process MLST {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_mlst.tsv"), emit: mlst_result

    script:
    """
    mlst ${assembly} > ${sample_id}_${asm_type}_mlst.tsv
    """
}

process MLST_WITH_SCHEME {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(asm_type), path(assembly), val(scheme)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_mlst.tsv"), emit: mlst_result

    script:
    """
    mlst --scheme ${scheme} ${assembly} > ${sample_id}_${asm_type}_mlst.tsv
    """
}

process PARSE_MLST {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(mlst_file)

    output:
    path "${sample_id}_${asm_type}_mlst_parsed.tsv", emit: for_multiqc

    script:
    def sample_name = "${sample_id}_${asm_type}"
    """
    parse_mlst.py --sample ${sample_name} --input ${mlst_file} --output ${sample_name}_mlst_parsed.tsv
    """
}

process GTDBTK {
    label 'process_max'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_gtdbtk"), emit: gtdbtk_dir
    path "${sample_id}_${asm_type}_gtdbtk_summary.tsv", emit: raw_summary, optional: true

    script:
    def outdir = "${sample_id}_${asm_type}_gtdbtk"
    def sample_name = "${sample_id}_${asm_type}"
    """
    mkdir -p genome_dir
    cp ${assembly} genome_dir/${sample_name}.fasta

    export GTDBTK_DATA_PATH=${params.gtdbtk_db}
    export TQDM_DISABLE=1

    gtdbtk classify_wf \
        --cpus ${task.cpus} \
        --pplacer_cpus ${task.cpus} \
        --genome_dir genome_dir \
        --extension fasta \
        --out_dir ${outdir}

    for f in ${outdir}/*.summary.tsv; do
        if [ -f "\$f" ]; then
            cp "\$f" "${sample_name}_gtdbtk_summary.tsv"
            break
        fi
    done
    """
}

process PARSE_GTDBTK {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(gtdbtk_dir)

    output:
    path "${sample_id}_${asm_type}_gtdbtk_parsed.tsv", emit: parsed_taxonomy

    script:
    def sample_name = "${sample_id}_${asm_type}"
    """
    parse_gtdbtk.py --sample ${sample_name} --gtdbtk-dir ${gtdbtk_dir} --output ${sample_name}_gtdbtk_parsed.tsv
    """
}

// sylph: alternative species engine (taxonomy_tool 'sylph' or 'both'); the assembly is the sample

// Species engine for this run.
def taxonomyTool(p) {
    def v = p.containsKey('taxonomy_tool') ? p.get('taxonomy_tool') : null
    return v ? v.toString().trim().toLowerCase() : 'gtdbtk'
}

process SYLPH_QUERY {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_sylph_query.tsv"), emit: hits

    script:
    """
    sylph query ${params.sylph_db} -r ${assembly} -t ${task.cpus} \\
        > ${sample_id}_${asm_type}_sylph_query.tsv
    """
}

process PARSE_SYLPH {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(query_tsv)
    path(taxonomy_files, stageAs: 'gtdb_taxonomy/*')
    val(out_suffix)

    output:
    path "${sample_id}_${asm_type}_${out_suffix}.tsv", emit: parsed_taxonomy

    script:
    // out_suffix is 'gtdbtk_parsed' when sylph is the species engine, so every downstream reader
    // of that filename works unchanged; 'sylph_parsed' in 'both' mode, where GTDB-Tk stays
    // authoritative and sylph is the comparison.
    def sample_name = "${sample_id}_${asm_type}"
    """
    parse_sylph.py \\
        --sample ${sample_name} \\
        --query ${query_tsv} \\
        --taxonomy gtdb_taxonomy/* \\
        --species-ani ${params.sylph_species_ani} \\
        --borderline-ani ${params.sylph_borderline_ani} \\
        --min-margin ${params.sylph_min_margin} \\
        --output ${sample_name}_${out_suffix}.tsv
    """
}

process AMRFINDERPLUS {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly), val(organism)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_amrfinder.tsv"), emit: amrfinder_tsv

    script:
    def sample_name = "${sample_id}_${asm_type}"
    // Without --organism AMRFinderPlus reports acquired genes only: point mutations are organism-
    // specific and are simply not called.
    def org_arg = organism ? "--organism ${organism}" : ""
    """
    amrfinder --plus \
        --name ${sample_name} \
        -n ${assembly} \
        ${org_arg} \
        --threads ${task.cpus} \
        -o ${sample_name}_amrfinder.tsv
    """
}

process RESFINDER {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly), val(species)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_resfinder"), emit: resfinder_dir
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_pheno_table.txt"), emit: pheno_table
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_pointfinder.txt"), emit: pointfinder, optional: true
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_resfinder_hits.txt"), emit: gene_hits, optional: true

    script:
    def sample_name = "${sample_id}_${asm_type}"
    def outdir = "${sample_name}_resfinder"
    // PointFinder needs a supported species; "other" restricts the run to acquired genes, which is
    // the only thing the phenotype table needs anyway.
    def resolved = species ?: 'other'
    // PointFinder takes its own length/identity thresholds; without -l_p/-t_p it falls back to its
    // defaults rather than inheriting -l/-t.
    def point_arg = species ? '-c -l_p 0.6 -t_p 0.8' : ''
    """
    mkdir -p ${outdir}
    python3 -m resfinder \
        -ifa ${assembly} \
        -o ${outdir} \
        -acq -l 0.6 -t 0.8 \
        ${point_arg} \
        --species "${resolved}"

    cp ${outdir}/pheno_table.txt ${sample_name}_pheno_table.txt

    # ResFinder_results_tab.txt carries a Contig column holding the full FASTA header, so its
    # genes can be placed on a chromosome or plasmid the same way AMRFinderPlus's are.
    if [ -s "${outdir}/ResFinder_results_tab.txt" ]; then
        cp ${outdir}/ResFinder_results_tab.txt ${sample_name}_resfinder_hits.txt
    fi

    # PointFinder only runs with -c against a supported species, and only writes a result file
    # when it finds a mutation.
    if [ -s "${outdir}/PointFinder_results.txt" ]; then
        cp ${outdir}/PointFinder_results.txt ${sample_name}_pointfinder.txt
    fi
    """
}

process PARSE_AMRFINDER {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(amrfinder_tsv)

    output:
    path "${sample_id}_${asm_type}_amrfinder_parsed.tsv", emit: for_multiqc

    script:
    def sample_name = "${sample_id}_${asm_type}"
    """
    parse_amrfinder.py --sample ${sample_name} --input ${amrfinder_tsv} --output ${sample_name}_amrfinder_parsed.tsv
    """
}

process PARSE_RESFINDER {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(pheno_table)

    output:
    path "${sample_id}_${asm_type}_resfinder_parsed.tsv", emit: for_multiqc

    script:
    def sample_name = "${sample_id}_${asm_type}"
    """
    parse_resfinder.py --sample ${sample_name} --input ${pheno_table} --output ${sample_name}_resfinder_parsed.tsv
    """
}

process MOB_RECON {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_mob_recon"), emit: mob_dir
    path "${sample_id}_${asm_type}_contig_report.txt", emit: contig_report, optional: true

    script:
    def outdir = "${sample_id}_${asm_type}_mob_recon"
    """
    mob_recon \
        --infile ${assembly} \
        --outdir ${outdir} \
        --num_threads ${task.cpus}

    if [ -f "${outdir}/contig_report.txt" ]; then
        cp "${outdir}/contig_report.txt" "${sample_id}_${asm_type}_contig_report.txt"
    fi
    """
}

process PARSE_MOB_RECON {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(mob_dir)

    output:
    path "${sample_id}_${asm_type}_mob_parsed.tsv", emit: for_multiqc

    script:
    def sample_name = "${sample_id}_${asm_type}"
    """
    parse_mob_recon.py --sample ${sample_name} --mob-dir ${mob_dir} --output ${sample_name}_mob_parsed.tsv
    """
}

process VIRULENCEFINDER {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_virulencefinder"), emit: virfinder_dir

    script:
    def outdir = "${sample_id}_${asm_type}_virulencefinder"
    """
    mkdir -p ${outdir}
    python -m virulencefinder \
        -ifa ${assembly} \
        --extended_output \
        -o ${outdir}
    """
}

process PARSE_VIRULENCEFINDER {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/annotation/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(virfinder_dir)

    output:
    path "${sample_id}_${asm_type}_virulencefinder_parsed.tsv", emit: for_multiqc

    script:
    def sample_name = "${sample_id}_${asm_type}"
    """
    parse_virulencefinder.py --sample ${sample_name} --virfinder-dir ${virfinder_dir} --output ${sample_name}_virulencefinder_parsed.tsv
    """
}