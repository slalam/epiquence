process MOSDEPTH {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/stats/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(bam), path(bai)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}.mosdepth.global.dist.txt"), path("${sample_id}_${asm_type}.mosdepth.summary.txt"), emit: stats

    script:
    """
    mosdepth -t 4 -n ${sample_id}_${asm_type} ${bam}
    """
}

process PARSE_MOSDEPTH {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/stats/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(global_dist), path(summary)

    output:
    path "${sample_id}_${asm_type}_coverage_mqc.tsv", emit: for_multiqc

    script:
    def prefix = "${sample_id}_${asm_type}"
    """
    parse_mosdepth.py --sample ${prefix} --global-dist ${global_dist} --summary ${summary} --output ${prefix}_coverage_mqc.tsv
    """
}

process QUAST {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/stats/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_quast"), emit: quast_dir

    script:
    """
    quast.py ${assembly} -o ${sample_id}_${asm_type}_quast -t ${task.cpus}
    """
}

process PARSE_QUAST {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/stats/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(quast_dir)

    output:
    path "${sample_id}_${asm_type}_quast", emit: quast_output
    path "${sample_id}_${asm_type}_quast_mqc.tsv", emit: for_multiqc

    script:
    def sample_name = "${sample_id}_${asm_type}"
    """
    parse_quast.py --sample ${sample_name} --quast-dir ${quast_dir} --output ${sample_name}_quast_mqc.tsv
    """
}

process BUSCO {
    label 'process_high'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/stats/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    path "${sample_id}_${asm_type}_busco", emit: busco_dir
    path "short_summary.specific.bacteria_odb12.${sample_id}_${asm_type}.txt", emit: for_multiqc

    script:
    def outdir = "${sample_id}_${asm_type}_busco"
    """
    busco -i ${assembly} \
        -o ${outdir} \
        -m genome \
        -l ${params.busco_db}/lineages/bacteria_odb12 \
        -c ${task.cpus} \
        --offline

    cp ${outdir}/short_summary.*.txt short_summary.specific.bacteria_odb12.${sample_id}_${asm_type}.txt
    """
}