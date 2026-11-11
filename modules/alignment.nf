process ALIGN_LONG {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"

    input:
    tuple val(sample_id), val(asm_type), path(assembly), path(long_reads)

    output:
    tuple val(sample_id), val(asm_type), path("aligned.sam"), emit: sam

    when:
    long_reads.name != 'NO_FILE'

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${assembly} ${long_reads} -o aligned.sam
    """
}

process ALIGN_SHORT {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"

    input:
    tuple val(sample_id), val(asm_type), path(assembly), path(r1), path(r2)

    output:
    tuple val(sample_id), val(asm_type), path("aligned.sam"), emit: sam

    when:
    r1.name != 'NO_FILE'

    script:
    """
    bwa-mem2 index ${assembly}
    bwa-mem2 mem -t ${task.cpus} ${assembly} ${r1} ${r2} > aligned.sam
    """
}

process SAMTOOLS_SORT {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/stats/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(sam)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}.bam"), path("${sample_id}_${asm_type}.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}_${asm_type}.bam ${sam}
    samtools index ${sample_id}_${asm_type}.bam
    """
}