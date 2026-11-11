/* Species-specific typing; add a new tool as a process here */

process ECTYPER {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/species_typing/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_ectyper"), emit: ectyper_dir

    script:
    def outdir = "${sample_id}_${asm_type}_ectyper"
    """
    ectyper \
        -i ${assembly} \
        -c ${task.cpus} \
        --pathotype \
        --verify \
        -o ${outdir}
    """
}

process KLEBORATE {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/species_typing/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly), val(preset)

    output:
    tuple val(sample_id), val(asm_type), val(preset), path("${sample_id}_${asm_type}_kleborate"), emit: kleborate_dir

    script:
    def outdir = "${sample_id}_${asm_type}_kleborate"
    """
    mkdir -p ${outdir}
    kleborate \
        -a ${assembly} \
        --preset ${preset} \
        --trim_headers \
        -o ${outdir}
    """
}

process EMMTYPER {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/species_typing/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_emmtyper.tsv"), emit: emmtyper_result

    script:
    """
    emmtyper ${assembly} -o ${sample_id}_${asm_type}_emmtyper.tsv -f verbose
    """
}
process SISTR {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/species_typing/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_sistr"), emit: sistr_dir

    script:
    def outdir = "${sample_id}_${asm_type}_sistr"
    """
    mkdir -p ${outdir}
    sistr --qc -vv -f csv -o ${outdir}/sistr-output ${assembly}
    """
}

process SEROBA {
    // Streptococcus pneumoniae serotyping from Illumina reads (short reads only).
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/species_typing/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(reads1), path(reads2)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_seroba"), emit: seroba_dir

    script:
    def prefix = "${sample_id}_${asm_type}_seroba"
    """
    # SeroBA matches the forward/reverse pair on FILENAMES and requires the differentiator (_1/_2)
    # immediately before the .fastq.gz extension.
    ln -sf ${reads1} ${sample_id}_1.fastq.gz
    ln -sf ${reads2} ${sample_id}_2.fastq.gz

    # Keep matplotlib's cache off the (read-only) container HOME.
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    mkdir -p "\$MPLCONFIGDIR"

    seroba runSerotyping /seroba/database/ ${sample_id}_1.fastq.gz ${sample_id}_2.fastq.gz ${prefix}
    """
}

process MENINGOTYPE {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/species_typing/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_meningotype.tsv"), emit: meningotype_result

    script:
    """
    meningotype --finetype ${assembly} > ${sample_id}_${asm_type}_meningotype.tsv
    """
}

process ELGATO {
    label 'process_medium'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/species_typing/${asm_type}", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_elgato"), emit: elgato_dir

    script:
    def outdir = "${sample_id}_${asm_type}_elgato"
    """
    el_gato.py --assembly ${assembly} --out ${outdir}
    """
}