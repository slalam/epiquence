process FLYE {
    label 'process_high'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(long_reads)

    output:
    tuple val(sample_id), path("assembly.fasta"), emit: assembly

    script:
    """
    flye --nano-hq ${long_reads} --out-dir . --threads ${task.cpus}
    """
}

process MEDAKA {
    label 'process_medium'
    tag "${sample_id}"
    stageInMode 'copy'

    input:
    tuple val(sample_id), path(assembly), path(long_reads)

    output:
    tuple val(sample_id), path("${sample_id}_medaka.fasta"), emit: assembly

    script:
    """
    # Remove any stale index files that may have been copied
    rm -f ${assembly}.fai ${assembly}.map-ont.mmi
    medaka_consensus -i ${long_reads} -d ${assembly} -o . -t ${task.cpus} -m ${params.medaka_model} --bacteria
    mv consensus.fasta ${sample_id}_medaka.fasta
    """
}

process UNICYCLER {
    label 'process_max'
    tag "${sample_id}_${asm_type}"

    input:
    tuple val(sample_id), val(asm_type), path(r1), path(r2), path(long_reads)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}_asm.fasta"), emit: assembly

    script:
    def long_opt = long_reads.name != 'NO_FILE' ? "-l ${long_reads}" : ""
    """
    unicycler -1 ${r1} -2 ${r2} ${long_opt} -o . -t ${task.cpus}
    mv assembly.fasta ${sample_id}_${asm_type}_asm.fasta
    """
}

process POLYPOLISH_ALIGN {
    label 'process_high'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(assembly), path(r1), path(r2)

    output:
    tuple val(sample_id), path(assembly), path("R1.sam"), path("R2.sam"), emit: aligned

    script:
    """
    bwa-mem2 index ${assembly}
    bwa-mem2 mem -t ${task.cpus} -a ${assembly} ${r1} > R1.sam
    bwa-mem2 mem -t ${task.cpus} -a ${assembly} ${r2} > R2.sam
    """
}

process POLYPOLISH {
    label 'process_high'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(assembly), path(sam1), path(sam2)

    output:
    tuple val(sample_id), path("${sample_id}_polypolish.fasta"), emit: assembly

    script:
    """
    polypolish filter --in1 ${sam1} --in2 ${sam2} --out1 R1_filt.sam --out2 R2_filt.sam
    polypolish polish ${assembly} R1_filt.sam R2_filt.sam > ${sample_id}_polypolish.fasta
    """
}

process DNAAPLER {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/assemblies", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}.fasta"), emit: assembly

    script:
    """
    set +e
    dnaapler all -i ${assembly} -o dnaapler_out -t ${task.cpus} -f
    set -e

    if [ -s "dnaapler_out/dnaapler_reoriented.fasta" ]; then
        mv dnaapler_out/dnaapler_reoriented.fasta "${sample_id}_${asm_type}.fasta"
    else
        cp ${assembly} "${sample_id}_${asm_type}.fasta"
    fi
    """
}

process RENAME_CONTIGS {
    label 'process_low'
    tag "${sample_id}_${asm_type}"

    input:
    tuple val(sample_id), val(asm_type), path(assembly, stageAs: 'input_assembly.fasta')

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}.fasta"), emit: assembly

    script:
    """
    rename_contigs.py --input input_assembly.fasta --output ${sample_id}_${asm_type}.fasta
    """
}

process ANNOTATE_PLASMIDS {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/${sample_id}/assemblies", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(asm_type), path(assembly, stageAs: 'input_assembly.fasta'), path(mob_dir)

    output:
    tuple val(sample_id), val(asm_type), path("${sample_id}_${asm_type}.fasta"), emit: assembly

    script:
    """
    annotate_plasmids.py --input input_assembly.fasta --mob-dir ${mob_dir} --output ${sample_id}_${asm_type}.fasta
    """
}

process VALIDATE_ASSEMBLY {
    label 'process_low'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path(assembly), emit: validated

    script:
    """
    # Check file is valid FASTA with at least one contig
    N_CONTIGS=\$(grep -c "^>" ${assembly} || true)
    if [ "\$N_CONTIGS" -eq 0 ]; then
        echo "ERROR: ${assembly} contains no contigs" >&2
        exit 1
    fi
    TOTAL_BP=\$(grep -v "^>" ${assembly} | tr -d '\\n\\r ' | wc -c)
    if [ "\$TOTAL_BP" -lt 100 ]; then
        echo "ERROR: ${assembly} too small (\$TOTAL_BP bp)" >&2
        exit 1
    fi
    echo "Validated: ${assembly} | \$N_CONTIGS contigs, \$TOTAL_BP bp"
    """
}

process DNAAPLER_ASSEMBLY {
    label 'process_low'
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/assemblies", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), val("A"), path("${sample_id}_A.fasta"), emit: assembly

    script:
    """
    # Check if already rotated (dnaapler or other markers in headers)
    NEEDS_ROTATION=1
    if grep -qiE 'rotated=True|rotated_gene=|dnaapler|\\[circularised\\]' ${assembly}; then
        echo "Assembly appears already rotated, skipping dnaapler"
        NEEDS_ROTATION=0
    fi

    if [ "\$NEEDS_ROTATION" -eq 1 ]; then
        set +e
        dnaapler all -i ${assembly} -o dnaapler_out -t ${task.cpus} -f
        set -e
        if [ -s "dnaapler_out/dnaapler_reoriented.fasta" ]; then
            cp dnaapler_out/dnaapler_reoriented.fasta ${sample_id}_A.fasta
        else
            cp ${assembly} ${sample_id}_A.fasta
        fi
    else
        cp ${assembly} ${sample_id}_A.fasta
    fi
    """
}