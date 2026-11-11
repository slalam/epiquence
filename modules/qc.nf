process VALIDATE_READS {
    label 'process_low'
    tag "${sample_id}_${read_type}"

    input:
    tuple val(sample_id), val(read_type), path(reads)

    output:
    tuple val(sample_id), path(reads), path("${sample_id}_${read_type}.validation"), emit: validated

    script:
    def input_files = reads instanceof List ? reads.join(' ') : reads
    """
    STATUS="PASS"
    ERRORS=""
    for f in ${input_files}; do
        [ "\$(basename \$f)" = "NO_FILE" ] && continue

        if command -v pigz &>/dev/null; then
            if ! pigz -t "\$f" 2>/dev/null; then
                ERRORS="\${ERRORS}\$f: failed gzip integrity check\\n"
                STATUS="FAIL"
                continue
            fi
        else
            if ! gzip -t "\$f" 2>/dev/null; then
                ERRORS="\${ERRORS}\$f: failed gzip integrity check\\n"
                STATUS="FAIL"
                continue
            fi
        fi

        HAS_NULL=\$(zcat "\$f" 2>/dev/null | tr -cd '\\0' | head -c1 | wc -c | tr -d ' ')
        if [ "\$HAS_NULL" -gt 0 ]; then
            ERRORS="\${ERRORS}\$f: contains null bytes; file is corrupt\\n"
            STATUS="FAIL"
        # Check FASTQ structure with seqkit stats -a (catches bad format for Illumina)
        elif ! seqkit stats "\$f" -aT 2>/tmp/seqkit_err_\$\$.log | tail -n +2 | grep -q .; then
            ERR_MSG=\$(grep -iE 'error|bad' /tmp/seqkit_err_\$\$.log 2>/dev/null | head -1)
            ERRORS="\${ERRORS}\$f: FASTQ validation failed: \${ERR_MSG:-no valid reads found}\\n"
            STATUS="FAIL"
        fi
        rm -f /tmp/seqkit_err_\$\$.log
    done

    echo "\$STATUS" > ${sample_id}_${read_type}.validation
    if [ "\$STATUS" = "FAIL" ]; then
        echo "VALIDATION FAILED for ${sample_id} (${read_type}):" >&2
        printf "\$ERRORS" >&2
        printf "\$ERRORS" >> ${sample_id}_${read_type}.validation
    fi
    """
}

process CHECK_VALIDATION {
    label 'process_low'
    errorStrategy 'terminate'

    input:
    path(validation_files)

    output:
    path "validation_passed"

    script:
    """
    FAILED=0
    echo "========================================" >&2
    echo "  Read Validation Summary" >&2
    echo "========================================" >&2
    for f in ${validation_files}; do
        SAMPLE=\$(basename "\$f" .validation)
        STATUS=\$(head -1 "\$f")
        if [ "\$STATUS" = "FAIL" ]; then
            FAILED=1
            echo "" >&2
            echo "  FAIL: \$SAMPLE" >&2
            tail -n +2 "\$f" | while read -r line; do
                [ -n "\$line" ] && echo "    - \$line" >&2
            done
        fi
    done

    if [ \$FAILED -eq 1 ]; then
        echo "" >&2
        echo "========================================" >&2
        echo "  Pipeline aborted: fix the above files" >&2
        echo "  or remove them from the sample sheet." >&2
        echo "========================================" >&2
        exit 1
    fi

    echo "All samples passed validation" >&2
    echo "========================================" >&2
    touch validation_passed
    """
}

process SEQKIT_STATS {
    label 'process_low'
    tag "${sample_id}_${read_type}"
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), val(read_type), path(reads)

    output:
    tuple val(sample_id), val(read_type), path("${sample_id}_${read_type}_seqkit_stats.tsv"), emit: stats

    script:
    def input_files = reads instanceof List ? reads.join(' ') : reads
    """
    seqkit stats ${input_files} -j ${task.cpus} -abT > ${sample_id}_${read_type}_seqkit_stats.tsv
    """
}

process PARSE_SEQKIT_STATS {
    label 'process_low'
    tag "${sample_id}_${read_type}"
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), val(read_type), path(stats_file)

    output:
    path "${sample_id}_${read_type}_seqkit_mqc.tsv", emit: for_multiqc

    script:
    def sample_name = "${sample_id}_${read_type}"
    """
    parse_seqkit_stats.py --sample ${sample_name} --input ${stats_file} --output ${sample_name}_seqkit_mqc.tsv
    """
}

process NANOPLOT {
    label 'process_medium'
    tag "${sample_id}_L"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${sample_id}/qc/nanoplot", mode: 'copy'

    input:
    tuple val(sample_id), path(long_reads)

    output:
    tuple val(sample_id), path("NanoStats.txt"), emit: stats

    script:
    """
    NanoPlot --fastq ${long_reads} -o . --no_static
    """
}

process PARSE_NANOPLOT {
    label 'process_low'
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(nanoplot_stats)

    output:
    path "${sample_id}_L_quality_mqc.tsv", emit: for_multiqc

    script:
    def sample_name = "${sample_id}_L"
    """
    parse_nanoplot.py --sample ${sample_name} --input ${nanoplot_stats} --output ${sample_name}_quality_mqc.tsv
    """
}

process FASTPLONG {
    label 'process_medium'
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy', pattern: "*_fastplong.json"

    input:
    tuple val(sample_id), path(long_reads)

    output:
    tuple val(sample_id), path("${sample_id}_clean_long.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_fastplong.json"), emit: json

    script:
    """
    fastplong -i ${long_reads} -o ${sample_id}_clean_long.fastq.gz --length_required 1000 --json ${sample_id}_fastplong.json
    """
}

process RASUSA {
    label 'process_low'
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(long_reads)

    output:
    tuple val(sample_id), path("${sample_id}_subsampled_long.fastq.gz"), emit: reads

    script:
    """
    rasusa reads \
        --coverage ${params.rasusa_coverage} \
        --genome-size ${params.rasusa_genome_size} \
        -s 1 \
        -O g \
        -o ${sample_id}_subsampled_long.fastq.gz \
        ${long_reads}
    """
}

process FASTP_SHORT {
    label 'process_medium'
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy', pattern: "*_fastp.json"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1_clean.fastq.gz"), path("${sample_id}_R2_clean.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json

    script:
    """
    fastp -i ${r1} -I ${r2} \
        -o ${sample_id}_R1_clean.fastq.gz \
        -O ${sample_id}_R2_clean.fastq.gz \
        -j ${sample_id}_fastp.json \
        --thread ${task.cpus} --detect_adapter_for_pe
    """
}

process PARSE_FASTP {
    label 'process_low'
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(json_file)

    output:
    path "${sample_id}_S_quality_mqc.tsv", emit: for_multiqc

    script:
    def sample_name = "${sample_id}_S"
    """
    parse_fastp.py --sample ${sample_name} --input ${json_file} --output ${sample_name}_quality_mqc.tsv
    """
}

process PARSE_FASTPLONG {
    label 'process_low'
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(json_file)

    output:
    path "${sample_id}_L_filterstats_mqc.tsv", emit: for_multiqc

    script:
    def sample_name = "${sample_id}_L"
    """
    parse_fastplong.py --sample ${sample_name} --input ${json_file} --output ${sample_name}_filterstats_mqc.tsv
    """
}