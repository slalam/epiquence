include { ALIGN_LONG; ALIGN_SHORT; SAMTOOLS_SORT } from '../modules/alignment'
include { MOSDEPTH; PARSE_MOSDEPTH; QUAST; PARSE_QUAST; BUSCO } from '../modules/stats'

workflow ALIGN_AND_STATS {
    take:
    final_assemblies
    clean_reads

    main:
    asm_with_reads = final_assemblies
        .map { id, asm_type, asm -> [id, [asm_type, asm]] }
        .combine(clean_reads.map { [it[0], it[1], it[2], it[3]] }, by: 0)
        .map { id, asm_data, long_r, r1, r2 -> [id, asm_data[0], asm_data[1], long_r, r1, r2] }

    long_read_coverage = asm_with_reads
        .filter { it[1] == 'L' || it[1] == 'HL' }
        .map { id, asm_type, asm, long_r, r1, r2 -> [id, asm_type, asm, long_r] }

    short_read_coverage = asm_with_reads
        .filter { it[1] == 'S' || it[1] == 'HS' }
        .map { id, asm_type, asm, long_r, r1, r2 -> [id, asm_type, asm, r1, r2] }

    ALIGN_LONG(long_read_coverage)
    ALIGN_SHORT(short_read_coverage)
    SAMTOOLS_SORT(ALIGN_LONG.out.sam.mix(ALIGN_SHORT.out.sam))
    MOSDEPTH(SAMTOOLS_SORT.out.bam)
    PARSE_MOSDEPTH(MOSDEPTH.out.stats)
    QUAST(final_assemblies)
    PARSE_QUAST(QUAST.out.quast_dir)
    BUSCO(final_assemblies)

    emit:
    mosdepth_for_multiqc = PARSE_MOSDEPTH.out.for_multiqc
    mosdepth_dist        = MOSDEPTH.out.stats
    quast                = PARSE_QUAST.out.quast_output
    quast_for_multiqc    = PARSE_QUAST.out.for_multiqc
    busco                = BUSCO.out.for_multiqc
}