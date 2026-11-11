include { FLYE; MEDAKA; UNICYCLER; POLYPOLISH_ALIGN; POLYPOLISH } from '../modules/assembly'

workflow LONG_READ_ASSEMBLY {
    take:
    clean_reads

    main:
    long_in = clean_reads.filter { it[1].name != 'NO_FILE' }.map { [it[0], it[1]] }
    FLYE(long_in)
    medaka_in = FLYE.out.assembly.join(long_in).map { id, asm, long_r -> [id, asm, long_r] }
    MEDAKA(medaka_in)

    emit:
    assembly      = MEDAKA.out.assembly
    flye_assembly = FLYE.out.assembly
    long_reads    = long_in
}

workflow SHORT_READ_ASSEMBLY {
    take:
    clean_reads

    main:
    short_in = clean_reads.filter { it[2].name != 'NO_FILE' }.map { id, l, r1, r2 -> [id, "S", r1, r2, file('NO_FILE')] }
    UNICYCLER(short_in)

    emit:
    assembly = UNICYCLER.out.assembly.map { id, asm_type, asm -> [id, asm] }
}

workflow HYBRID_LONG_FIRST {
    take:
    flye_assembly
    long_reads
    clean_reads

    main:
    // Run MEDAKA (fasta) for polishing before Polypolish
    medaka_in = flye_assembly.join(long_reads).map { id, asm, long_r -> [id, asm, long_r] }
    MEDAKA(medaka_in)
    poly_in = MEDAKA.out.assembly.join(clean_reads).filter { it[3].name != 'NO_FILE' }.map { id, asm, l, r1, r2 -> [id, asm, r1, r2] }
    POLYPOLISH_ALIGN(poly_in)
    POLYPOLISH(POLYPOLISH_ALIGN.out.aligned)

    emit:
    assembly = POLYPOLISH.out.assembly
}

workflow HYBRID_SHORT_FIRST {
    take:
    clean_reads

    main:
    hybrid_in = clean_reads.filter { it[1].name != 'NO_FILE' && it[2].name != 'NO_FILE' }.map { id, l, r1, r2 -> [id, "HS", r1, r2, l] }
    UNICYCLER(hybrid_in)

    emit:
    assembly = UNICYCLER.out.assembly.map { id, asm_type, asm -> [id, asm] }
}