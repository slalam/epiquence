/* Routes each assembly to typing tools by its GTDB species */

include { ECTYPER; KLEBORATE; EMMTYPER; SISTR; SEROBA; MENINGOTYPE; ELGATO } from '../modules/species_typing'

workflow SPECIES_TYPING {
    take:
    assemblies_with_species  // [sample_id, asm_type, assembly, species]
    short_reads              // [sample_id, reads1, reads2]  (for reads-based tools; NO_FILE when absent)

    main:
    // Species lists from config (comma-separated strings)
    def ectyper_species      = (params.ectyper_species ?: '').tokenize(',').collect { it.trim().toLowerCase() }
    def kpsc_species         = (params.kleborate_kpsc_species ?: '').tokenize(',').collect { it.trim().toLowerCase() }
    def kosc_species         = (params.kleborate_kosc_species ?: '').tokenize(',').collect { it.trim().toLowerCase() }
    def emmtyper_species     = (params.emmtyper_species ?: '').tokenize(',').collect { it.trim().toLowerCase() }
    def sistr_species        = (params.sistr_species ?: '').tokenize(',').collect { it.trim().toLowerCase() }
    def seroba_species       = (params.seroba_species ?: '').tokenize(',').collect { it.trim().toLowerCase() }
    def meningotype_species  = (params.meningotype_species ?: '').tokenize(',').collect { it.trim().toLowerCase() }
    def elgato_species       = (params.elgato_species ?: '').tokenize(',').collect { it.trim().toLowerCase() }

    // Assembly types that carry Illumina/short reads (SeroBA is reads-based, Illumina-only).
    def short_read_asm_types = ['S', 'HS', 'HL']

    // Route samples to tools based on species
    ectyper_in = assemblies_with_species
        .filter { sample_id, asm_type, assembly, species ->
            ectyper_species.any { sp -> species?.toLowerCase()?.contains(sp) }
        }
        .map { sample_id, asm_type, assembly, species -> [sample_id, asm_type, assembly] }

    kleborate_kpsc_in = assemblies_with_species
        .filter { sample_id, asm_type, assembly, species ->
            kpsc_species.any { sp -> species?.toLowerCase()?.contains(sp) }
        }
        .map { sample_id, asm_type, assembly, species -> [sample_id, asm_type, assembly, 'kpsc'] }

    kleborate_kosc_in = assemblies_with_species
        .filter { sample_id, asm_type, assembly, species ->
            kosc_species.any { sp -> species?.toLowerCase()?.contains(sp) }
        }
        .map { sample_id, asm_type, assembly, species -> [sample_id, asm_type, assembly, 'kosc'] }

    emmtyper_in = assemblies_with_species
        .filter { sample_id, asm_type, assembly, species ->
            emmtyper_species.any { sp -> species?.toLowerCase()?.contains(sp) }
        }
        .map { sample_id, asm_type, assembly, species -> [sample_id, asm_type, assembly] }

    sistr_in = assemblies_with_species
        .filter { sample_id, asm_type, assembly, species ->
            sistr_species.any { sp -> species?.toLowerCase()?.contains(sp) }
        }
        .map { sample_id, asm_type, assembly, species -> [sample_id, asm_type, assembly] }

    meningotype_in = assemblies_with_species
        .filter { sample_id, asm_type, assembly, species ->
            meningotype_species.any { sp -> species?.toLowerCase()?.contains(sp) }
        }
        .map { sample_id, asm_type, assembly, species -> [sample_id, asm_type, assembly] }

    elgato_in = assemblies_with_species
        .filter { sample_id, asm_type, assembly, species ->
            elgato_species.any { sp -> species?.toLowerCase()?.contains(sp) }
        }
        .map { sample_id, asm_type, assembly, species -> [sample_id, asm_type, assembly] }

    // SeroBA: reads-based, Illumina only.
    seroba_in = assemblies_with_species
        .filter { sample_id, asm_type, assembly, species ->
            seroba_species.any { sp -> species?.toLowerCase()?.contains(sp) } &&
            short_read_asm_types.contains(asm_type)
        }
        .map { sample_id, asm_type, assembly, species -> [sample_id, asm_type] }
        .combine(short_reads, by: 0)   // [sample_id, asm_type, reads1, reads2]
        .filter { sample_id, asm_type, r1, r2 -> r1.name != 'NO_FILE' && r2.name != 'NO_FILE' }

    ECTYPER(ectyper_in)
    KLEBORATE(kleborate_kpsc_in.mix(kleborate_kosc_in))
    EMMTYPER(emmtyper_in)
    SISTR(sistr_in)
    SEROBA(seroba_in)
    MENINGOTYPE(meningotype_in)
    ELGATO(elgato_in)

    emit:
    ectyper_dirs        = ECTYPER.out.ectyper_dir              // [sample_id, asm_type, dir]
    kleborate_dirs      = KLEBORATE.out.kleborate_dir          // [sample_id, asm_type, preset, dir]
    emmtyper_results    = EMMTYPER.out.emmtyper_result         // [sample_id, asm_type, tsv]
    sistr_dirs          = SISTR.out.sistr_dir                  // [sample_id, asm_type, dir]
    seroba_dirs         = SEROBA.out.seroba_dir                // [sample_id, asm_type, dir]
    meningotype_results = MENINGOTYPE.out.meningotype_result   // [sample_id, asm_type, tsv]
    elgato_dirs         = ELGATO.out.elgato_dir                // [sample_id, asm_type, dir]
}