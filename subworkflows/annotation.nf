include { MLST; MLST_WITH_SCHEME; PARSE_MLST; GTDBTK; PARSE_GTDBTK; SYLPH_QUERY; PARSE_SYLPH; AMRFINDERPLUS; RESFINDER; PARSE_AMRFINDER; PARSE_RESFINDER; MOB_RECON; PARSE_MOB_RECON; VIRULENCEFINDER; PARSE_VIRULENCEFINDER } from '../modules/annotation'
include { ANNOTATION_SUMMARY } from '../modules/summary'

// Resolve a GTDB-Tk species call to a tool-specific organism name.
def resolve_organism(species, map_) {
    if (!species) return ''
    def s = species.toLowerCase()
    def hit = map_.keySet().sort { -it.length() }.find { s.contains(it.toLowerCase()) }
    return hit ? map_[hit] : ''
}

workflow ANNOTATION {
    take:
    final_assemblies
    taxonomy_tool       // 'gtdbtk' | 'sylph' | 'both' (resolved and validated in main.nf)

    main:
    // Run initial MLST (auto-detect scheme) and species identification in parallel
    MLST(final_assemblies)

    // ── Species identification ──
    // Whichever engine is authoritative produces {sample_id}_{asm_type}_gtdbtk_parsed.tsv with the
    // same leading columns, so everything below and downstream is engine-agnostic.
    taxonomy_parsed     = channel.empty()
    taxonomy_comparison = channel.empty()
    gtdbtk_raw_summary  = channel.empty()

    if (taxonomy_tool in ['gtdbtk', 'both']) {
        GTDBTK(final_assemblies)
        PARSE_GTDBTK(GTDBTK.out.gtdbtk_dir)
        taxonomy_parsed    = PARSE_GTDBTK.out.parsed_taxonomy
        gtdbtk_raw_summary = GTDBTK.out.raw_summary
    }

    if (taxonomy_tool in ['sylph', 'both']) {
        sylph_taxonomy_files = channel.value(
            params.sylph_taxonomy.toString().tokenize(',').collect { p -> file(p.trim()) })
        SYLPH_QUERY(final_assemblies)
        PARSE_SYLPH(SYLPH_QUERY.out.hits, sylph_taxonomy_files,
                    taxonomy_tool == 'sylph' ? 'gtdbtk_parsed' : 'sylph_parsed')
        if (taxonomy_tool == 'sylph') {
            taxonomy_parsed = PARSE_SYLPH.out.parsed_taxonomy
        } else {
            taxonomy_comparison = PARSE_SYLPH.out.parsed_taxonomy
        }
    }

    // Species per sample from the authoritative taxonomy file.
    gtdbtk_species_ch = taxonomy_parsed
        .map { f ->
            def lines = f.readLines()
            if (lines.size() > 1) {
                def parts = lines[1].split('\t')
                def sample_full = parts[0]       // e.g. SAMPLE-04_HL
                def species = parts[7] ?: ''      // species column (8th, 0-indexed 7)
                def sample_parts = sample_full.split('_')
                def asm_type = sample_parts[-1]
                def sample_id = sample_parts[0..-2].join('_')
                return [sample_id, asm_type, species]
            }
            return null
        }
        .filter { it != null }

    // ── MLST scheme override ──
    // Join assemblies with GTDB species and initial MLST results, and determine which samples need
    // MLST re-run with a forced scheme.
    assemblies_keyed = final_assemblies.map { id, asm, assembly -> [id, asm, assembly] }

    needs_rerun_ch = assemblies_keyed
        .join(gtdbtk_species_ch, by: [0, 1])      // [sample_id, asm_type, assembly, species]
        .join(MLST.out.mlst_result, by: [0, 1])   // + mlst_file
        .map { sample_id, asm_type, assembly, species, mlst_file ->
            def correct_scheme = null
            params.species_scheme_map.each { sp, scheme ->
                if (species?.toLowerCase()?.contains(sp.toLowerCase())) {
                    correct_scheme = scheme
                }
            }
            if (correct_scheme) {
                def mlst_content = mlst_file.text.trim()
                def mlst_parts = mlst_content.split('\t')
                def current_scheme = mlst_parts.size() > 1 ? mlst_parts[1] : ''
                if (current_scheme != correct_scheme) {
                    return [sample_id, asm_type, assembly, correct_scheme]
                }
            }
            return null
        }
        .filter { it != null }

    rerun_keys = needs_rerun_ch.map { id, asm, assembly, scheme -> [id, asm] }

    MLST_WITH_SCHEME(needs_rerun_ch)

    original_mlst = MLST.out.mlst_result
        .join(rerun_keys.map { id, asm -> [id, asm, true] }, by: [0, 1], remainder: true)
        .filter { it[-1] == null }  // Keep only those NOT in rerun list
        .map { it[0..2] }           // [sample_id, asm_type, mlst_file]

    all_mlst = original_mlst.mix(MLST_WITH_SCHEME.out.mlst_result)

    PARSE_MLST(all_mlst)

    // ── Resistance ──
    // Both tools are organism-aware, so they consume the assembly joined to its GTDB-Tk species
    // rather than the bare assembly channel.
    assemblies_with_species = final_assemblies.join(gtdbtk_species_ch, by: [0, 1])

    amrfinder_in = assemblies_with_species.map { sample_id, asm_type, assembly, species ->
        [sample_id, asm_type, assembly, resolve_organism(species, params.amrfinder_organism_map)]
    }
    resfinder_in = assemblies_with_species.map { sample_id, asm_type, assembly, species ->
        [sample_id, asm_type, assembly, resolve_organism(species, params.resfinder_species_map)]
    }

    AMRFINDERPLUS(amrfinder_in)
    RESFINDER(resfinder_in)
    PARSE_AMRFINDER(AMRFINDERPLUS.out.amrfinder_tsv)
    PARSE_RESFINDER(RESFINDER.out.pheno_table)

    MOB_RECON(final_assemblies)
    PARSE_MOB_RECON(MOB_RECON.out.mob_dir)
    VIRULENCEFINDER(final_assemblies)
    PARSE_VIRULENCEFINDER(VIRULENCEFINDER.out.virfinder_dir)

    mlst_files       = all_mlst.map { id, asm_type, f -> f }.collect().ifEmpty(file('NO_MLST'))
    amrfinder_files  = AMRFINDERPLUS.out.amrfinder_tsv.map { id, asm_type, f -> f }.collect().ifEmpty(file('NO_AMRFINDER'))
    pheno_files      = RESFINDER.out.pheno_table.map { id, asm_type, f -> f }.collect().ifEmpty(file('NO_PHENO'))
    mob_reports      = MOB_RECON.out.contig_report.collect().ifEmpty(file('NO_MOB'))
    gtdbtk_summaries = gtdbtk_raw_summary.flatten().collect().ifEmpty(file('NO_GTDBTK'))

    ANNOTATION_SUMMARY(mlst_files, amrfinder_files, pheno_files, mob_reports, gtdbtk_summaries)

    emit:
    mlst_for_multiqc      = PARSE_MLST.out.for_multiqc
    gtdbtk_parsed         = taxonomy_parsed          // authoritative engine, *_gtdbtk_parsed.tsv
    taxonomy_comparison   = taxonomy_comparison      // 'both' mode only, *_sylph_parsed.tsv
    amrfinder_for_multiqc = PARSE_AMRFINDER.out.for_multiqc
    resfinder_for_multiqc = PARSE_RESFINDER.out.for_multiqc
    amrfinder_tsv         = AMRFINDERPLUS.out.amrfinder_tsv   // [id, asm_type, tsv]
    pheno_table           = RESFINDER.out.pheno_table         // [id, asm_type, txt]
    pointfinder           = RESFINDER.out.pointfinder         // [id, asm_type, txt] (optional)
    resfinder_hits        = RESFINDER.out.gene_hits           // [id, asm_type, txt] (optional)
    mob_for_multiqc       = PARSE_MOB_RECON.out.for_multiqc
    mob_recon_dir         = MOB_RECON.out.mob_dir
    mob_raw_dir           = MOB_RECON.out.mob_dir
    virfinder_for_multiqc = PARSE_VIRULENCEFINDER.out.for_multiqc
    virfinder_raw_dir     = VIRULENCEFINDER.out.virfinder_dir
}