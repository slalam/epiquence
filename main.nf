nextflow.enable.dsl=2

/*
 * EPIquence | Sequence to Surveillance
 * Entry workflow: inputs, species engine and run setup, module wiring.
 */

include { VALIDATE_READS; CHECK_VALIDATION; SEQKIT_STATS; PARSE_SEQKIT_STATS; NANOPLOT; PARSE_NANOPLOT; FASTPLONG; RASUSA; FASTP_SHORT; PARSE_FASTP; PARSE_FASTPLONG } from './modules/qc'
include { DNAAPLER; RENAME_CONTIGS; RENAME_CONTIGS as RENAME_CONTIGS_ASM; ANNOTATE_PLASMIDS; VALIDATE_ASSEMBLY; DNAAPLER_ASSEMBLY } from './modules/assembly'
include { COLLECT_ASSEMBLIES; COLLECT_PLASMIDS; COLLECT_CHROMOSOMES_SUMMARY; GENERATE_MULTIQC_CONFIG; MULTIQC } from './modules/summary'
include { BUILD_HGT_REPORT; RENDER_HGT_TREES; UMAP_PLASMIDS; AMR_PLASMID_MATRIX } from './modules/plasmid_hgt'
include { LONG_READ_ASSEMBLY; SHORT_READ_ASSEMBLY; HYBRID_LONG_FIRST; HYBRID_SHORT_FIRST } from './subworkflows/assembly'
include { ALIGN_AND_STATS } from './subworkflows/align_and_stats'
include { ANNOTATION } from './subworkflows/annotation'

include { OUTBREAK_DETECTION } from './subworkflows/outbreak'
include { CORE_GENOME_ANALYSIS } from './subworkflows/core_genome'
include { BUILD_OUTBREAK_REPORT } from './modules/outbreak'
include { PLASMID_HGT_DETECTION } from './subworkflows/plasmid_hgt'
include { STORE_RESULTS_DB; FINALIZE_WORKFLOW_LOG; METABASE_RESCAN } from './modules/database'
include { SURVEILLANCE } from './subworkflows/surveillance'
include { SPECIES_TYPING } from './subworkflows/species_typing'
include { runSettingsPath; writeRunSettings; deleteRunSettings; urlHasUserInfo } from './modules/utils'
include { taxonomyTool } from './modules/annotation'

workflow {
    def asm_input = params.assembly_type.toString().trim().toUpperCase()
    def asm_types = (asm_input == 'ALL' || asm_input == 'HYBRID')
        ? ['L', 'S', 'HL', 'HS'] as Set
        : asm_input.split(',').collect { it.trim().toUpperCase() } as Set

    def valid_types = ['L', 'S', 'HL', 'HS'] as Set
    def invalid = asm_types - valid_types
    if (invalid) {
        error "Invalid assembly type(s): ${invalid.join(', ')}. Valid values: L, S, HL, HS, ALL"
    }

    def run_L  = 'L'  in asm_types
    def run_S  = 'S'  in asm_types
    def run_HL = 'HL' in asm_types
    def run_HS = 'HS' in asm_types
    def need_long_assembly = run_L || run_HL

    // Check if assembly inputs are provided
    def has_assemblies = params.assemblies_dir && file(params.assemblies_dir).exists()
    def has_reads = params.sample_sheet && file(params.sample_sheet).exists()

    if (!has_reads && !has_assemblies) {
        error "Either --sample_sheet or --assemblies_dir (or both) must be provided"
    }

    // Count all samples (reads + assemblies)
    def read_sample_count = has_reads ? file(params.sample_sheet).readLines().size() - 1 : 0
    def asm_sample_count = 0
    if (has_assemblies) {
        def asm_dir = file(params.assemblies_dir)
        asm_sample_count = asm_dir.listFiles().findAll { it.name =~ /\.(fasta|fna|fa)$/ }.size()
    }
    def sample_count = read_sample_count + asm_sample_count

    // Add 'A' to asm_types if assemblies are provided
    if (has_assemblies && asm_sample_count > 0) {
        asm_types = asm_types + ['A'] as Set
    }

    // --- Species identification engine --- Validated before any work starts: a wrong path here
    // would otherwise only surface after QC and assembly, hours into the run.
    def taxonomy_tool = taxonomyTool(params)
    if (!(taxonomy_tool in ['gtdbtk', 'sylph', 'both'])) {
        error "taxonomy_tool must be 'gtdbtk', 'sylph' or 'both' (got '${taxonomy_tool}')"
    }
    if (taxonomy_tool in ['gtdbtk', 'both'] && !(params.gtdbtk_db && file(params.gtdbtk_db).exists())) {
        error "taxonomy_tool '${taxonomy_tool}' needs gtdbtk_db; '${params.gtdbtk_db}' does not exist"
    }
    if (taxonomy_tool in ['sylph', 'both']) {
        if (!(params.sylph_db && file(params.sylph_db).exists())) {
            error "taxonomy_tool '${taxonomy_tool}' needs sylph_db (a .syldb file); '${params.sylph_db}' does not exist"
        }
        if (!params.sylph_taxonomy) {
            error "taxonomy_tool '${taxonomy_tool}' needs sylph_taxonomy (GTDB taxonomy TSVs, comma-separated)"
        }
        params.sylph_taxonomy.toString().tokenize(',').each { t ->
            if (!file(t.trim()).exists()) {
                error "sylph_taxonomy file '${t.trim()}' does not exist"
            }
        }
    }
    // GTDB-Tk's placement step needs ~90+ GB of RAM; the low_capacity profile reserves far less,
    // so the run would be killed mid-classification.
    if (taxonomy_tool in ['gtdbtk', 'both'] && workflow.profile?.tokenize(',')?.contains('low_capacity')) {
        log.warn "taxonomy_tool '${taxonomy_tool}' runs GTDB-Tk, which needs ~90+ GB of RAM; " +
                 "the low_capacity profile is meant for taxonomy_tool = 'sylph'."
    }

    log.info """
    ========================================
    EPIquence v${workflow.manifest.version}
    Species engine:      ${taxonomy_tool}${workflow.profile && workflow.profile != 'standard' ? '   (profile: ' + workflow.profile + ')' : ''}
    Assembly types:      ${asm_types.sort().join(', ')}
    Read samples:        ${read_sample_count}
    Assembly samples:    ${asm_sample_count}
    Total samples:       ${sample_count}
    Outbreak detection:  ${params.outbreak}
    Plasmid HGT:         ${params.plasmid_hgt}
    HGT Maps:            ${params.hgt_maps}
    Metadata:            ${params.metadata ?: 'none'}
    ========================================
    """.stripIndent()

    log.info "  Sample count:      ${sample_count}"

    // Database and Metabase settings go to a private per-run file; tasks only get its path.
    if (params.surveillance) {
        def db_url  = params.surveillance_db_url ? params.surveillance_db_url.toString() : ''
        def api_key = params.metabase_api_key ? params.metabase_api_key.toString() : ''
        def settings_path = runSettingsPath(workflow.launchDir)
        writeRunSettings(settings_path, db_url, api_key)
        // Remove it when the run ends, successful or not; by then every database/Metabase task has
        // finished.
        workflow.onComplete {
            deleteRunSettings(settings_path)
        }

        // conf/env.config now holds secrets: make sure only the owner can read it.
        def env_cfg = file("${projectDir}/conf/env.config")
        if ((urlHasUserInfo(db_url) || api_key) && env_cfg.exists()) {
            try {
                def perms = java.nio.file.Files.getPosixFilePermissions(env_cfg)
                def loose = perms.findAll { perm -> perm.name().startsWith('GROUP') || perm.name().startsWith('OTHERS') }
                if (loose) {
                    log.warn "conf/env.config contains a password or API key and is readable by group/others; run: chmod 600 ${env_cfg}"
                }
            } catch (UnsupportedOperationException _e) {
                // Non-POSIX filesystem: nothing to check.
            }
        }
    }

    metadata_ch = params.metadata ? file(params.metadata) : file('NO_METADATA')

    // ── Read-based samples ──
    if (has_reads) {
        samples_ch = Channel.fromPath(params.sample_sheet)
            .splitCsv(header: true, strip: true)
            .map { row ->
                def id       = row.sample_id
                def long_r   = row.long_R?.trim()   ? file(row.long_R)   : file('NO_FILE')
                def short_r1 = row.short_R1?.trim() ? file(row.short_R1) : file('NO_FILE')
                def short_r2 = row.short_R2?.trim() ? file(row.short_R2) : file('NO_FILE')
                [id, long_r, short_r1, short_r2]
            }
    } else {
        samples_ch = Channel.empty()
    }

    // ── Assembly-only samples ──
    if (has_assemblies) {
        assembly_input_ch = Channel.fromPath("${params.assemblies_dir}/*.{fasta,fna,fa}")
            .map { fasta ->
                def base = fasta.baseName.replaceAll(/\.(fasta|fna|fa)$/, '')
                [base, fasta]
            }

        VALIDATE_ASSEMBLY(assembly_input_ch)
        DNAAPLER_ASSEMBLY(VALIDATE_ASSEMBLY.out.validated)
        RENAME_CONTIGS_ASM(DNAAPLER_ASSEMBLY.out.assembly)
    }

    long_ch  = samples_ch.filter { it[1].name != 'NO_FILE' }.map { [it[0], it[1]] }
    short_ch = samples_ch.filter { it[2].name != 'NO_FILE' }.map { [it[0], it[2], it[3]] }

    // Determine which read types are needed
    def need_long  = run_L || run_HL || run_HS
    def need_short = run_S || run_HL || run_HS

    // Fork channels for multiple consumers
    long_ch.multiMap { id, reads ->
        for_validate: [id, "L", reads]
        for_seqkit: [id, "L", reads]
        for_nanoplot: [id, reads]
        for_fastplong: [id, reads]
    }.set { long_forked }

    short_ch.multiMap { id, r1, r2 ->
        for_validate: [id, "S", [r1, r2]]
        for_seqkit: [id, "S", [r1, r2]]
        for_fastp: [id, r1, r2]
    }.set { short_forked }

    // --- Validate all read files before processing ---
    validate_input = Channel.empty()
    if (need_long)  { validate_input = validate_input.mix(long_forked.for_validate) }
    if (need_short) { validate_input = validate_input.mix(short_forked.for_validate) }
    VALIDATE_READS(validate_input)

    // Collect all validation results and fail if any are invalid
    CHECK_VALIDATION(VALIDATE_READS.out.validated.map { id, reads, status -> status }.collect())

    // Run QC only for needed read types
    if (need_long) {
        NANOPLOT(long_forked.for_nanoplot)
        PARSE_NANOPLOT(NANOPLOT.out.stats)
        FASTPLONG(long_forked.for_fastplong)
        PARSE_FASTPLONG(FASTPLONG.out.json)
        if (!params.skip_rasusa) {
            RASUSA(FASTPLONG.out.reads)
        }
    }

    if (need_short) {
        FASTP_SHORT(short_forked.for_fastp)
        PARSE_FASTP(FASTP_SHORT.out.json)
    }

    // SeqKit stats: mix only the read types we're processing
    seqkit_input = Channel.empty()
    if (need_long)  { seqkit_input = seqkit_input.mix(long_forked.for_seqkit) }
    if (need_short) { seqkit_input = seqkit_input.mix(short_forked.for_seqkit) }
    SEQKIT_STATS(seqkit_input)
    PARSE_SEQKIT_STATS(SEQKIT_STATS.out.stats)

    if (params.stop_at == 'qc') { return }

    // Build clean reads channels: use filtered + subsampled reads where available
    long_reads_ch = need_long
        ? (params.skip_rasusa ? FASTPLONG.out.reads : RASUSA.out.reads)
        : Channel.empty()
    long_clean = need_long
        ? long_reads_ch.mix(samples_ch.filter { it[1].name == 'NO_FILE' }.map { [it[0], it[1]] })
        : samples_ch.map { [it[0], it[1]] }
    short_clean = need_short
        ? FASTP_SHORT.out.reads.mix(samples_ch.filter { it[2].name == 'NO_FILE' }.map { [it[0], it[2], it[3]] })
        : samples_ch.map { [it[0], it[2], it[3]] }
    clean_reads = long_clean.join(short_clean).map { id, l, r1, r2 -> [id, l, r1, r2] }

    assemblies_ch = Channel.empty()

    if (need_long_assembly) {
        LONG_READ_ASSEMBLY(clean_reads)
        if (run_L) {
            assemblies_ch = assemblies_ch.mix(LONG_READ_ASSEMBLY.out.assembly.map { id, asm -> [id, "L", asm] })
        }
    }
    if (run_S) {
        SHORT_READ_ASSEMBLY(clean_reads)
        assemblies_ch = assemblies_ch.mix(SHORT_READ_ASSEMBLY.out.assembly.map { id, asm -> [id, "S", asm] })
    }
    if (run_HL) {
        HYBRID_LONG_FIRST(LONG_READ_ASSEMBLY.out.flye_assembly, LONG_READ_ASSEMBLY.out.long_reads, clean_reads)
        assemblies_ch = assemblies_ch.mix(HYBRID_LONG_FIRST.out.assembly.map { id, asm -> [id, "HL", asm] })
    }
    if (run_HS) {
        HYBRID_SHORT_FIRST(clean_reads)
        assemblies_ch = assemblies_ch.mix(HYBRID_SHORT_FIRST.out.assembly.map { id, asm -> [id, "HS", asm] })
    }

    if (params.stop_at == 'assembly') { return }

    DNAAPLER(assemblies_ch)
    RENAME_CONTIGS(DNAAPLER.out.assembly)

    // Merge read-based and assembly-input samples
    all_renamed_assemblies = RENAME_CONTIGS.out.assembly
    if (has_assemblies && asm_sample_count > 0) {
        all_renamed_assemblies = all_renamed_assemblies.mix(RENAME_CONTIGS_ASM.out.assembly)
    }

    if (params.stop_at == 'polish') { return }

    ALIGN_AND_STATS(all_renamed_assemblies, clean_reads)
    ANNOTATION(all_renamed_assemblies, taxonomy_tool)

    annotate_input = all_renamed_assemblies
        .join(ANNOTATION.out.mob_recon_dir, by: [0, 1])
        .map { id, asm_type, asm, mob_dir -> [id, asm_type, asm, mob_dir] }
    ANNOTATE_PLASMIDS(annotate_input)

    COLLECT_ASSEMBLIES(ANNOTATE_PLASMIDS.out.assembly.map { id, asm_type, asm -> asm }.collect())
    COLLECT_PLASMIDS(ANNOTATION.out.mob_recon_dir)
    COLLECT_CHROMOSOMES_SUMMARY(ANNOTATION.out.mob_recon_dir)

    // --- Outbreak Detection (conditional) - runs after GENERATE_MULTIQC_CONFIG for chromosomal
    // AMR data ---

    // --- Plasmid HGT Detection (conditional, requires 2+ samples) ---
    if (params.plasmid_hgt && sample_count >= 2) {
        PLASMID_HGT_DETECTION(
            ANNOTATION.out.mob_recon_dir,
            ANNOTATION.out.amrfinder_tsv,
            asm_types.toList()
        )
    } else if (params.plasmid_hgt && sample_count < 2) {
        log.info "  Plasmid HGT skipped: requires 2+ samples (found ${sample_count})"
    }

    // --- Species-Based Typing --- Join assemblies with GTDB-Tk species for routing to typing
    // tools
    gtdbtk_species_for_typing = ANNOTATION.out.gtdbtk_parsed
        .map { f ->
            def lines = f.readLines()
            if (lines.size() > 1) {
                def parts = lines[1].split('\t')
                def sample_full = parts[0]
                def species = parts[7] ?: ''
                def sample_parts = sample_full.split('_')
                def asm_type = sample_parts[-1]
                def sample_id = sample_parts[0..-2].join('_')
                return [sample_id, asm_type, species]
            }
            return null
        }
        .filter { it != null }

    // Use all_renamed_assemblies so assembly-input (_A) samples are typed too.
    assemblies_with_species = all_renamed_assemblies
        .join(gtdbtk_species_for_typing, by: [0, 1])
        // [sample_id, asm_type, assembly, species]

    // short_clean = [sample_id, reads1, reads2]; used by reads-based tools (SeroBA).
    SPECIES_TYPING(assemblies_with_species, short_clean)

    seqkit_mqc   = PARSE_SEQKIT_STATS.out.for_multiqc.collect().ifEmpty(file('NO_SEQKIT'))
    coverage_mqc = ALIGN_AND_STATS.out.mosdepth_for_multiqc.collect().ifEmpty(file('NO_COVERAGE'))
    mlst_mqc     = ANNOTATION.out.mlst_for_multiqc.collect().ifEmpty(file('NO_MLST'))
    // In 'both' mode the sylph comparison files ride along; only the report reads them.
    gtdbtk_mqc   = ANNOTATION.out.gtdbtk_parsed.mix(ANNOTATION.out.taxonomy_comparison).collect().ifEmpty(file('NO_GTDBTK'))
    quast_mqc    = ALIGN_AND_STATS.out.quast_for_multiqc.collect().ifEmpty(file('NO_QUAST'))
    amrfinder_mqc = ANNOTATION.out.amrfinder_for_multiqc.collect().ifEmpty(file('NO_AMRFINDER'))
    resfinder_mqc = ANNOTATION.out.resfinder_for_multiqc.collect().ifEmpty(file('NO_RESFINDER'))
    mob_mqc      = ANNOTATION.out.mob_for_multiqc.collect().ifEmpty(file('NO_MOB'))
    quality_mqc = Channel.empty()
    if (need_long)  { quality_mqc = quality_mqc.mix(PARSE_NANOPLOT.out.for_multiqc) }
    if (need_short) { quality_mqc = quality_mqc.mix(PARSE_FASTP.out.for_multiqc) }
    quality_mqc = quality_mqc.collect().ifEmpty(file('NO_QUALITY'))

    filterstats_mqc = Channel.empty()
    if (need_short) { filterstats_mqc = filterstats_mqc.mix(PARSE_FASTP.out.for_multiqc) }
    if (need_long)  { filterstats_mqc = filterstats_mqc.mix(PARSE_FASTPLONG.out.for_multiqc) }
    filterstats_mqc = filterstats_mqc.collect().ifEmpty(file('NO_FILTERSTATS'))

    // Raw dirs for resistance marker joining
    amrfinder_files = ANNOTATION.out.amrfinder_tsv.map { id, asm_type, f -> f }.collect().ifEmpty(file('NO_AMRFINDER_FILES'))
    pheno_files     = ANNOTATION.out.pheno_table.map { id, asm_type, f -> f }.collect().ifEmpty(file('NO_PHENO_FILES'))
    pointfinder_files = ANNOTATION.out.pointfinder.map { id, asm_type, f -> f }.collect().ifEmpty(file('NO_POINTFINDER_FILES'))
    resfinder_hit_files = ANNOTATION.out.resfinder_hits.map { id, asm_type, f -> f }.collect().ifEmpty(file('NO_RESFINDER_HITS'))
    mob_dirs     = ANNOTATION.out.mob_raw_dir.map { id, asm_type, dir -> dir }.collect().ifEmpty(file('NO_MOB_DIRS'))
    virfinder_dirs = ANNOTATION.out.virfinder_raw_dir.map { id, asm_type, dir -> dir }.collect().ifEmpty(file('NO_VIRFINDER_DIRS'))

    ectyper_dirs   = SPECIES_TYPING.out.ectyper_dirs.map { id, asm_type, dir -> dir }.collect().ifEmpty(file('NO_ECTYPER_DIRS'))
    kleborate_dirs = SPECIES_TYPING.out.kleborate_dirs.map { id, asm_type, preset, dir -> dir }.collect().ifEmpty(file('NO_KLEBORATE_DIRS'))
    emmtyper_files = SPECIES_TYPING.out.emmtyper_results.map { id, asm_type, tsv -> tsv }.collect().ifEmpty(file('NO_EMMTYPER'))
    sistr_dirs        = SPECIES_TYPING.out.sistr_dirs.map { id, asm_type, dir -> dir }.collect().ifEmpty(file('NO_SISTR_DIRS'))
    seroba_dirs       = SPECIES_TYPING.out.seroba_dirs.map { id, asm_type, dir -> dir }.collect().ifEmpty(file('NO_SEROBA_DIRS'))
    meningotype_files = SPECIES_TYPING.out.meningotype_results.map { id, asm_type, tsv -> tsv }.collect().ifEmpty(file('NO_MENINGOTYPE'))
    elgato_dirs       = SPECIES_TYPING.out.elgato_dirs.map { id, asm_type, dir -> dir }.collect().ifEmpty(file('NO_ELGATO_DIRS'))

    // Raw mosdepth dist files for cumulative coverage plot
    mosdepth_dist = ALIGN_AND_STATS.out.mosdepth_dist.map { id, asm_type, dist, summary -> dist }.collect().ifEmpty(file('NO_MOSDEPTH_DIST'))

    GENERATE_MULTIQC_CONFIG(seqkit_mqc, coverage_mqc, mlst_mqc, gtdbtk_mqc, quast_mqc, amrfinder_mqc, resfinder_mqc, mob_mqc, quality_mqc, amrfinder_files, pheno_files, pointfinder_files, resfinder_hit_files, mob_dirs, filterstats_mqc, mosdepth_dist, virfinder_dirs, ectyper_dirs, kleborate_dirs, emmtyper_files, sistr_dirs, seroba_dirs, meningotype_files, elgato_dirs)

    // --- Outbreak Detection (conditional, requires 2+ samples) ---
    if (params.outbreak && sample_count >= 2) {
        OUTBREAK_DETECTION(
            ANNOTATION.out.mob_recon_dir,
            ANNOTATION.out.gtdbtk_parsed,
            metadata_ch,
            ANNOTATION.out.mlst_for_multiqc,
            GENERATE_MULTIQC_CONFIG.out.chromosomal_amr
        )

        // --- Core Genome Analysis (optional, separate subworkflow) ---
        cga_tree_svgs = file('NO_CGA_TREES')
        cga_assembled = file('NO_CGA_RESULTS')

        if (params.core_genome) {
            CORE_GENOME_ANALYSIS(
                OUTBREAK_DETECTION.out.cluster_list,
                OUTBREAK_DETECTION.out.clusters_dir,
                metadata_ch,
                OUTBREAK_DETECTION.out.summary,
                ANNOTATION.out.mlst_for_multiqc,
                GENERATE_MULTIQC_CONFIG.out.chromosomal_amr
            )
            cga_tree_svgs = CORE_GENOME_ANALYSIS.out.cga_tree_svgs
            cga_assembled = CORE_GENOME_ANALYSIS.out.cga_results
        }

        BUILD_OUTBREAK_REPORT(
            OUTBREAK_DETECTION.out.tree_svgs,
            OUTBREAK_DETECTION.out.ska_results,
            OUTBREAK_DETECTION.out.summary,
            file('NO_REF_FILES'),
            OUTBREAK_DETECTION.out.mlst_collected,
            GENERATE_MULTIQC_CONFIG.out.chromosomal_amr,
            metadata_ch,
            cga_tree_svgs,
            cga_assembled
        )
    }

    // Outbreak report (empty file if not enabled or single sample)
    outbreak_html = (params.outbreak && sample_count >= 2)
        ? BUILD_OUTBREAK_REPORT.out.report_yaml
        : Channel.value(file('NO_OUTBREAK_REPORT'))

    outbreak_summary_table = (params.outbreak && sample_count >= 2)
        ? BUILD_OUTBREAK_REPORT.out.summary_table_yaml
        : Channel.value(file('NO_OUTBREAK_SUMMARY_TABLE'))

    // HGT reports: RENDER_HGT_TREES + UMAP_PLASMIDS + BUILD_HGT_REPORT run after
    // GENERATE_MULTIQC_CONFIG
    if (params.plasmid_hgt && sample_count >= 2) {
        // Combine assembled results with metadata
        hgt_with_metadata = PLASMID_HGT_DETECTION.out.assembled_results
            .combine(GENERATE_MULTIQC_CONFIG.out.hgt_metadata)

        // Collect GTDB-Tk parsed files for species annotation on trees
        gtdbtk_collected = ANNOTATION.out.gtdbtk_parsed
            .collect()
            .ifEmpty(file('NO_GTDBTK'))

        // Combine tree files with metadata by asm_type
        trees_with_metadata = PLASMID_HGT_DETECTION.out.cluster_trees
            .combine(GENERATE_MULTIQC_CONFIG.out.hgt_metadata)

        // Render trees with toytree (conda env)
        RENDER_HGT_TREES(trees_with_metadata, gtdbtk_collected)

        // UMAP embedding: combine mash matrix + clusters + metadata by asm_type
        umap_input = PLASMID_HGT_DETECTION.out.mash_matrices
            .join(PLASMID_HGT_DETECTION.out.cluster_assignments, by: 0)
            .combine(GENERATE_MULTIQC_CONFIG.out.hgt_metadata)
        UMAP_PLASMIDS(umap_input)

        // AMR-only Mash distance matrix: same population as the AMR plasmid tree, hierarchically
        // clustered and coloured like the outbreak SNP heatmap.
        amr_matrix_input = PLASMID_HGT_DETECTION.out.mash_matrices
            .join(PLASMID_HGT_DETECTION.out.amr_plasmid_dirs, by: 0)
            .combine(GENERATE_MULTIQC_CONFIG.out.hgt_metadata)
        AMR_PLASMID_MATRIX(amr_matrix_input)

        // Build report: join assembled results + metadata + tree SVGs + UMAP HTML (+ comparison maps
        // if enabled) by asm_type
        if (params.hgt_maps) {
            BUILD_HGT_REPORT(
                hgt_with_metadata
                    .join(RENDER_HGT_TREES.out.svg_dir, by: 0)
                    // remainder: true because maps are only drawn for AMR-bearing plasmid groups.
                    .join(PLASMID_HGT_DETECTION.out.comparison_maps, by: 0, remainder: true)
                    .map { items -> items.collect { it == null ? file('NO_COMPARISON_MAPS') : it } }
                    .join(UMAP_PLASMIDS.out.umap_html, by: 0)
                    .join(AMR_PLASMID_MATRIX.out.matrix_html, by: 0)
            )
        } else {
            BUILD_HGT_REPORT(
                hgt_with_metadata
                    .join(RENDER_HGT_TREES.out.svg_dir, by: 0)
                    .map { it + [file('NO_COMPARISON_MAPS')] }
                    .join(UMAP_PLASMIDS.out.umap_html, by: 0)
                    .join(AMR_PLASMID_MATRIX.out.matrix_html, by: 0)
            )
        }
        hgt_reports = BUILD_HGT_REPORT.out.table_yaml
            .mix(BUILD_HGT_REPORT.out.plots_yaml)
            .collect()
            .ifEmpty(file('NO_HGT_REPORT'))

        hgt_summary_table = BUILD_HGT_REPORT.out.summary_yaml
            .collect()
            .ifEmpty(file('NO_HGT_SUMMARY'))
    } else {
        hgt_reports = Channel.value(file('NO_HGT_REPORT'))
        hgt_summary_table = Channel.value(file('NO_HGT_SUMMARY'))
    }

    qc_files = SEQKIT_STATS.out.stats.map { id, rt, f -> f }
        .mix(ALIGN_AND_STATS.out.quast)
        .mix(ALIGN_AND_STATS.out.busco)
        .collect()

    MULTIQC(
        qc_files,
        GENERATE_MULTIQC_CONFIG.out.config,
        GENERATE_MULTIQC_CONFIG.out.general_stats,
        GENERATE_MULTIQC_CONFIG.out.mlst_table,
        GENERATE_MULTIQC_CONFIG.out.taxonomy_table,
        GENERATE_MULTIQC_CONFIG.out.point_mutation_profile.ifEmpty(file("NO_POINT_MUTATION_PROFILE")),
        GENERATE_MULTIQC_CONFIG.out.pointfinder_profile.ifEmpty(file("NO_POINTFINDER_PROFILE")),
        GENERATE_MULTIQC_CONFIG.out.plasmid_table,
        GENERATE_MULTIQC_CONFIG.out.plasmid_resistance,
        GENERATE_MULTIQC_CONFIG.out.chromosome_resistance,
        GENERATE_MULTIQC_CONFIG.out.resistance_profile,
        GENERATE_MULTIQC_CONFIG.out.chromo_gene_profile.ifEmpty(file("NO_CHROMO_GENE_PROFILE")),
        GENERATE_MULTIQC_CONFIG.out.plasmid_gene_profile.ifEmpty(file("NO_PLASMID_GENE_PROFILE")),
        GENERATE_MULTIQC_CONFIG.out.filtered_reads_plot,
        GENERATE_MULTIQC_CONFIG.out.coverage_plot,
        GENERATE_MULTIQC_CONFIG.out.custom_css,
        outbreak_html,
        outbreak_summary_table,
        hgt_reports,
        hgt_summary_table,
        GENERATE_MULTIQC_CONFIG.out.ectyper_mqc,
        GENERATE_MULTIQC_CONFIG.out.kleborate_kpsc_mqc,
        GENERATE_MULTIQC_CONFIG.out.kleborate_kosc_mqc,
        GENERATE_MULTIQC_CONFIG.out.emmtyper_mqc,
        GENERATE_MULTIQC_CONFIG.out.sistr_mqc,
        GENERATE_MULTIQC_CONFIG.out.seroba_mqc,
        GENERATE_MULTIQC_CONFIG.out.meningotype_mqc,
        GENERATE_MULTIQC_CONFIG.out.elgato_mqc
    )

    // --- Store results in database (optional) ---
    if (params.surveillance) {
        // Distinct placeholder names per slot.
        pairwise_ch = (params.outbreak && sample_count >= 2)
            ? BUILD_OUTBREAK_REPORT.out.pairwise_tsv
            : file('NO_PAIRWISE_DISTANCES')
        cga_pairwise_ch = (params.outbreak && params.core_genome && sample_count >= 2)
            ? BUILD_OUTBREAK_REPORT.out.cga_pairwise_tsv.ifEmpty(file('NO_CGA_PAIRWISE_DISTANCES'))
            : file('NO_CGA_PAIRWISE_DISTANCES')

        STORE_RESULTS_DB(
            GENERATE_MULTIQC_CONFIG.out.combined_general_stats,
            GENERATE_MULTIQC_CONFIG.out.combined_taxonomy,
            GENERATE_MULTIQC_CONFIG.out.combined_mlst,
            GENERATE_MULTIQC_CONFIG.out.combined_amrfinder,
            GENERATE_MULTIQC_CONFIG.out.combined_resfinder,
            GENERATE_MULTIQC_CONFIG.out.combined_point_mutations,
            GENERATE_MULTIQC_CONFIG.out.combined_pointfinder,
            GENERATE_MULTIQC_CONFIG.out.combined_resfinder_genes,
            GENERATE_MULTIQC_CONFIG.out.combined_plasmid_resistance,
            GENERATE_MULTIQC_CONFIG.out.combined_chromosomal_resistance,
            pairwise_ch,
            cga_pairwise_ch,
            GENERATE_MULTIQC_CONFIG.out.combined_virulence,
            GENERATE_MULTIQC_CONFIG.out.combined_ectyper,
            GENERATE_MULTIQC_CONFIG.out.combined_kleborate,
            GENERATE_MULTIQC_CONFIG.out.combined_emmtyper,
            GENERATE_MULTIQC_CONFIG.out.combined_sistr,
            GENERATE_MULTIQC_CONFIG.out.combined_seroba,
            GENERATE_MULTIQC_CONFIG.out.combined_meningotype,
            GENERATE_MULTIQC_CONFIG.out.combined_elgato
        )
    }

    // --- Surveillance (optional): SKA2 + ReporTree + PostgreSQL ---
    if (params.surveillance) {
        // Get species from GTDB-Tk parsed files (emitted as plain paths).
        gtdbtk_species = ANNOTATION.out.gtdbtk_parsed
            .map { species_file ->
                def fname = species_file.name.toString().replace('_gtdbtk_parsed.tsv', '')
                def parts = fname.split('_')
                def asm_type = parts[-1]
                def sample_id = parts[0..-2].join('_')
                def lines = species_file.text.trim().split('\n')
                def species = 'Unknown'
                if (lines.size() > 1) {
                    // Data row (skip header)
                    def cols = lines[1].split('\t')
                    // species is column 7 (0-indexed), genus is 6
                    species = (cols.size() > 7 && cols[7]) ? cols[7] : ((cols.size() > 6 && cols[6]) ? cols[6] : 'Unknown')
                }
                [sample_id, asm_type, species]
            }

        // Chromosome bins, not whole assemblies: surveillance measures clonal relatedness, so plasmid
        // content must not contribute to SKA distances.
        SURVEILLANCE(
            ANNOTATION.out.mob_recon_dir,
            gtdbtk_species,
            GENERATE_MULTIQC_CONFIG.out.combined_general_stats,
            GENERATE_MULTIQC_CONFIG.out.chromosomal_amr
        )

        // Runs last: gated on MULTIQC and STORE_RESULTS_DB so completed_at is taken once everything
        // else has finished.
        FINALIZE_WORKFLOW_LOG(
            MULTIQC.out.report,
            STORE_RESULTS_DB.out.log
        )

        // Optional Metabase refresh (sync_schema + rescan_values).
        METABASE_RESCAN(
            FINALIZE_WORKFLOW_LOG.out.log,
            SURVEILLANCE.out.load_log.map { species, log -> species }.collect().ifEmpty([])
        )
    }
}
