/*
 * Cumulative surveillance per species
 * Pages for all isolates and for each group.
 */

include { EXTRACT_CHROMOSOME; SKA_BUILD; SKA_MERGE; SKA_DISTANCE; REPORTREE; PREPARE_ISOLATES_TSV; UPDATE_CUMULATIVE_ISOLATES; SURVEILLANCE_TREE; SURVEILLANCE_DISTANCE_MATRIX; LOAD_SURVEILLANCE_DB } from '../modules/surveillance'

// Folder name for a group's views.
def scopeDir(String groupId) {
    def g = (groupId ?: '').trim()
    if (!g) return 'ungrouped'
    def safe = g.replaceAll(/[^A-Za-z0-9._-]/, '_')
    return safe == 'all' ? 'group_all' : safe
}

// Distinct group_id values in a cumulative isolates TSV ('' = no group recorded).
def groupsIn(cumulative) {
    def lines = cumulative.readLines()
    if (lines.size() < 2) return []
    def header = lines[0].split('\t', -1) as List
    def col = header.indexOf('group_id')
    if (col < 0) return ['']
    return lines.drop(1)
        .findAll { line -> line.trim() }
        .collect { line ->
            def cells = line.split('\t', -1)
            cells.size() > col ? cells[col].trim() : ''
        }
        .unique()
        .sort()
}

workflow SURVEILLANCE {
    take:
    mob_recon_dirs      // [sample_id, asm_type, mob_dir]: MOB-recon output dirs
    gtdbtk_species      // [sample_id, asm_type, species]: species assignments
    general_stats_tsv   // general_statistics.tsv path
    chromosomal_amr_tsv // chromosomal_amr.tsv path (Sample, Chromosomal_AMR_Genes)

    main:
    // Chromosome bins only: samples whose bin is missing or implausibly small are dropped here
    // (optional output) and never enter the cumulative sketch.
    EXTRACT_CHROMOSOME(mob_recon_dirs)

    // Join chromosomes with species: every sample (id + asm_type) is a distinct entity
    species_assemblies = EXTRACT_CHROMOSOME.out.chromosome
        .map { id, asm, fasta -> [id + '_' + asm, fasta] }
        .join(
            gtdbtk_species.map { id, asm, species -> [id + '_' + asm, species] }
        )
        .map { sample, fasta, species -> [species, sample, fasta] }

    // Build SKA sketches: sample name includes assembly type suffix (e.g., SAMPLE-01_L)
    SKA_BUILD(species_assemblies)

    // Group sketches by species
    sketches_by_species = SKA_BUILD.out.sketch
        .map { species, sample, skf -> [species, skf] }
        .groupTuple()

    // Find previous merged sketch for each species
    prev_merged = sketches_by_species
        .map { species, sketches ->
            def species_safe = species.replaceAll(/\s+/, '_')
            def prev_path = "${params.surveillance_sketch_dir}/${species_safe}/merged/${species_safe}_merged.skf"
            def prev_file = file(prev_path)
            [species, sketches, prev_file.exists() ? prev_file : file('NO_PREV_MERGED')]
        }

    SKA_MERGE(prev_merged)

    // Compute pairwise distances from merged sketch
    SKA_DISTANCE(SKA_MERGE.out.merged)

    // Find previous partitions for stable nomenclature
    reportree_in = SKA_DISTANCE.out.distances
        .map { species, distances ->
            def species_safe = species.replaceAll(/\s+/, '_')
            def prev_path = "${params.surveillance_sketch_dir}/${species_safe}/merged/reportree/partitions.tsv"
            def prev_file = file(prev_path)
            [species, distances, prev_file.exists() ? prev_file : file('NO_PREV_PARTITIONS')]
        }

    REPORTREE(reportree_in)

    // Prepare isolates metadata TSV per species species_assemblies emits [species, sample, fasta]
    // where sample = id_asmtype (e.g., SAMPLE-01_L)
    isolate_info = species_assemblies
        .map { species, sample, fasta -> [species, sample] }
        .groupTuple()

    PREPARE_ISOLATES_TSV(isolate_info, general_stats_tsv)

    // Build merged sketch path per species (keyed by species name)
    merged_paths = SKA_MERGE.out.merged
        .map { species, merged ->
            def species_safe = species.replaceAll(/\s+/, '_')
            [species, "${params.surveillance_sketch_dir}/${species_safe}/merged/${species_safe}_merged.skf"]
        }

    db_input = PREPARE_ISOLATES_TSV.out.isolates_tsv
        .join(SKA_DISTANCE.out.distances, remainder: true)
        .join(REPORTREE.out.partitions, remainder: true)
        .join(REPORTREE.out.nomenclature_changes, remainder: true)
        .join(merged_paths, remainder: true)
        .map { items ->
            def species    = items[0]
            def isolates   = items[1]
            def distances  = items[2] ?: file('NO_DISTANCES')
            def partitions = items[3] ?: file('NO_PARTITIONS')
            def changes    = items[4] ?: file('NO_CHANGES')
            def sketch_path = items[5] ?: ''
            [species, isolates, distances, partitions, changes, sketch_path]
        }

    LOAD_SURVEILLANCE_DB(db_input)

    // ── Cumulative isolate metadata sidecar ──
    // PREPARE_ISOLATES_TSV covers only this run; the tree spans every isolate ever sketched.
    prev_cumulative = PREPARE_ISOLATES_TSV.out.isolates_tsv
        .map { species, iso ->
            def species_safe = species.replaceAll(/\s+/, '_')
            def prev_path = "${params.surveillance_sketch_dir}/${species_safe}/${species_safe}_isolates_cumulative.tsv"
            def prev_file = file(prev_path)
            [species, iso, prev_file.exists() ? prev_file : file('NO_PREV_ISOLATES')]
        }

    // .first() makes the AMR table a value channel so every species can read it; a plain queue
    // channel with one item would let only one species run.
    UPDATE_CUMULATIVE_ISOLATES(prev_cumulative, chromosomal_amr_tsv.first())

    // ── Views per species and scope ──
    // Scope 'all' covers every isolate of the species; one extra scope per group found in the
    // cumulative TSV covers that group's isolates.
    scopes = UPDATE_CUMULATIVE_ISOLATES.out.cumulative
        .flatMap { species, cumulative ->
            def out = [[species, 'all', '__ALL__']]
            if (params.surveillance_group_views) {
                groupsIn(cumulative).each { g -> out << [species, scopeDir(g), g] }
            }
            out
        }

    // ── Cumulative interactive tree ──
    tree_input = SKA_DISTANCE.out.distances
        .join(UPDATE_CUMULATIVE_ISOLATES.out.cumulative)
        .join(REPORTREE.out.partitions, remainder: true)
        .map { items ->
            [items[0], items[1], items[2], items[3] ?: file('NO_PARTITIONS')]
        }

    // [species, scope_dir, group, distances, cumulative, partitions]
    SURVEILLANCE_TREE(scopes.combine(tree_input, by: 0))

    // ── Cumulative distance matrix ──
    matrix_input = SKA_DISTANCE.out.distances
        .join(UPDATE_CUMULATIVE_ISOLATES.out.cumulative)

    // [species, scope_dir, group, distances, cumulative]
    SURVEILLANCE_DISTANCE_MATRIX(scopes.combine(matrix_input, by: 0))

    emit:
    load_log = LOAD_SURVEILLANCE_DB.out.log
    partitions = LOAD_SURVEILLANCE_DB.out.saved_partitions
}