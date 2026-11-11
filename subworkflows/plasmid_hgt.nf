/* Plasmid HGT detection */

include { COLLECT_PLASMID_FASTAS; MASH_DIST_PLASMIDS; CLUSTER_PLASMIDS; UMAP_PLASMIDS; SELECT_AMR_PLASMIDS; MASHTREE_AMR_PLASMIDS; PLASMID_ANI_AF; SELECT_PLASMID_REF; BAKTA_PLASMID_REF; COMPARISON_MAP } from '../modules/plasmid_hgt'

workflow PLASMID_HGT_DETECTION {
    take:
    mob_recon_dirs    // channel: [sample_id, asm_type, mob_dir]
    amrfinder_tsvs    // channel: [sample_id, asm_type, amrfinder_tsv]
    asm_types_list    // list of assembly types to process

    main:
    mob_dirs_collected = mob_recon_dirs
        .map { id, asm_type, dir -> dir }
        .collect()

    amrfinder_collected = amrfinder_tsvs
        .map { id, asm_type, f -> f }
        .collect()
        .ifEmpty(file('NO_AMRFINDER_FILES'))

    asm_types_ch = Channel.from(asm_types_list)

    COLLECT_PLASMID_FASTAS(asm_types_ch, mob_dirs_collected, amrfinder_collected)
    MASH_DIST_PLASMIDS(COLLECT_PLASMID_FASTAS.out.plasmid_dir)
    CLUSTER_PLASMIDS(MASH_DIST_PLASMIDS.out.mash_matrix)

    // Parse clusters and prepare per-cluster channels
    cluster_ch = CLUSTER_PLASMIDS.out.clusters
        .combine(MASH_DIST_PLASMIDS.out.plasmid_dir_out, by: 0)
        .flatMap { asm_type, clusters_file, plasmid_dir ->
            def rows = []
            clusters_file.text.eachLine { line, idx ->
                if (idx == 0) return
                def parts = line.split('\t')
                if (parts.size() >= 2) {
                    rows << [asm_type, parts[0], parts[1], plasmid_dir]
                }
            }
            def grouped = rows.groupBy { it[1] }
            grouped.collect { cid, members ->
                [members[0][0], cid, members.collect { it[2] }, members[0][3]]
            }
        }
        .filter { it[2].size() >= 2 }

    cluster_fastas = cluster_ch.map { asm_type, cluster_id, plasmid_names, plasmid_dir ->
        def fastas = plasmid_names.collect { name ->
            file("${plasmid_dir}/${name}.fasta")
        }
        [asm_type, cluster_id, fastas]
    }

    // For COMPARISON_MAP: also stage this assembly type's AMRFinderPlus hits
    cluster_fastas_with_amr = cluster_ch.map { asm_type, cluster_id, plasmid_names, plasmid_dir ->
        def fastas = plasmid_names.collect { name ->
            file("${plasmid_dir}/${name}.fasta")
        }
        def amrfinder_files = []
        plasmid_dir.eachFileMatch(~/.+_amrfinder\.tsv/) { f ->
            amrfinder_files << f
        }
        [asm_type, cluster_id, fastas, amrfinder_files]
    }

    // Run ANI/AF and mashtree in parallel per cluster
    PLASMID_ANI_AF(cluster_fastas)
    // One Mash tree per assembly type across ALL AMR-carrying plasmids, independent of HGT cluster
    // membership (singletons included).
    SELECT_AMR_PLASMIDS(COLLECT_PLASMID_FASTAS.out.plasmid_dir)
    MASHTREE_AMR_PLASMIDS(SELECT_AMR_PLASMIDS.out.selected)

    // Comparison maps only when --hgt_maps is enabled
    if (params.hgt_maps) {
        SELECT_PLASMID_REF(cluster_fastas_with_amr)

        // Flatten: read manifest JSON, emit one tuple per group.
        per_group_ch = SELECT_PLASMID_REF.out.selected
            .flatMap { asm_type, cluster_id, ref_dir, plasmid_fastas, amrfinder_files ->
                def manifest = new groovy.json.JsonSlurper().parse(file("${ref_dir}/groups_manifest.json"))
                def wanted = manifest.groups.findAll { it.has_amr }
                if (!wanted) {
                    log.info "No AMR-bearing plasmid group in ${asm_type}_${cluster_id}; no comparison map drawn"
                }
                wanted.collect { grp ->
                    def suffix = grp.suffix ? "_${grp.suffix}" : ''
                    def ref_fasta = file("${ref_dir}/${grp.ref_fasta}")
                    def amrfinder_tsv = grp.amrfinder_tsv ?: ''
                    [asm_type, cluster_id, suffix, ref_fasta, grp.display_name, amrfinder_tsv, plasmid_fastas, amrfinder_files]
                }
            }

        BAKTA_PLASMID_REF(per_group_ch)
        COMPARISON_MAP(BAKTA_PLASMID_REF.out.annotated)
        maps_by_asm = COMPARISON_MAP.out.map_html_flat
            .groupTuple(by: 0)
    } else {
        maps_by_asm = Channel.empty()
    }

    // Assemble ANI/AF results per assembly type
    ASSEMBLE_HGT_RESULTS(
        PLASMID_ANI_AF.out.ani_af.groupTuple(by: 0)
    )

    // One Newick per assembly type (already keyed by asm_type)
    trees_by_asm = MASHTREE_AMR_PLASMIDS.out.tree

    // AMR-only plasmid set, keyed by asm_type.
    amr_plasmid_dirs = SELECT_AMR_PLASMIDS.out.selected
        .map { asm_type, amr_dir, fallback_nwk -> [asm_type, amr_dir] }

    emit:
    assembled_results = ASSEMBLE_HGT_RESULTS.out.results  // tuple(asm_type, results_dir)
    cluster_trees = trees_by_asm                            // tuple(asm_type, amr_plasmids.nwk)
    comparison_maps = maps_by_asm                           // tuple(asm_type, [html_files])
    mash_matrices = MASH_DIST_PLASMIDS.out.mash_matrix      // tuple(asm_type, mash_dist_matrix.tsv)
    cluster_assignments = CLUSTER_PLASMIDS.out.clusters     // tuple(asm_type, plasmid_clusters.tsv)
    amr_plasmid_dirs = amr_plasmid_dirs                     // tuple(asm_type, amr_plasmids/)
}

process ASSEMBLE_HGT_RESULTS {
    label 'process_low'
    tag "${asm_type}"

    input:
    tuple val(asm_type), val(cluster_ids), path(ani_files)

    output:
    tuple val(asm_type), path("hgt_results_${asm_type}"), emit: results

    script:
    """
    mkdir -p hgt_results_${asm_type}

    for f in ${ani_files}; do
        if [ -f "\$f" ]; then
            cluster=\$(basename "\$f" _ani_af.tsv)
            mkdir -p "hgt_results_${asm_type}/\$cluster"
            cp "\$f" "hgt_results_${asm_type}/\$cluster/"
        fi
    done

    echo "Assembled HGT results for ${asm_type}:"
    find hgt_results_${asm_type} -type f | sort
    """
}