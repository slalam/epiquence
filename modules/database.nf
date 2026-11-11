/* Results database, workflow log and Metabase refresh */

include { runSettingsArg } from './utils'

process METABASE_RESCAN {
    label 'process_low'
    tag "metabase"
    cache false

    input:
    // Gates only: this must run after every write to the database, so it depends on the workflow
    // log finalizer and on the per-species surveillance loads.
    val(finalize_log)
    val(surveillance_species)

    output:
    path "metabase_rescan.log", emit: log

    script:
    // Triggers sync_schema then rescan_values.
    def key_arg = runSettingsArg(workflow.launchDir)
    """
    metabase_rescan.py \
        --base-url '${params.metabase_base_url ?: ''}' \
        --db-id '${params.metabase_db_id ?: ''}' \
        ${key_arg} \
        --timeout ${params.metabase_timeout_seconds ?: 120} \
        > metabase_rescan.log 2>&1
    """
}

process FINALIZE_WORKFLOW_LOG {
    label 'process_low'
    tag "workflow_log"
    cache false

    input:
    // Both inputs exist purely to place this last in the DAG: MULTIQC runs after STORE_RESULTS_DB,
    // so gating on both means completed_at is taken once everything else has finished.
    path(multiqc_report)
    path(db_store_log)

    output:
    path "finalize_workflow_log.txt", emit: log

    script:
    def db_arg = runSettingsArg(workflow.launchDir)
    """
    finalize_workflow_log.py \
        ${db_arg} \
        --workflow-id "${params.workflow_id ?: ''}" \
        --workflow-start "${workflow.start}" \
        2>&1 | tee finalize_workflow_log.txt
    """
}

process STORE_RESULTS_DB {
    label 'process_low'
    tag "surveillance"
    cache false

    input:
    path(general_statistics)
    path(taxonomy)
    path(mlst)
    path(amrfinder)
    path(resfinder_phenotype)
    path(point_mutations)
    path(pointfinder)
    path(resfinder_genes)
    path(plasmid_resistance)
    path(chromosomal_resistance)
    path(pairwise_distances)
    path(cga_pairwise_distances)
    path(virulence_factors)
    path(ectyper)
    path(kleborate)
    path(emmtyper)
    path(sistr)
    path(seroba)
    path(meningotype)
    path(elgato)

    output:
    path "db_store.log", emit: log

    script:
    def gs_arg  = general_statistics.name != 'NO_FILE' ? "--general-statistics ${general_statistics}" : ""
    def tax_arg = taxonomy.name != 'NO_FILE' ? "--taxonomy ${taxonomy}" : ""
    def mlst_arg = mlst.name != 'NO_FILE' ? "--mlst ${mlst}" : ""
    def amrf_arg = amrfinder.name != 'NO_FILE' ? "--amrfinder ${amrfinder}" : ""
    def rfp_arg = resfinder_phenotype.name != 'NO_FILE' ? "--resfinder-phenotype ${resfinder_phenotype}" : ""
    def pm_arg  = point_mutations.name != 'NO_FILE' ? "--point-mutations ${point_mutations}" : ""
    def pf_arg  = pointfinder.name != 'NO_FILE' ? "--pointfinder ${pointfinder}" : ""
    def rfg_arg = resfinder_genes.name != 'NO_FILE' ? "--resfinder-genes ${resfinder_genes}" : ""
    def pr_arg  = plasmid_resistance.name != 'NO_FILE' ? "--plasmid-resistance ${plasmid_resistance}" : ""
    def cr_arg  = chromosomal_resistance.name != 'NO_FILE' ? "--chromosomal-resistance ${chromosomal_resistance}" : ""
    // These two use their own placeholder names rather than NO_FILE: both can be absent in the
    // same run, and two inputs sharing a staged filename collide.
    def pw_arg  = pairwise_distances.name != 'NO_PAIRWISE_DISTANCES' ? "--pairwise-distances ${pairwise_distances}" : ""
    def cga_arg = cga_pairwise_distances.name != 'NO_CGA_PAIRWISE_DISTANCES' ? "--cga-pairwise-distances ${cga_pairwise_distances}" : ""
    def vf_arg  = virulence_factors.name != 'NO_FILE' ? "--virulence-factors ${virulence_factors}" : ""
    def ec_arg  = ectyper.name != 'NO_FILE' ? "--ectyper ${ectyper}" : ""
    def kleb_arg = kleborate.name != 'NO_FILE' ? "--kleborate ${kleborate}" : ""
    def emm_arg = emmtyper.name != 'NO_FILE' ? "--emmtyper ${emmtyper}" : ""
    def sistr_arg = sistr.name != 'NO_FILE' ? "--sistr ${sistr}" : ""
    def seroba_arg = seroba.name != 'NO_FILE' ? "--seroba ${seroba}" : ""
    def mening_arg = meningotype.name != 'NO_FILE' ? "--meningotype ${meningotype}" : ""
    def elgato_arg = elgato.name != 'NO_FILE' ? "--elgato ${elgato}" : ""
    def db_arg = runSettingsArg(workflow.launchDir)
    """
    store_results_db.py \
        ${db_arg} \
        --workflow-id "${params.workflow_id ?: ''}" \
        --workflow-description "${params.workflow_description ?: ''}" \
        --group-id "${params.group_id ?: ''}" \
        --group-description "${params.group_description ?: ''}" \
        --assembly-types "${params.assembly_type ?: ''}" \
        --nextflow-run-name "${workflow.runName ?: ''}" \
        --workflow-start "${workflow.start}" \
        ${gs_arg} \
        ${tax_arg} \
        ${mlst_arg} \
        ${amrf_arg} \
        ${rfp_arg} \
        ${pm_arg} \
        ${pf_arg} \
        ${rfg_arg} \
        ${pr_arg} \
        ${cr_arg} \
        ${pw_arg} \
        ${cga_arg} \
        ${vf_arg} \
        ${ec_arg} \
        ${kleb_arg} \
        ${emm_arg} \
        ${sistr_arg} \
        ${seroba_arg} \
        ${mening_arg} \
        ${elgato_arg} \
        2>&1 | tee db_store.log
    """
}