include { taxonomyTool } from './annotation'

process ANNOTATION_SUMMARY {
    label 'process_low'
    tag "annotation_summary"
    publishDir "${params.outdir}/summaries/raw", mode: 'copy'

    input:
    path('mlst/*')
    path('amrfinder/*')
    path('pheno/*')
    path('mob/*')
    path('gtdbtk/*')

    output:
    path "mlst.tsv", emit: mlst_combined
    path "amrfinder.tsv", emit: amrfinder_combined, optional: true
    path "resfinder_phenotype.tsv", emit: resfinder_combined, optional: true
    path "taxonomy.tsv", emit: taxonomy_combined, optional: true
    path "plasmids.tsv", emit: plasmids_combined, optional: true

    script:
    """
    echo -e "File\\tScheme\\tST\\tAllele1\\tAllele2\\tAllele3\\tAllele4\\tAllele5\\tAllele6\\tAllele7" > mlst.tsv
    for f in mlst/*_mlst.tsv; do
        [ -f "\$f" ] && cat "\$f" >> mlst.tsv
    done

    # AMRFinderPlus already carries a Name column, so a plain header-once concatenation keeps
    # every row attributable to its sample.
    first=true
    for f in amrfinder/*_amrfinder.tsv; do
        if [ -f "\$f" ]; then
            if [ "\$first" = true ]; then cat "\$f" > amrfinder.tsv; first=false
            else tail -n +2 "\$f" >> amrfinder.tsv; fi
        fi
    done

    # pheno_table.txt has no sample column and a commented preamble, so the sample name is taken
    # from the filename and prepended to every data row.
    first=true
    for f in pheno/*_pheno_table.txt; do
        if [ -f "\$f" ]; then
            sample=\$(basename "\$f" _pheno_table.txt)
            if [ "\$first" = true ]; then
                echo -e "Sample\\tAntimicrobial\\tClass\\tWGS-predicted phenotype\\tMatch\\tGenetic background" > resfinder_phenotype.tsv
                first=false
            fi
            grep -v '^#' "\$f" | grep -v '^[[:space:]]*\$' | sed "s/^/\${sample}\\t/" >> resfinder_phenotype.tsv
        fi
    done

    first=true
    for f in mob/*_contig_report.txt; do
        if [ -f "\$f" ]; then
            if [ "\$first" = true ]; then cat "\$f" > plasmids.tsv; first=false
            else tail -n +2 "\$f" >> plasmids.tsv; fi
        fi
    done

    first=true
    for f in gtdbtk/*_gtdbtk_summary.tsv; do
        if [ -f "\$f" ]; then
            if [ "\$first" = true ]; then cat "\$f" > taxonomy.tsv; first=false
            else tail -n +2 "\$f" >> taxonomy.tsv; fi
        fi
    done
    """
}

process COLLECT_ASSEMBLIES {
    label 'process_low'
    tag "collect_assemblies"
    publishDir "${params.outdir}/summaries/assemblies", mode: 'copy', saveAs: { fn -> fn.split('/')[-1] }

    input:
    path(assemblies)

    output:
    path "out/*.fasta"

    script:
    """
    mkdir -p out
    for f in ${assemblies}; do
        [ -f "\$f" ] && cp "\$f" out/
    done
    """
}

process COLLECT_CHROMOSOMES_SUMMARY {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    // Sits alongside summaries/assemblies and summaries/plasmids.
    publishDir "${params.outdir}/summaries/chromosomes", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(mob_dir)

    output:
    path "${sample_id}_${asm_type}.fasta", optional: true

    script:
    """
    if [ -s "${mob_dir}/chromosome.fasta" ]; then
        cp "${mob_dir}/chromosome.fasta" "${sample_id}_${asm_type}.fasta"
    else
        echo "WARN: no chromosome.fasta in ${mob_dir}; nothing to collect for ${sample_id}_${asm_type}" >&2
    fi
    """
}

process COLLECT_PLASMIDS {
    label 'process_low'
    tag "${sample_id}_${asm_type}"
    publishDir "${params.outdir}/summaries/plasmids", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(mob_dir)

    output:
    path "${asm_type}/*.fasta", optional: true

    script:
    """
    #!/usr/bin/env python3
    import csv, os, shutil

    asm_type = "${asm_type}"
    sample_id = "${sample_id}"
    mob_dir = "${mob_dir}"
    os.makedirs(asm_type, exist_ok=True)

    contig_report = os.path.join(mob_dir, "contig_report.txt")
    if not os.path.exists(contig_report):
        exit(0)

    # Map cluster_id -> list of contig_ids
    cluster_contigs = {}
    with open(contig_report) as f:
        for row in csv.DictReader(f, delimiter='\\t'):
            if row.get('molecule_type', '') == 'plasmid':
                cid = row.get('primary_cluster_id', '-').strip()
                if cid and cid != '-':
                    ctg = row.get('contig_id', '').strip().split()[0]
                    cluster_contigs.setdefault(cid, []).append(ctg)

    # For each plasmid fasta, find which cluster it belongs to
    for pf in sorted(os.listdir(mob_dir)):
        if not pf.startswith('plasmid_') or not pf.endswith('.fasta'):
            continue
        pf_path = os.path.join(mob_dir, pf)
        with open(pf_path) as fh:
            headers = [l.strip().lstrip('>').split()[0] for l in fh if l.startswith('>')]

        for cluster_id, ctgs in cluster_contigs.items():
            if any(h in ctgs for h in headers):
                dst = os.path.join(asm_type, f"p{cluster_id}_{sample_id}_{asm_type}.fasta")
                shutil.copy2(pf_path, dst)
                break
    """
}

process GENERATE_MULTIQC_CONFIG {
    label 'process_low'
    tag "multiqc_config"
    cache false
    publishDir "${params.outdir}/multiqc", mode: 'copy', pattern: '*.yaml'
    publishDir "${params.outdir}/multiqc", mode: 'copy', pattern: '*.css'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'software_versions.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'general_statistics.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'taxonomy.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'mlst.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'amrfinder.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'resfinder_phenotype.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'point_mutations.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'pointfinder.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'resfinder_genes.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'plasmid_resistance_markers.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'chromosomal_resistance_markers.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'virulence_factors.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'ectyper.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'kleborate.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'emmtyper.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'sistr.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'seroba.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'meningotype.tsv'
    publishDir "${params.outdir}/summaries", mode: 'copy', pattern: 'elgato.tsv'

    input:
    path('seqkit/*')
    path('coverage/*')
    path('mlst/*')
    path('gtdbtk/*')
    path('quast/*')
    path('amrfinder/*')
    path('resfinder/*')
    path('mob/*')
    path('quality/*')
    path('amrfinder_files/*')
    path('pheno_files/*')
    path('pointfinder_files/*')
    path('resfinder_hit_files/*')
    path('mob_dirs/*')
    path('filterstats/*')
    path('mosdepth_dist/*')
    path('virfinder_dirs/*')
    path('ectyper_dirs/*')
    path('kleborate_dirs/*')
    path('emmtyper_files/*')
    path('sistr_dirs/*')
    path('seroba_dirs/*')
    path('meningotype_files/*')
    path('elgato_dirs/*')

    output:
    path "general_stats_mqc.yaml", emit: general_stats
    path "mlst_table_mqc.yaml", emit: mlst_table
    path "taxonomy_table_mqc.yaml", emit: taxonomy_table
    path "point_mutation_profile_mqc.yaml", emit: point_mutation_profile, optional: true
    path "pointfinder_profile_mqc.yaml", emit: pointfinder_profile
    path "plasmid_table_mqc.yaml", emit: plasmid_table
    path "plasmid_resistance_mqc.yaml", emit: plasmid_resistance
    path "chromosome_resistance_mqc.yaml", emit: chromosome_resistance
    path "resistance_profile_mqc.yaml", emit: resistance_profile
    path "chromo_gene_profile_mqc.yaml", emit: chromo_gene_profile, optional: true
    path "plasmid_gene_profile_mqc.yaml", emit: plasmid_gene_profile, optional: true
    path "filtered_reads_plot_mqc.yaml", emit: filtered_reads_plot
    path "coverage_plot_mqc.yaml", emit: coverage_plot
    path "multiqc_config.yaml", emit: config
    path "custom_styles.css", emit: custom_css
    path "software_versions.tsv", emit: versions
    path "general_statistics.tsv", emit: combined_general_stats
    path "taxonomy.tsv", emit: combined_taxonomy
    path "mlst.tsv", emit: combined_mlst
    path "amrfinder.tsv", emit: combined_amrfinder
    path "resfinder_phenotype.tsv", emit: combined_resfinder
    path "point_mutations.tsv", emit: combined_point_mutations
    path "pointfinder.tsv", emit: combined_pointfinder
    path "resfinder_genes.tsv", emit: combined_resfinder_genes
    path "plasmid_resistance_markers.tsv", emit: combined_plasmid_resistance
    path "chromosomal_resistance_markers.tsv", emit: combined_chromosomal_resistance
    path "plasmid_hgt_metadata.tsv", emit: hgt_metadata
    path "chromosomal_amr.tsv", emit: chromosomal_amr
    path "virulence_factors.tsv", emit: combined_virulence
    path "ectyper.tsv", emit: combined_ectyper
    path "kleborate.tsv", emit: combined_kleborate
    path "emmtyper.tsv", emit: combined_emmtyper
    path "ectyper_mqc.yaml", emit: ectyper_mqc
    path "kleborate_kpsc_mqc.yaml", emit: kleborate_kpsc_mqc
    path "kleborate_kosc_mqc.yaml", emit: kleborate_kosc_mqc
    path "emmtyper_mqc.yaml", emit: emmtyper_mqc
    path "sistr.tsv", emit: combined_sistr
    path "seroba.tsv", emit: combined_seroba
    path "meningotype.tsv", emit: combined_meningotype
    path "elgato.tsv", emit: combined_elgato
    path "sistr_mqc.yaml", emit: sistr_mqc
    path "seroba_mqc.yaml", emit: seroba_mqc
    path "meningotype_mqc.yaml", emit: meningotype_mqc
    path "elgato_mqc.yaml", emit: elgato_mqc

    script:
    def title = params.multiqc_title ?: "Hybrid Assembly Pipeline Report"
    def subtitle = params.multiqc_subtitle ?: ""
    def intro = params.multiqc_intro ?: "This report summarizes quality control metrics, assembly statistics, and annotation results."
    def analysis_name = params.workflow_description ?: ""
    def workflow_id = params.workflow_id ?: ""
    def workflow_desc = params.workflow_description ?: ""
    def group_id = params.group_id ?: "GP001"
    def group_desc = params.group_description ?: "Healthcare-Associated Infections"
    def tax_tool = taxonomyTool(params)
    """
    generate_multiqc_config.py \
        --title "${title}" \
        --subtitle "${subtitle}" \
        --intro "${intro}" \
        --analysis-name "${analysis_name}" \
        --workflow-id "${workflow_id}" \
        --workflow-description "${workflow_desc}" \
        --group-id "${group_id}" \
        --group-description "${group_desc}" \
        --workflow-timestamp "\$(date '+%Y-%m-%d %H:%M:%S')" \
        --taxonomy-tool ${tax_tool}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    // Central archive of every run's final report, one file per workflow.
    publishDir "${params.analysis_report_dir}", mode: 'copy', overwrite: true,
               pattern: 'multiqc_report.html',
               enabled: params.analysis_report_dir as boolean,
               saveAs: { "${params.workflow_id ?: file(params.outdir).name}_report.html" }
    cache false

    input:
    path('inputs/*')
    path(config)
    path(general_stats)
    path(mlst_table)
    path(taxonomy_table)
    path(point_mutation_profile)
    path(pointfinder_profile)
    path(plasmid_table)
    path(plasmid_resistance)
    path(chromosome_resistance)
    path(resistance_profile)
    path(chromo_gene_profile)
    path(plasmid_gene_profile)
    path(filtered_reads_plot)
    path(coverage_plot)
    path(custom_css)
    path(outbreak_report)
    path(outbreak_summary_table)
    path('hgt_inputs/*')
    path(hgt_summary_table)
    path(ectyper_mqc)
    path(kleborate_kpsc_mqc)
    path(kleborate_kosc_mqc)
    path(emmtyper_mqc)
    path(sistr_mqc)
    path(seroba_mqc)
    path(meningotype_mqc)
    path(elgato_mqc)

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    def ai_enabled = (params.ai_summary || params.ai_summary_full) ? 'true' : 'false'
    """
    mkdir -p inputs

    # Remove any empty or near-empty YAML files that Nextflow may have staged
    find inputs/ -name '*_mqc.yaml' -size -11c -delete 2>/dev/null || true

    # Copy YAML files only if they have meaningful content (>10 bytes)
    for yaml_file in ${general_stats} ${mlst_table} ${taxonomy_table} ${plasmid_table} \
                     ${plasmid_resistance} ${chromosome_resistance} ${resistance_profile} \
                     ${filtered_reads_plot} ${coverage_plot}; do
        if [ -f "\$yaml_file" ] && [ \$(stat -c%s "\$yaml_file") -gt 10 ]; then
            cp "\$yaml_file" inputs/
        fi
    done
    if [ -f "${point_mutation_profile}" ] && [ \$(stat -c%s "${point_mutation_profile}") -gt 10 ]; then cp ${point_mutation_profile} inputs/; fi
    if [ -f "${pointfinder_profile}" ] && [ \$(stat -c%s "${pointfinder_profile}") -gt 10 ]; then cp ${pointfinder_profile} inputs/; fi
    if [ -f "${chromo_gene_profile}" ] && [ \$(stat -c%s "${chromo_gene_profile}") -gt 10 ]; then cp ${chromo_gene_profile} inputs/; fi
    if [ -f "${plasmid_gene_profile}" ] && [ \$(stat -c%s "${plasmid_gene_profile}") -gt 10 ]; then cp ${plasmid_gene_profile} inputs/; fi
    if [ -f "${outbreak_report}" ] && [ "${outbreak_report}" != "NO_OUTBREAK_REPORT" ] && [ \$(stat -c%s "${outbreak_report}") -gt 10 ]; then cp ${outbreak_report} inputs/; fi
    if [ -f "${outbreak_summary_table}" ] && [ "${outbreak_summary_table}" != "NO_OUTBREAK_SUMMARY_TABLE" ] && [ \$(stat -c%s "${outbreak_summary_table}") -gt 10 ]; then cp ${outbreak_summary_table} inputs/; fi
    for hgt in hgt_inputs/*; do
        if [ -f "\$hgt" ] && [ "\$(basename \$hgt)" != "NO_HGT_REPORT" ] && [ \$(stat -c%s "\$hgt") -gt 10 ]; then cp "\$hgt" inputs/; fi
    done
    if [ -f "${hgt_summary_table}" ] && [ "${hgt_summary_table}" != "NO_HGT_SUMMARY" ] && [ \$(stat -c%s "${hgt_summary_table}") -gt 10 ]; then cp ${hgt_summary_table} inputs/; fi
    # Copy species typing YAMLs if they have content
    for st in ${ectyper_mqc} ${kleborate_kpsc_mqc} ${kleborate_kosc_mqc} ${emmtyper_mqc} \
              ${sistr_mqc} ${seroba_mqc} ${meningotype_mqc} ${elgato_mqc}; do
        if [ -f "\$st" ] && [ \$(stat -c%s "\$st") -gt 10 ]; then cp "\$st" inputs/; fi
    done

    # Final safety: remove any remaining empty/tiny YAML files from inputs
    find inputs/ -name '*_mqc.yaml' -size -11c -delete 2>/dev/null || true

    export HOME="\$PWD"
    export TMPDIR="\$PWD/tmp"; mkdir -p "\$TMPDIR"

    if [ "${ai_enabled}" = "true" ]; then
        # Local OpenAI-compatible endpoint (e.g. Ollama); a key is required but ignored.
        export OPENAI_API_KEY="\${OPENAI_API_KEY:-ollama}"
        {
            echo "ai_summary: ${params.ai_summary}"
            echo "ai_summary_full: ${params.ai_summary_full}"
            echo "ai_provider: custom"
            echo "ai_model: ${params.ai_model}"
            echo "ai_custom_endpoint: ${params.ai_endpoint}"
            echo "ai_custom_context_window: ${params.ai_context_window}"
        } > ai_config.yaml

        if multiqc inputs/ --config ${config} --config ai_config.yaml -f -n multiqc_report; then
            echo "MultiQC: report generated with AI summary (${params.ai_model})" >&2
        else
            echo "WARN: MultiQC AI summary failed (is ${params.ai_endpoint} reachable?); regenerating without AI" >&2
            multiqc inputs/ --config ${config} -f -n multiqc_report
        fi

        # Wrap the AI summary in a native <details> so it collapses without JS.
        ai_collapse.py multiqc_report.html || echo "WARN: AI collapse wrap skipped" >&2
    else
        multiqc inputs/ --config ${config} -f -n multiqc_report
    fi

    cat > _fix_theme.py << 'PYEOF'
with open('multiqc_report.html', 'r') as f:
    html = f.read()
script = '''<script>(function(){var h=document.documentElement;h.setAttribute("data-bs-theme","light");new MutationObserver(function(){if(h.getAttribute("data-bs-theme")!=="light"){h.setAttribute("data-bs-theme","light")}}).observe(h,{attributes:true,attributeFilter:["data-bs-theme"]})})()</script>'''
idx = html.rfind('</body>')
if idx != -1:
    html = html[:idx] + script + html[idx:]
with open('multiqc_report.html', 'w') as f:
    f.write(html)
PYEOF
    python3 _fix_theme.py
    rm -f _fix_theme.py
    """
}