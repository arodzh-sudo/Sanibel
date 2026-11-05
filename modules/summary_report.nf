process summary_report {
    tag "summary"
    publishDir { "${params.output}" }, mode: 'copy', pattern: '*report.txt'

    input:
        val  sample_ids
        path serotype_files
        path assembly_stats_files
        path read_metrics_files
        path prokka_txt_files
        path mlst_files
        path kraken_reports
        path pmga_files
        path skani_files
        path resolved_species_files
        path blast16s_files
        path amrfinder_files

    output:
        path "sum_report.txt",    emit: summary
        path "amr_report.txt",    emit: amr,        optional: true
        path "*_mqc.tsv",         emit: mqc_tables, optional: true

    script:
    """
    summary_report.py --samples ${sample_ids.join(' ')}
    """
}
