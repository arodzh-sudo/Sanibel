process quast {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/assembly" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("${meta.id}_quast_report.tsv"), emit: report

    script:
    """
    quast.py --threads ${task.cpus} -o quast_results ${assembly}
    mv quast_results/report.tsv ${meta.id}_quast_report.tsv
    """
}
