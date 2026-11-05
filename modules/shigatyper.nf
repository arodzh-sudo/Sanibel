process shigatyper {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/shigatyper" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("${meta.id}_shigatyper.txt"), emit: result

    script:
    """
    shigatyper --R1 ${reads[0]} --R2 ${reads[1]} > ${meta.id}_shigatyper.txt
    """
}
