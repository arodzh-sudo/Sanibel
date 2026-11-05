process legsta {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/legsta" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("${meta.id}_legsta.txt"), emit: result

    script:
    """
    legsta ${assembly} > ${meta.id}_legsta.txt
    """
}
