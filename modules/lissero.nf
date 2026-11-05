process lissero {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/lissero" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("${meta.id}_lissero.txt"), emit: result

    script:
    """
    lissero ${assembly} > ${meta.id}_lissero.txt
    """
}
