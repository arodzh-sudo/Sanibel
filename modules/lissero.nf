process lissero {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("lissero"), emit: results

    script:
    """
    mkdir lissero
    lissero ${assembly} > lissero/lissero_output.txt
    """
}
