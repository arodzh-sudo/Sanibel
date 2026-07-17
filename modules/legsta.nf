process legsta {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("legsta"), emit: results

    script:
    """
    mkdir legsta
    legsta ${assembly} > legsta/legsta_output.txt
    """
}
