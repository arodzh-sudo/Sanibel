process pasty {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("pasty"), emit: results

    script:
    """
    mkdir pasty
    pasty --input ${assembly} --prefix ${meta.id} --outdir pasty
    """
}
