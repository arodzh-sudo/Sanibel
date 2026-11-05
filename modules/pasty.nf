process pasty {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/pasty" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("${meta.id}.tsv"), emit: result
        path("${meta.id}.details.tsv")
        path("${meta.id}.blastn.tsv")

    script:
    """
    pasty --input ${assembly} --prefix ${meta.id} --outdir .
    """
}
