process plasmidfinder {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("plasmidfinder"), emit: results

    script:
    """
    mkdir plasmidfinder
    python -m plasmidfinder -i ${reads[0]} ${reads[1]} -o plasmidfinder
    """
}
