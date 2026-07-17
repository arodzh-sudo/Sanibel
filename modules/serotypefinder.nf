process serotypefinder {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("serotypefinder"), emit: results

    script:
    """
    mkdir serotypefinder
    serotypefinder.py -i ${reads[0]} ${reads[1]} -o serotypefinder
    """
}
