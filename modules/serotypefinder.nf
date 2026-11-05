process serotypefinder {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/serotypefinder" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        path("data.json")
        tuple val(meta), path("${meta.id}_serotypefinder.json"), emit: result

    script:
    """
    serotypefinder.py -i ${reads[0]} ${reads[1]} -o ./
    cp data.json ${meta.id}_serotypefinder.json
    """
}
