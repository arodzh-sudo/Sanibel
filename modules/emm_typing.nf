process emm_typing {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/emm_typing" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("${meta.id}_emm_typing.txt"), emit: result
        path("groupAstrep_output/")

    script:
    """
    mkdir -p groupAstrep_output
    emm_typing.py --fastq_1 ${reads[0]} --fastq_2 ${reads[1]} \\
        -m /db/ -o groupAstrep_output > ${meta.id}_emm_typing.txt
    """
}
