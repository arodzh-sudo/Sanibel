process seqsero2 {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("seqsero2"), emit: results

    script:
    """
    SeqSero2_package.py -p ${task.cpus} -t 2 -i ${reads[0]} ${reads[1]}
    mkdir seqsero2
    mv SeqSero_result_* seqsero2/
    """
}
