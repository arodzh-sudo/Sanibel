process emm_typing {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'
    containerOptions "--bind ${task.workDir}:/EMBOSS-6.6.0/emboss/.libs"

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("emm_typing"), emit: results

    script:
    """
    mkdir -p emm_typing/groupAstrep_output
    emm_typing.py --fastq_1 ${reads[0]} --fastq_2 ${reads[1]} \\
        -m /db/ -o emm_typing/groupAstrep_output > emm_typing/groupAstrep_result.txt
    """
}
