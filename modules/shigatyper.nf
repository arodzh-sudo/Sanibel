process shigatyper {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("shigatyper"), emit: results

    script:
    """
    mkdir shigatyper
    shigatyper --R1 ${reads[0]} --R2 ${reads[1]} > shigatyper/shigatyper_output.txt
    """
}
