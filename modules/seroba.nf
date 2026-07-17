process seroba {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("seroba"), emit: results

    script:
    """
    ln -s ${reads[0]} ${meta.id}_1.fastq.gz
    ln -s ${reads[1]} ${meta.id}_2.fastq.gz
    seroba runSerotyping /seroba/database/ ${meta.id}_1.fastq.gz ${meta.id}_2.fastq.gz ${meta.id}
    mkdir seroba
    mv ${meta.id} seroba/
    """
}
