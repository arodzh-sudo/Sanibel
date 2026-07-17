process kleborate {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("kleborate"), emit: results

    script:
    """
    mkdir kleborate
    kleborate -a ${assembly} -o kleborate/kleborate_out -p ${meta.kleborate_preset} --trim_headers
    """
}
