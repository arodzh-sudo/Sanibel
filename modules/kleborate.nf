process kleborate {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/kleborate" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        path("kleborate_out/*_output.txt")
        tuple val(meta), path("${meta.id}_kleborate_*_output.txt"), emit: result

    script:
    """
    kleborate -a ${assembly} -o kleborate_out -p ${meta.kleborate_preset} --trim_headers
    for f in kleborate_out/*_output.txt; do cp "\$f" "${meta.id}_kleborate_\$(basename "\$f")"; done
    """
}
