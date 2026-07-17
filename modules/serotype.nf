process serotype {
    tag "${meta.id}:${tool}"

    input:
        tuple val(meta), val(tool), path(tool_dir)
    output:
        tuple val(meta), path("${meta.id}_serotype.tsv"), emit: out

    script:
    """
    parse_serotype.py ${tool} ${meta.id} > ${meta.id}_serotype.tsv
    """
}
