process parse_assembly {
    tag "${meta.id}"

    input:
        tuple val(meta), path(distances)
    output:
        tuple val(meta), path("${meta.id}_mash_tophit.tsv"), emit: tophit

    script:
    """
    parse_assembly.py ${distances} > ${meta.id}_mash_tophit.tsv
    """
}
