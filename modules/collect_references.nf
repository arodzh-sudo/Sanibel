process collect_references {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/reference_genomes" }, mode: 'copy', pattern: '*_reference_genomes.txt'

    input:
        tuple val(meta), path(fnas)
    output:
        tuple val(meta), path("${meta.id}_references/"),           emit: references
        tuple val(meta), path("${meta.id}_reference_genomes.txt"), emit: manifest

    script:
    """
    mkdir -p ${meta.id}_references
    cp *.fna ${meta.id}_references/ 2>/dev/null || true

    ls ${meta.id}_references/*.fna 2>/dev/null | xargs -r -n1 basename \\
        | sed 's/\\.fna\$//; s/__/_/' | sort > ${meta.id}_reference_genomes.txt

    n=\$(find ${meta.id}_references -name '*.fna' | wc -l)
    [ "\$n" -gt 0 ] || { echo "No reference genomes for ${meta.id}" >&2; exit 1; }
    """
}
