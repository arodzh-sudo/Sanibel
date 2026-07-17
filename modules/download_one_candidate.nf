process download_one_candidate {
    tag "${meta.id}:${species}"

    input:
        tuple val(meta), val(species), val(accession)
    output:
        tuple val(meta), path("*.fna"), optional: true, emit: fna

    script:
    def n_refs = 5
    def safe   = species.replaceAll(' ', '_')
    """
    download_candidate.sh "${species}" "${accession}" "${safe}" ${n_refs}
    """
}
