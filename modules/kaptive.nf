process kaptive {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
        val variant
    output:
        tuple val(meta), path("kaptive_${variant}"), emit: results

    script:
    def cfg = [
        ab: [ k: [ db: 'ab_k', out: "${meta.id}_ab_k.txt"  ],
              o: [ db: 'ab_o', out: "${meta.id}_ab_oc.txt" ] ],
        vp: [ k: [ db: '/kaptive/reference_database/VibrioPara_Kaptivedb_K.gbk', out: "${meta.id}_vp_k.txt" ],
              o: [ db: '/kaptive/reference_database/VibrioPara_Kaptivedb_O.gbk', out: "${meta.id}_vp_o.txt" ] ],
    ][variant]
    """
    mkdir kaptive_${variant}
    kaptive assembly ${cfg.k.db} ${assembly} -o kaptive_${variant}/${cfg.k.out}
    kaptive assembly ${cfg.o.db} ${assembly} -o kaptive_${variant}/${cfg.o.out}
    """
}
