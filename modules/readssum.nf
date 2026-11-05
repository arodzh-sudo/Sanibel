process readssum {
    tag "${meta.id}"

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("${meta.id}_readMetrics.txt"), emit: out

    script:
    def prefix = meta.id
    """
    seqkit stats -a -T -j ${task.cpus} ${reads[0]} ${reads[1]} > stats.tsv
    awk -F'\\t' -v g=${meta.genome_size} '
        NR == 1 { for (i = 1; i <= NF; i++) col[\$i] = i; next }
        { n += \$col["num_seqs"]; b += \$col["sum_len"]; q += \$col["num_seqs"] * \$col["AvgQual"] }
        END {
            if (!n || !("AvgQual" in col)) exit 1
            print "avgReadLength avgQuality numReads coverage"
            printf "%.2f %.2f %d %.2f\\n", b / n, q / n, n, b / g
        }' stats.tsv > ${prefix}_readMetrics.txt
    """
}
