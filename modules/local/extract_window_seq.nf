process EXTRACT_WINDOW_SEQ {
    tag "$window_id"

    input:
    tuple val(chrom), val(start), val(end), val(window_id), val(samples), path(consensus_fastas)

    output:
    tuple val(window_id), path("${window_id}.samples.fa"), emit: seq

    script:
    def region = "${chrom}:${start.toInteger() + 1}-${end}"
    """
    set -o pipefail
    : > ${window_id}.samples.fa
    sample_list=('${samples.join("' '")}')
    fasta_list=('${consensus_fastas.join("' '")}')
    for index in "\${!sample_list[@]}"; do
        samtools faidx "\${fasta_list[\$index]}" '${region}' \\
            | sed "s/^>.*/>\${sample_list[\$index]}/" >> ${window_id}.samples.fa
    done
    """
}
