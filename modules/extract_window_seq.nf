process EXTRACT_WINDOW_SEQ {
    tag "$window_id"

    input:
    tuple val(chrom), val(start), val(end), val(window_id), val(samples), path(consensus_fastas)
    path reference
    path reference_indexes

    output:
    tuple val(window_id), path("${window_id}.filtered.fasta"), emit: alignment

    script:
    def region = "${chrom}:${start.toInteger() + 1}-${end}"
    def reference_command = params.include_reference ? """
    samtools faidx '${reference}' '${region}' \\
        | sed 's/^>.*/>reference/' >> ${window_id}.samples.fa
    """ : ''
    """
    set -o pipefail
    : > ${window_id}.samples.fa
    sample_list=('${samples.join("' '")}')
    fasta_list=('${consensus_fastas.join("' '")}')
    for index in "\${!sample_list[@]}"; do
        samtools faidx "\${fasta_list[\$index]}" '${region}' \\
            | sed "s/^>.*/>\${sample_list[\$index]}/" >> ${window_id}.samples.fa
    done
    ${reference_command}
    python3 '${projectDir}/bin/filter_alignment.py' \\
        --input ${window_id}.samples.fa \\
        --output ${window_id}.filtered.fasta \\
        --max-site-missing ${params.max_site_missing} \\
        --max-sample-missing ${params.max_sample_missing} \\
        --min-samples ${params.min_samples} \\
        --min-length ${params.min_length}
    """
}
