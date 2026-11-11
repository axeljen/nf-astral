process EXTRACT_WINDOW_SEQ {
    tag "$window_id"

    input:
    tuple val(chrom), val(start), val(end), val(window_id), val(samples), path(consensus_fastas)
    path reference
    path reference_indexes

    output:
    tuple val(window_id), path("${window_id}.filtered.fasta"), emit: alignment, optional: true
    tuple val(window_id), path("${window_id}.tree.nwk"), emit: tree, optional: true
    path "${window_id}.filter_summary.txt", emit: filter_summary

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
        --min-length ${params.min_length} \\
        --stats-output ${window_id}.filter_stats.txt
    
    # Check if alignment passed filtering
    if [ -s ${window_id}.filtered.fasta ]; then
        # Alignment passed filtering - run iqtree
        iqtree -s ${window_id}.filtered.fasta -st DNA -m ${params.iqtree_model} \\
            -B ${params.bootstrap} -bnni -nt AUTO -pre ${window_id} -quiet
        nw_ed ${window_id}.treefile 'i & b<=${params.min_support}' o > ${window_id}.tree.nwk
        PASS="PASS"
    else
        # Alignment did not pass filtering
        touch ${window_id}.tree.nwk
        PASS="FAIL"
    fi
    
    # Create summary file
    if [ -f ${window_id}.filter_stats.txt ]; then
        read -r N_SEQS_ORIG N_SITES_ORIG N_SEQS_FILTERED N_SITES_FILTERED < ${window_id}.filter_stats.txt
    else
        N_SEQS_ORIG=0
        N_SITES_ORIG=0
        N_SEQS_FILTERED=0
        N_SITES_FILTERED=0
    fi
    
    echo -e "${window_id}\t\${N_SEQS_ORIG}\t\${N_SITES_ORIG}\t\${N_SEQS_FILTERED}\t\${N_SITES_FILTERED}\t\${PASS}" > ${window_id}.filter_summary.txt
    """
}
