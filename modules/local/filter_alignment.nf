process FILTER_ALIGNMENT {
    tag "$window_id"

    input:
    tuple val(window_id), path(alignment)

    output:
    tuple val(window_id), path("${window_id}.filtered.fasta"), emit: alignment

    script:
    """
    filter_alignment.py \\
        --input ${alignment} \\
        --output ${window_id}.filtered.fasta \\
        --max-site-missing ${params.max_site_missing} \\
        --max-sample-missing ${params.max_sample_missing} \\
        --min-samples ${params.min_samples} \
        --min-length ${params.min_length}
    """
}
