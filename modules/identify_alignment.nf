process EXTRACT_SAMPLE_INFO {
    tag "$alignment"
    label "small"

    input:
    path alignment

    output:
    tuple path(alignment), path("sample.txt"), emit: alignment

    script:
    """
    samtools view -H '${alignment}' > header.sam
    python3 '${projectDir}/bin/validate_sample_names.py' alignment header.sam > sample.txt
    """
}
