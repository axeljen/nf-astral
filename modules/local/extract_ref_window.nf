process EXTRACT_REF_WINDOW {
    tag "${window_id}:reference"

    input:
    tuple val(chrom), val(start), val(end), val(window_id)
    path reference
    path reference_indexes

    output:
    tuple val(window_id), path("${window_id}.reference.fa"), emit: seq

    script:
    def region = "${chrom}:${start.toInteger() + 1}-${end}"
    """
    samtools faidx ${reference} ${region} > tmp.fa
    { echo ">reference"; grep -v '^>' tmp.fa; } > ${window_id}.reference.fa
    """
}
