process BUILD_ALIGNMENT {
    tag "$window_id"

    input:
    tuple val(window_id), path(fastas)

    output:
    tuple val(window_id), path("${window_id}.aln.fasta"), emit: alignment

    script:
    // sequences are all derived from the same reference coordinates via SNP
    // substitution only (no indels), so simple concatenation already yields
    // a positionally-aligned pseudo-alignment
    """
    cat ${fastas} > ${window_id}.aln.fasta
    """
}
