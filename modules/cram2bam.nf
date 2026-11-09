
process CRAM2BAM {
    tag "$sample"

    input:
    tuple val(sample), path(cram)
    path reference
    path reference_indexes

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam

    script:
    """
    samtools view -T '${reference}' -b -o ${sample}.bam '${cram}'
    samtools index ${sample}.bam
    """
}