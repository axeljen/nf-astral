process PREPARE_MASK {
    tag "$sample"
    label "small"

    input:
    tuple val(sample), path(mask_bed)
    path genome_file

    output:
    tuple val(sample), path("${sample}.noncallable.bed"), emit: noncallable

    script:
    // complement of the positive (callable) mask = regions to hardmask as N
    """
    bedtools sort -i ${mask_bed} -g ${genome_file} > ${sample}.sorted.bed
    bedtools complement -i ${sample}.sorted.bed -g ${genome_file} > ${sample}.noncallable.bed
    """
}
