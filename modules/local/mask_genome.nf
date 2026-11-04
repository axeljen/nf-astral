process MASK_GENOME {
    tag "$sample"

    input:
    tuple val(sample), path(consensus), path(consensus_indexes), path(noncallable_bed)

    output:
    tuple val(sample), path("${sample}.consensus.fa.gz"), path("${sample}.consensus.fa.gz.fai"), emit: masked

    script:
    // def mask_input = consensus.name.endsWith('.gz') ? "${sample}.unmasked.fa" : consensus

    """
    bedtools maskfasta -fi <(bgzip -dc '${consensus}') -bed '${noncallable_bed}' -fo ${sample}.consensus.fa
    bgzip ${sample}.consensus.fa
    samtools faidx ${sample}.consensus.fa.gz
    """
}
