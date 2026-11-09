process ANGSD_CONSENSUS {
    tag "$sample"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.unmasked.fa.gz"), path("${sample}.unmasked.fa.gz.{fai,gzi}"), emit: consensus

    script:
    def dofasta_arg = params.angsd_mode == 'random' ? "1" : "2"
    """
    angsd \
        -i '${bam}' \
        -doFasta ${dofasta_arg} \
        -doCounts 1 \
        -explode 1 \
        -seed ${params.angsd_seed} \
        -out ${sample}.unmasked
    samtools faidx ${sample}.unmasked.fa.gz
    """
}