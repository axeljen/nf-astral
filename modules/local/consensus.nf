process CONSENSUS {
    tag "$sample"

    input:
    val sample
    path reference
    path reference_indexes
    path vcf
    path vcf_tbi

    output:
    tuple val(sample), path("${sample}.unmasked.fa.gz"), path("${sample}.unmasked.fa.gz.fai"), emit: consensus

    script:
    def hap_flag = params.consensus_mode == 'iupac' ? '-I' : "-H ${params.consensus_mode}"
    """
    bcftools query -s ${sample} -i 'GT="mis"' \
        -f '%CHROM\t%POS0\t%END\n' ${vcf} > ${sample}.missing.bed
    bcftools consensus -f ${reference} -m ${sample}.missing.bed \
        -s ${sample} ${hap_flag} ${vcf} | bgzip > ${sample}.unmasked.fa.gz
    samtools faidx ${sample}.unmasked.fa.gz
    """
}
