process INDEX_VCF {
    tag "$vcf"

    input:
    path vcf

    output:
    path "input.vcf.gz",     emit: vcf
    path "input.vcf.gz.tbi", emit: tbi

    script:
    """
    if [[ "${vcf}" == *.gz ]]; then
        cp -L ${vcf} input.vcf.gz
    else
        bgzip -c ${vcf} > input.vcf.gz
    fi
    bcftools index -t input.vcf.gz
    """
}
