process INDEX_VCF {
    tag "$vcf_name"
    label "small"

    input:
    tuple val(vcf_name), path(vcf_files)

    output:
    path "input.vcf.gz",         emit: vcf
    path "input.vcf.gz.{tbi,csi}", emit: index

    script:
    def vcf = vcf_files.find { it.name == vcf_name }
    """
    if [[ '${vcf}' == *.gz ]]; then
        cp -L '${vcf}' input.vcf.gz
        if [[ -f '${vcf}.tbi' ]]; then
            cp -L '${vcf}.tbi' input.vcf.gz.tbi
        elif [[ -f '${vcf}.csi' ]]; then
            cp -L '${vcf}.csi' input.vcf.gz.csi
        else
            bcftools index -t input.vcf.gz
        fi
    else
        bgzip -c '${vcf}' > input.vcf.gz
        bcftools index -t input.vcf.gz
    fi
    """
}
