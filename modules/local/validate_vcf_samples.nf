process VALIDATE_VCF_SAMPLES {
    tag "VCF sample names"

    input:
    path vcf
    path tbi
    val samples

    output:
    path vcf, emit: vcf
    path tbi, emit: tbi

    script:
    def sample_args = samples.collect { "'${it}'" }.join(' ')
    """
    bcftools query -l '${vcf}' > vcf_samples.txt
    python3 '${projectDir}/bin/validate_sample_names.py' vcf vcf_samples.txt ${sample_args}
    """
}
