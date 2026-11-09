process INDEX_REFERENCE {
    tag "$reference_name"
    label "small"

    input:
    tuple val(reference_name), path(reference_files)

    output:
    path "${reference_name}",           emit: fasta
    path "${reference_name}.{fai,gzi}", emit: indexes
    path "genome.txt",                  emit: genome_file

    script:
    def reference = reference_files.find { it.name == reference_name }
    """
    if [[ ! -f '${reference}.fai' ]]; then
        samtools faidx '${reference}'
    fi
    cut -f1,2 '${reference}.fai' > genome.txt
    """
}
