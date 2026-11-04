process INDEX_REFERENCE {
    tag "$reference"

    input:
    path reference

    output:
    path "${reference}",        emit: fasta
    path "${reference}.{fai,gzi}", emit: indexes
    path "genome.txt",          emit: genome_file

    script:
    // BGZF references need both indexes staged together to prevent downstream rebuilding.
    """

    samtools faidx ${reference}
    cut -f1,2 ${reference}.fai > genome.txt
    """
}
