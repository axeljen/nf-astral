process COMBINE_TREES {
    publishDir "${params.outdir}/gene_trees", mode: 'copy'
    label "small"

    input:
    path trees

    output:
    path "gene_trees.nwk", emit: trees

    script:
    """
    cat ${trees} > gene_trees.nwk
    """
}
