process IQTREE {
    tag "$window_id"

    input:
    tuple val(window_id), path(alignment)

    output:
    path "${window_id}.collapsed.nwk", emit: tree

    script:
    """
    iqtree -s ${alignment} -st DNA -m ${params.iqtree_model} \\
        -B ${params.bootstrap} -bnni -nt AUTO -pre ${window_id} -quiet
    nw_ed ${window_id}.treefile 'i & b<=${params.min_support}' o > ${window_id}.collapsed.nwk
    """
}
