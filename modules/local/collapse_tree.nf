process COLLAPSE_TREE {
    tag "$window_id"

    input:
    tuple val(window_id), path(tree)

    output:
    path "${window_id}.collapsed.nwk", emit: tree

    script:
    // collapse internal branches with bootstrap support <= min_support into polytomies
    """
    nw_ed ${tree} 'i & b<=${params.min_support}' o > ${window_id}.collapsed.nwk
    """
}
