process IQTREE {
    tag "$window_id"

    input:
    tuple val(window_id), path(alignment)

    output:
    tuple val(window_id), path("${window_id}.treefile"), emit: tree

    script:
    """
    iqtree -s ${alignment} -st DNA -m ${params.iqtree_model} \\
        -B ${params.bootstrap} -bnni -nt AUTO -pre ${window_id} -quiet
    """
}
