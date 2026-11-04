process ASTRAL {

    input:
    path gene_trees

    output:
    path "species_tree.nwk", emit: tree
    path "astral.log",       emit: log

    script:
    def cmd = params.astral_jar ? "java -jar ${params.astral_jar}" : "astral"
    """
    ${cmd} -i ${gene_trees} -o species_tree.nwk 2> astral.log
    """
}
