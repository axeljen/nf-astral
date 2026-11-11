process SUMMARIZE_SAMPLE_SURVIVAL {
    label "small"

    input:
    val sample_names
    path alignments

    output:
    path "sample_survival_summary.tsv", emit: summary

    script:
    def sample_list = sample_names.join(",")
    def alignment_list = alignments.join(" ")
    """
    python3 << 'EOF'
    import os
    from pathlib import Path

    samples = '${sample_list}'.split(',')
    alignment_files = '${alignment_list}'.split()
    
    sample_retention = {sample: 0 for sample in samples}
    
    for alignment_file in alignment_files:
        with open(alignment_file) as f:
            for line in f:
                if line.startswith('>'):
                    sample_id = line[1:].strip()
                    if sample_id in sample_retention:
                        sample_retention[sample_id] += 1
    
    with open('sample_survival_summary.tsv', 'w') as f:
        f.write('sample\\tnumber_of_trees\\n')
        for sample in samples:
            f.write(f'{sample}\\t{sample_retention[sample]}\\n')
    EOF
    """
}
