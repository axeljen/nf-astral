process SUMMARIZE_FILTERING {
    label "small"

    input:
    path filter_summaries

    output:
    path "alignment_filtering_summary.tsv", emit: summary

    script:
    """
    (
        echo -e "window_id\tseq_count_orig\tsite_count_orig\tseq_count_retained\tsite_count_retained\ttree_generated"
        cat ${filter_summaries} | sort
    ) > alignment_filtering_summary.tsv
    """
}
