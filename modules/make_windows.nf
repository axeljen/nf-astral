process MAKE_WINDOWS {
    tag "windows"
    label "small"

    input:
    path genome_file

    output:
    path "windows.bed", emit: bed

    script:
    // tile the genome in `window_step`-sized bins, then truncate each bin down
    // to `window_size` so windows are `window_size` long and `window_step` apart
    """
    bedtools makewindows -g ${genome_file} -w ${params.window_step} -s ${params.window_step} \\
        | awk -v OFS='\\t' -v wsize=${params.window_size} -v minlen=${params.min_length} '
            { end = \$2 + wsize; if (end > \$3) end = \$3;
              if (end > \$2 && end - \$2 >= minlen) print \$1, \$2, end, \$1"_"\$2"_"end }' \\
        > windows.bed
    """
}
