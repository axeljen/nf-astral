#!/usr/bin/env python3
"""Filter a multi-FASTA pseudo-alignment on site-wise and sample-wise missingness.

Missing characters are treated as: N, n, -, ? (case-insensitive gap/ambiguous-N).
Filtering order: sites (columns) first, then samples (sequences), matching the
common convention of removing uninformative columns before assessing how much
missing data each remaining sample carries.
"""
import argparse
import sys

MISSING_CHARS = set("Nn-?")


def read_fasta(path):
    records = []
    header = None
    seq_chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            records.append((header, "".join(seq_chunks)))
    return records


def write_fasta(records, path, wrap=80):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="input multi-fasta alignment")
    parser.add_argument("--output", required=True, help="output filtered multi-fasta alignment")
    parser.add_argument("--max-site-missing", type=float, default=0.5,
                         help="max fraction of missing bases allowed per alignment column")
    parser.add_argument("--max-sample-missing", type=float, default=0.5,
                         help="max fraction of missing bases allowed per sequence (after site filtering)")
    parser.add_argument("--min-samples", type=int, default=4,
                         help="minimum number of sequences required to keep the alignment")
    parser.add_argument("--min-length", type=int, default=1,
                         help="minimum alignment length required after site filtering")
    args = parser.parse_args()

    records = read_fasta(args.input)
    if not records:
        open(args.output, "w").close()
        return

    lengths = {len(seq) for _, seq in records}
    if len(lengths) > 1:
        sys.exit(f"ERROR: sequences in {args.input} have unequal lengths: {lengths}")
    aln_len = lengths.pop()

    # site-wise filtering
    keep_cols = []
    n_seqs = len(records)
    for col in range(aln_len):
        n_missing = sum(1 for _, seq in records if seq[col] in MISSING_CHARS)
        if n_missing / n_seqs <= args.max_site_missing:
            keep_cols.append(col)

    trimmed = [
        (header, "".join(seq[c] for c in keep_cols))
        for header, seq in records
    ]

    # sample-wise filtering
    trimmed_len = len(keep_cols)
    kept_records = []
    for header, seq in trimmed:
        if trimmed_len == 0:
            frac_missing = 1.0
        else:
            n_missing = sum(1 for c in seq if c in MISSING_CHARS)
            frac_missing = n_missing / trimmed_len
        if frac_missing <= args.max_sample_missing:
            kept_records.append((header, seq))

    if len(kept_records) < args.min_samples or trimmed_len < args.min_length:
        # not enough data left to build a meaningful tree for this window
        open(args.output, "w").close()
        return

    write_fasta(kept_records, args.output)


if __name__ == "__main__":
    main()
