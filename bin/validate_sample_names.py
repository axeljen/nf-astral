#!/usr/bin/env python3
"""Validate alignment SM tags or requested VCF sample names."""
import argparse
from pathlib import Path
import sys


def alignment_sample(header):
    samples = set()
    for line in header.splitlines():
        if not line.startswith("@RG\t"):
            continue
        tags = [field[3:] for field in line.split("\t")[1:] if field.startswith("SM:")]
        if len(tags) != 1 or not tags[0]:
            raise ValueError("Every alignment @RG record must have exactly one non-empty SM tag")
        samples.add(tags[0])
    if len(samples) != 1:
        raise ValueError(
            f"Alignment must contain exactly one unique SM sample; found {sorted(samples)}"
        )
    return samples.pop()


def validate_vcf_samples(available, requested):
    missing = sorted(set(requested) - set(available))
    if missing:
        raise ValueError(f"Samplesheet samples missing from VCF: {', '.join(missing)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["alignment", "vcf"])
    parser.add_argument("input", type=Path)
    parser.add_argument("samples", nargs="*")
    args = parser.parse_args()
    try:
        content = args.input.read_text()
        if args.mode == "alignment":
            print(alignment_sample(content))
        else:
            validate_vcf_samples(content.splitlines(), args.samples)
    except (OSError, ValueError) as error:
        sys.exit(f"ERROR: {error}")


if __name__ == "__main__":
    main()
