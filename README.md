# nf-astral

A Nextflow pipeline that builds a genome-wide, multi-species coalescent based phylogeny directly from mapped bam/cram files or a vcf file with filtered snps, together with sample-specific mappability masks. 

## Get started

Requirements:
- [Nextflow](https://www.nextflow.io/)
- [Singularity/apptainer](https://apptainer.org/)

With these software requirements installed, you can run the pipeline as described below:

```bash
# clone repository
git clone https://github.com/axeljen/nf-astral.git
cd nf-astral
```

Prepare a sample sheet with two columns, sample and mask file (see `assets/samplesheet.csv`) for an example.

To run with a vcf file as source:

```bash
nextflow run main.nf \
   --vcf snps.vcf.gz \
   --samplesheet assets/samplesheet.csv \
   --reference reference.fasta
```

To run with BAM/CRAM files as source:

```bash
nextflow run main.nf \
   --bamfiles list/of/bamfiles.txt \
   --samplesheet assets/samplesheet.csv \
   --reference reference.fasta
```

## Workflow

1. **Windows**: `bedtools makewindows` tiles the genome in `window_step`-sized
   bins and each bin is truncated to `window_size`, giving windows of
   `window_size` bp spaced `window_step` bp apart (defaults: 10 kb windows,
   200 kb apart) to reduce linkage between loci used for gene trees.
   Windows shorter than `--min_length` are excluded from the generated BED,
   before any per-window tasks are submitted.
2. **Per-sample consensus** (computed once per sample, not once per window,
   for efficiency):
   - **CRAM**: `CRAM2BAM` decodes each CRAM with the supplied reference,
     then follows the BAM path. A CRAM index is not required.
   - **BAM**: `ANGSD_CONSENSUS` runs `angsd -doFasta` with random read-base
     sampling by default (`--angsd_mode random`, `-doFasta 1`).
     ANGSD's native BGZF-compressed FASTA is indexed directly with
     `samtools faidx` in the same process, producing `.fai` and `.gzi`
     indexes without recompression. `--angsd_mode majority` selects the most frequent base
     (`-doFasta 2`).
   - **VCF**: `bcftools consensus` applies each requested sample's genotypes
     to the reference, masking missing genotypes to `N`. Non-variant sites
     inside the callable region are implicitly homozygous reference.
   All branches then use the same `MASK_GENOME` process to hardmask the
   complement of each sample's callable BED to `N`. This process decompresses
   BGZF input before running `bedtools maskfasta`. Variants outside the
   callable mask are also masked (unlike the previous mask-before-consensus
   VCF behavior).
3. **Window slicing**: for every window, the corresponding region is sliced
   out of each sample's consensus genome (`samtools faidx`) and, optionally,
   from the raw reference (`--include_reference`). Because all sequences are
   derived from the same reference coordinates via substitution only (no
   indels), concatenating them per window already yields a positionally
   aligned pseudo-alignment. Extraction, optional reference inclusion, and
   missingness filtering run in a single `EXTRACT_WINDOW_SEQ` task per window.
4. **Missingness filtering**: a small bundled script (`bin/filter_alignment.py`)
   drops alignment columns exceeding `--max_site_missing`, then drops
   sequences exceeding `--max_sample_missing` (computed after site
   filtering). Windows left with fewer than `--min_samples` sequences are
   dropped, as are windows shorter than `--min_length` after site filtering.
   Swap this step out for your own script, or for a tool such as
   [trimAl](http://trimal.cgenomics.org/) (`-gt`) / [ClipKit](https://github.com/JLSteenwyk/ClipKit)
   for the site-wise part if preferred.
5. **Gene trees**: `iqtree2 -bb <bootstrap>` (default 1000 ultrafast
   bootstrap replicates) builds a tree per window.
6. **Branch collapsing**: `nw_ed '<tree> "i & b<=<min_support>" o'` (Newick
   Utilities) collapses internal branches with support at or below
   `min_support` (default 20%) into polytomies. Tree inference and collapsing
   run in a single `IQTREE` task per retained window.
7. **Species tree**: all collapsed gene trees are concatenated and passed to
   [ASTRAL](https://github.com/smirarab/ASTRAL) to infer the species tree.

## Usage

```bash
nextflow run main.nf \
    --reference reference.fasta \
    --vcf samples.vcf.gz \
    --samplesheet samplesheet.csv \
    --outdir results
```

For alignment input, replace `--vcf samples.vcf.gz` with
`--bamfiles alignments.txt`. The list may contain BAM and CRAM paths in any
combination; CRAM files are converted to BAM before consensus generation.
Supply **exactly one** of these two options, always with `--samplesheet` and
`--reference`.

Each list contains one file path per line, for example:

```text
/path/to/sampleA.bam
/path/to/sampleB.bam
```

Each path must end in `.bam` or `.cram`. Blank lines and lines beginning with
`#` are ignored. Relative paths in lists and samplesheets are resolved from
the directory where Nextflow is launched. Alignments must be coordinate-sorted;
the reference must match the one used to create the alignments, especially
for CRAM decoding.

### Samplesheet (`--samplesheet`)

CSV with a header row and one line per sample (see
[assets/samplesheet.csv](assets/samplesheet.csv)):

```csv
sample,mask
sampleA,/path/to/sampleA.callable.bed
sampleB,/path/to/sampleB.callable.bed
```

`sample` must match a VCF sample name or the `SM` tag in the BAM/CRAM `@RG`
header records. Every alignment file must contain exactly one unique, non-empty `SM` sample name;
multiple read groups for that same sample are allowed. 

`mask` is a BED file of trusted regions (positive/callable mask), required
and applied in **all** input modes. All branches converge on the same masking,
window extraction, alignment filtering, and tree construction.

### Key parameters

| Parameter               | Default   | Description                                              |
|--------------------------|-----------|-----------------------------------------------------------|
| `--window_size`          | `10000`   | Length (bp) of each window                                 |
| `--window_step`          | `200000`  | Distance (bp) between window starts                        |
| `--consensus_mode`       | `iupac`   | `iupac` for ambiguity codes, or a haplotype number (`1`/`2`) for `bcftools consensus -H` |
| `--angsd_mode`           | `random`  | BAM/CRAM consensus: `random` read base or `majority` base |
| `--angsd_seed`           | `1`       | Random seed for ANGSD consensus |
| `--include_reference`    | `false`   | Add the raw reference sequence to each window alignment    |
| `--max_site_missing`     | `0.5`     | Max fraction missing per alignment column                  |
| `--max_sample_missing`   | `0.5`     | Max fraction missing per sequence (after site filtering)   |
| `--min_samples`          | `4`       | Minimum sequences required to keep/tree a window            |
| `--iqtree_model`         | `MFP`     | Substitution model passed to IQ-TREE. Default is to detect the best-fit model automatically                        |
| `--bootstrap`            | `1000`    | Ultrafast bootstrap replicates (`iqtree -B`)               |
| `--min_support`          | `20`      | Collapse branches with support at/below this value (%)     |
| `--astral_jar`           | `null`    | Path to `astral.jar`; otherwise the `astral` executable on `$PATH` is used |

## Requirements

Singularity is enabled by default, with automatic bind mounts. 

The `EXTRACT_WINDOW_SEQ` image contains samtools and Python; the `IQTREE`
image contains both IQ-TREE and Newick Utilities. Custom container overrides
must provide the same combined toolsets.

## Notes / assumptions

- The reference may be plain FASTA or BGZF-compressed FASTA. An existing
  `.fai` index next to the reference is reused; otherwise the reference is
  indexed. Any existing `.gzi` index is staged with the reference.
- Uncompressed VCFs are bgzipped and indexed automatically. Existing
  `.tbi`/`.csi` indexes next to bgzipped VCFs are reused; otherwise the VCF is
  indexed.
- Only SNPs are expected; the pipeline does not handle indels in the
  consensus step, since `bcftools consensus` on indel-containing VCFs would
  break the fixed-coordinate concatenation used to build the pseudo-alignment.

## Tests

With Nextflow, bedtools, samtools, IQ-TREE (`iqtree`), and Newick Utilities
(`nw_ed`) on `PATH`, run the focused integration tests:

```bash
python3 -m unittest discover -s tests -v
```
