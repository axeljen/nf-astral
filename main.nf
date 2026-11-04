#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INDEX_REFERENCE   } from './modules/index_reference'
include { INDEX_VCF         } from './modules/index_vcf'
include { MAKE_WINDOWS      } from './modules/make_windows'
include { PREPARE_MASK      } from './modules/prepare_mask'
include { MASK_GENOME       } from './modules/mask_genome'
include { CONSENSUS         } from './modules/consensus'
include { ANGSD_CONSENSUS   } from './modules/angsd_consensus'
include { CRAM2BAM          } from './modules/cram2bam'
include { IDENTIFY_ALIGNMENT } from './modules/identify_alignment'
include { VALIDATE_VCF_SAMPLES } from './modules/validate_vcf_samples'
include { EXTRACT_WINDOW_SEQ } from './modules/extract_window_seq'
include { EXTRACT_REF_WINDOW } from './modules/extract_ref_window'
include { BUILD_ALIGNMENT   } from './modules/build_alignment'
include { FILTER_ALIGNMENT  } from './modules/filter_alignment'
include { IQTREE            } from './modules/iqtree'
include { COLLAPSE_TREE     } from './modules/collapse_tree'
include { COMBINE_TREES     } from './modules/combine_trees'
include { ASTRAL            } from './modules/astral'

workflow {
    main:

    if (!params.reference)    error "Please provide --reference <path to reference fasta>"
    if ([params.bamfiles, params.cramfiles, params.vcf].count { it } != 1) {
        error "Please provide exactly one of --bamfiles, --cramfiles, or --vcf"
    }
    if (!params.samplesheet) error "Please provide --samplesheet <CSV with sample,mask columns> for every input mode"
    if (!(params.angsd_mode in ['random', 'majority'])) {
        error "--angsd_mode must be random or majority"
    }
    if (params.vcf && !(params.consensus_mode.toString() in ['iupac', '1', '2'])) {
        error "--consensus_mode must be iupac, 1, or 2 for VCF input"
    }

    ch_reference = Channel.fromPath(params.reference, checkIfExists: true)
    ch_sample_rows = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.mask) error "Samplesheet must contain non-empty sample and mask columns"
            if (!(row.sample ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) {
                error "Invalid samplesheet sample name '${row.sample}': use letters, digits, underscores, dots, or hyphens"
            }
            tuple(row.sample, file(row.mask, checkIfExists: true))
        }
        .toList()
        .map { rows ->
            if (!rows) error "Samplesheet must contain at least one sample"
            def names = rows.collect { it[0] }
            if (names.unique(false).size() != names.size()) error "Samplesheet contains duplicate sample names"
            rows
        }
    ch_samples = ch_sample_rows.flatMap { it }
    ch_sample_names = ch_sample_rows.map { rows -> rows.collect { it[0] } }

    // ---- reference indexing ----
    INDEX_REFERENCE(ch_reference)
    PREPARE_MASK(ch_samples, INDEX_REFERENCE.out.genome_file.first())

    // ---- windows across the genome ----
    MAKE_WINDOWS(INDEX_REFERENCE.out.genome_file)

    ch_windows = MAKE_WINDOWS.out.bed
        .splitCsv(sep: '\t')
        .map { chrom, start, end, window_id -> tuple(chrom, start, end, window_id) }

    // ---- per-sample consensus genome (computed once per sample) ----
    if (params.bamfiles || params.cramfiles) {
        def extension = params.cramfiles ? '.cram' : '.bam'
        ch_alignments = Channel
            .fromPath(params.cramfiles ?: params.bamfiles, checkIfExists: true)
            .splitText()
            .map { line -> line.trim() }
            .filter { line -> line && !line.startsWith('#') }
            .map { alignment_path ->
                if (!alignment_path.toLowerCase().endsWith(extension)) {
                    error "Alignment list entry must end in ${extension}: ${alignment_path}"
                }
                file(alignment_path, checkIfExists: true)
            }
            .toList()
            .map { alignments ->
                if (!alignments) error "Alignment list must contain at least one ${extension} file"
                if (alignments.unique(false).size() != alignments.size()) error "Alignment list contains duplicate paths"
                alignments
            }
            .flatMap { it }

        IDENTIFY_ALIGNMENT(ch_alignments)
        ch_identified = IDENTIFY_ALIGNMENT.out.alignment
            .map { alignment, sample_file -> tuple(sample_file.text.trim(), alignment) }
            .toList()
            .map { alignments -> [alignments] }
            .combine(ch_sample_names.map { samples -> [samples] })
            .flatMap { alignments, samples ->
                def names = alignments.collect { it[0] }
                if (names.unique(false).size() != names.size()) {
                    error "Multiple alignment files have the same SM tag; provide one file per sample"
                }
                def missing = samples - names
                def unexpected = names - samples
                if (missing || unexpected) {
                    error "Alignment SM tags do not match samplesheet: missing=${missing}, unexpected=${unexpected}"
                }
                alignments
            }

        if (params.cramfiles) {
            CRAM2BAM(ch_identified, INDEX_REFERENCE.out.fasta.first(), INDEX_REFERENCE.out.indexes.first())
            ch_bams = CRAM2BAM.out.bam.map { sample, bam, bai -> tuple(sample, bam) }
        } else {
            ch_bams = ch_identified
        }
        ANGSD_CONSENSUS(ch_bams)
        ch_consensus = ANGSD_CONSENSUS.out.consensus
    } else {
        ch_vcf = Channel.fromPath(params.vcf, checkIfExists: true)

        INDEX_VCF(ch_vcf)
        VALIDATE_VCF_SAMPLES(INDEX_VCF.out.vcf, INDEX_VCF.out.tbi, ch_sample_names)
        CONSENSUS(ch_samples.map { sample, mask -> sample },
            INDEX_REFERENCE.out.fasta.first(), INDEX_REFERENCE.out.indexes.first(),
            VALIDATE_VCF_SAMPLES.out.vcf.first(), VALIDATE_VCF_SAMPLES.out.tbi.first())
        ch_consensus = CONSENSUS.out.consensus
    }

    MASK_GENOME(ch_consensus.join(PREPARE_MASK.out.noncallable, failOnDuplicate: true, failOnMismatch: true))

    // ---- slice all sample consensus genomes once per window ----
    ch_consensus_bundle = MASK_GENOME.out.masked
    .toList()
    .map { consensuses ->
        def sorted = consensuses.sort { a, b -> a[0] <=> b[0] }

        tuple(
            sorted.collect { sample, _fa, _fai -> sample },
            sorted.collect { _sample, fa, _fai -> fa }
        )
    }

    ch_window_consensus = ch_windows.combine(ch_consensus_bundle)

    EXTRACT_WINDOW_SEQ(ch_window_consensus)

    ch_window_seqs = EXTRACT_WINDOW_SEQ.out.seq

    if (params.include_reference) {
        EXTRACT_REF_WINDOW(ch_windows, INDEX_REFERENCE.out.fasta.first(), INDEX_REFERENCE.out.indexes.first())
        ch_window_seqs = ch_window_seqs.mix(EXTRACT_REF_WINDOW.out.seq)
    }

    // ---- build one pseudo-alignment per window ----
    ch_alignment_input = ch_window_seqs
        .groupTuple(by: 0)

    BUILD_ALIGNMENT(ch_alignment_input)
    FILTER_ALIGNMENT(BUILD_ALIGNMENT.out.alignment)

    // drop windows that no longer have enough sequences after filtering
    ch_filtered = FILTER_ALIGNMENT.out.alignment
        .filter { window_id, fasta -> fasta.size() > 0 }

    // ---- gene trees ----
    IQTREE(ch_filtered)
    COLLAPSE_TREE(IQTREE.out.tree)

    // ---- species tree ----
    COMBINE_TREES(COLLAPSE_TREE.out.tree.collect())
    ASTRAL(COMBINE_TREES.out.trees)

    publish: 
        species_tree = ASTRAL.out.tree
        locus_trees = COMBINE_TREES.out.trees

    }

    output {
        species_tree {
            path "${params.name}_species_tree.nwk"
        }
        locus_trees {
            path "${params.name}_locus_trees.nwk"
        }
     }
    
