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
include { EXTRACT_SAMPLE_INFO } from './modules/identify_alignment'
include { VALIDATE_VCF_SAMPLES } from './modules/validate_vcf_samples'
include { EXTRACT_WINDOW_SEQ } from './modules/extract_window_seq'
include { COMBINE_TREES     } from './modules/combine_trees'
include { ASTRAL            } from './modules/astral'
include { SUMMARIZE_FILTERING } from './modules/summarize_filtering'
include { SUMMARIZE_SAMPLE_SURVIVAL } from './modules/summarize_sample_retention'

workflow {
    main:

    if (!params.reference)    error "Please provide --reference <path to reference fasta>"
    if ([params.bamfiles, params.vcf].count { it } != 1) {
        error "Please provide exactly one of --bamfiles or --vcf"
    }
    if (!params.samplesheet) error "Please provide --samplesheet <CSV with sample,mask columns> for every input mode"
    if (!(params.angsd_mode in ['random', 'majority'])) {
        error "--angsd_mode must be random or majority"
    }
    if (params.vcf && !(params.consensus_mode.toString() in ['iupac', '1', '2'])) {
        error "--consensus_mode must be iupac, 1, or 2 for VCF input"
    }

    ch_reference = Channel
        .fromPath(params.reference, checkIfExists: true)
        .map { reference ->
            def reference_files = [reference]
            ['.fai', '.gzi'].each { suffix ->
                def index = file("${reference}${suffix}")
                if (index.exists()) reference_files << index
            }
            tuple(reference.name, reference_files)
        }
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
    if (params.bamfiles) {
        ch_alignments = Channel
            .fromPath(params.bamfiles, checkIfExists: true)
            .splitText()
            .map { line -> line.trim() }
            .filter { line -> line && !line.startsWith('#') }
            .map { alignment_path ->
                def extension = alignment_path.toLowerCase()
                if (!(extension.endsWith('.bam') || extension.endsWith('.cram'))) {
                    error "Alignment list entries must end in .bam or .cram: ${alignment_path}"
                }
                file(alignment_path, checkIfExists: true)
            }
            .toList()
            .map { alignments ->
                if (!alignments) error "Alignment list must contain at least one BAM or CRAM file"
                if (alignments.unique(false).size() != alignments.size()) error "Alignment list contains duplicate paths"
                alignments
            }
            .flatMap { it }

        EXTRACT_SAMPLE_INFO(ch_alignments)
        ch_identified = EXTRACT_SAMPLE_INFO.out.alignment
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

        ch_crams = ch_identified.filter { sample, alignment -> alignment.name.toLowerCase().endsWith('.cram') }
        ch_bam_files = ch_identified.filter { sample, alignment -> alignment.name.toLowerCase().endsWith('.bam') }
        CRAM2BAM(ch_crams, INDEX_REFERENCE.out.fasta.first(), INDEX_REFERENCE.out.indexes.first())
        ch_bams = ch_bam_files.mix(
            CRAM2BAM.out.bam.map { sample, bam, bai -> tuple(sample, bam) }
        )
        ANGSD_CONSENSUS(ch_bams)
        ch_consensus = ANGSD_CONSENSUS.out.consensus
    } else {
        ch_vcf = Channel
            .fromPath(params.vcf, checkIfExists: true)
            .map { vcf ->
                def vcf_files = [vcf]
                if (vcf.name.toLowerCase().endsWith('.gz')) {
                    ['.tbi', '.csi'].each { suffix ->
                        def index = file("${vcf}${suffix}")
                        if (index.exists()) vcf_files << index
                    }
                }
                tuple(vcf.name, vcf_files)
            }

        INDEX_VCF(ch_vcf)
        VALIDATE_VCF_SAMPLES(INDEX_VCF.out.vcf, INDEX_VCF.out.index, ch_sample_names)
        CONSENSUS(ch_samples.map { sample, mask -> sample },
            INDEX_REFERENCE.out.fasta.first(), INDEX_REFERENCE.out.indexes.first(),
            VALIDATE_VCF_SAMPLES.out.vcf.first(), VALIDATE_VCF_SAMPLES.out.index.first())
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

    EXTRACT_WINDOW_SEQ(ch_window_consensus,
        params.include_reference ? INDEX_REFERENCE.out.fasta.first() : Channel.value([]),
        params.include_reference ? INDEX_REFERENCE.out.indexes.first() : Channel.value([]))

    // Collect trees (only from alignments that passed filtering)
    ch_trees = EXTRACT_WINDOW_SEQ.out.tree
        .filter { window_id, tree -> tree.size() > 0 }
        .map { window_id, tree -> tree }

    // ---- species tree ----
    COMBINE_TREES(ch_trees.collect())
    ASTRAL(COMBINE_TREES.out.trees)

    // ---- summary files ----
    SUMMARIZE_FILTERING(EXTRACT_WINDOW_SEQ.out.filter_summary.collect())

    ch_collected_alignments = EXTRACT_WINDOW_SEQ.out.alignment
        .filter { window_id, fasta -> fasta.size() > 0 }
        .map { window_id, fasta -> fasta }
        .toList()

    SUMMARIZE_SAMPLE_SURVIVAL(ch_sample_names, ch_collected_alignments)

    publish: 
        species_tree = ASTRAL.out.tree
        locus_trees = COMBINE_TREES.out.trees
        filtering_summary = SUMMARIZE_FILTERING.out.summary
        sample_survival = SUMMARIZE_SAMPLE_SURVIVAL.out.summary

    }

    output {
        species_tree {
            path "${params.name}"
        }
        locus_trees {
            path "${params.name}"
        }
        filtering_summary {
            path "${params.name}"
        }
        sample_survival {
            path "${params.name}"
        }
     }
    
