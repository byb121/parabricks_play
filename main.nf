#!/usr/bin/env nextflow

/*
================================================================================
    Parabricks Somatic Variant Calling Pipeline
================================================================================
    A Nextflow pipeline that:
    1. Takes paired-end FASTQ files for tumor and normal samples
    2. Aligns them to reference genome using PARABRICKS fq2bam
    3. Calls somatic variants using MUTECT2_TUMOR_NORMAL_SOMATIC_GPU workflow
================================================================================
*/

include { PARABRICKS_FQ2BAM                          } from './modules/nf-core/parabricks/fq2bam/main'
include { MUTECT2_TUMOR_NORMAL_SOMATIC_GPU           } from './subworkflows/local/mutect2_tumor_normal_somatic_gpu'
include { PARABRICKS_APPLYBQSR                       } from './modules/nf-core/parabricks/applybqsr/main'

workflow {

    // read in sample sheet and create channels for tumor and normal FASTQ files
    ch_samples = Channel.fromPath(params.input_samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .flatMap { row ->
            def tumor_meta = [id: row.tumor_id, sample_type: 'tumor', sample_name: row.sample_name, single_end: false]
            def normal_meta = [id: row.normal_id, sample_type: 'normal', sample_name: row.sample_name, single_end: false]
            def tumor_fastq_files = [tumor_meta, [file(row.tumor_fastq_1), file(row.tumor_fastq_2)]]
            def normal_fastq_files = [normal_meta, [file(row.normal_fastq_1), file(row.normal_fastq_2)]]
            [tumor_fastq_files, normal_fastq_files]
        }

    // Reference genome
    ch_fasta = Channel.value([
        [id: 'genome'],
        file(params.fasta, checkIfExists: true, glob: false)
    ])
    ch_interval_file = Channel.value([[id: 'intervals'], []]) // No intervals for now, can be set to a channel of interval files if needed
    ch_bwa_index = Channel.value([
        [id: 'genome'],
        file(params.bwa_index, checkIfExists: true, glob: false)
    ])
    ch_mutect_interval = Channel.value([
        [id: 'mutect_intervals'],
        file(params.mutect_interval, checkIfExists: true, glob: false)
    ])
    ch_dbsnp = Channel.value([
        [id: 'dbsnp'],
        file(params.dbsnp, checkIfExists: true, glob: false)
    ])
    ch_dnsnp_tbi = Channel.value([
        [id: 'dbsnp_index'],
        file(params.dbsnp_tbi, checkIfExists: true, glob: false)
    ])
    ch_known_indels = Channel.value([
        [id: 'known_indels'],
        file(params.known_indels, checkIfExists: true, glob: false)
    ])
    ch_known_indels_tbi = Channel.value([
        [id: 'known_indels_index'],
        file(params.known_indels_tbi, checkIfExists: true, glob: false)
    ])
    ch_mills_and_1000G_gold_standard = Channel.value([
        [id: 'mills_and_1000G_gold_standard'],
        file(params.mills_and_1000G_gold_standard, checkIfExists: true, glob: false)
    ])
    ch_mills_and_1000G_gold_standard_tbi = Channel.value([
        [id: 'mills_and_1000G_gold_standard_index'],
        file(params.mills_and_1000G_gold_standard_tbi, checkIfExists: true, glob: false)
    ])

    ch_known_sites = ch_dbsnp
        .mix(ch_dnsnp_tbi)
        .mix(ch_known_indels)
        .mix(ch_known_indels_tbi)
        .mix(ch_mills_and_1000G_gold_standard)
        .mix(ch_mills_and_1000G_gold_standard_tbi)
        .collect(flat: false)
        .map { dbsnp, dnsnp_tbi, indels, indels_tbi, gold_standard, gold_standard_tbi ->
            [[id: 'known_sites'], [dbsnp[1], dnsnp_tbi[1], indels[1], indels_tbi[1], gold_standard[1], gold_standard_tbi[1]]]
        }

    // Align sample fastqs with fq2bam
    PARABRICKS_FQ2BAM(
        ch_samples,
        ch_fasta,
        ch_bwa_index,
        ch_interval_file,
        ch_known_sites,
        'bam'
    )

    ch_apply_bqsr_input = PARABRICKS_FQ2BAM.out.bam
        .join(PARABRICKS_FQ2BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)
        .join(PARABRICKS_FQ2BAM.out.bqsr_table, failOnDuplicate: true, failOnMismatch: true)

    // apply BQSR with PARABRICKS_APPLYBQSR
    PARABRICKS_APPLYBQSR(   
        ch_apply_bqsr_input,
        ch_interval_file,
        ch_fasta
    )

    // reference files for variant calling
    ch_fasta_fai = Channel.value([
        [id: 'genome'],
        file(params.fasta_fai, checkIfExists: true, glob: false)
    ])
    ch_fasta_dict = Channel.value([
        [id: 'genome'],
        file(params.dict, checkIfExists: true, glob: false)
    ])
    ch_alleles = Channel.value([])
    ch_alleles_tbi = Channel.value([])
    ch_germline_resource = Channel.value(file(params.germline_resource, checkIfExists: true, glob: false))
    ch_germline_resource_tbi = Channel.value(file(params.germline_resource_tbi, checkIfExists: true, glob: false))
    ch_panel_of_normals = Channel.value(file(params.panel_of_normals, checkIfExists: true, glob: false))
    ch_panel_of_normals_tbi = Channel.value(file(params.panel_of_normals_tbi, checkIfExists: true, glob: false))

    // Set output name for tumor BAM
    ch_recalibrated_bam = PARABRICKS_APPLYBQSR.out.bam
        .join(PARABRICKS_APPLYBQSR.out.bai, failOnDuplicate: true, failOnMismatch: true)
        .branch { meta, bam, bai ->
            def new_meta = [id: meta.sample_name]
            tumor: meta.sample_type == 'tumor'
                return [new_meta, meta, bam, bai]
            normal: true
                return [new_meta, meta, bam, bai]
        }

    ch_mutect2_input = ch_recalibrated_bam.tumor
        .join(ch_recalibrated_bam.normal, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
            def new_meta = [id: meta.id, normal_id: normal_meta.id, tumor_id: tumor_meta.id]
            [new_meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }
    
    MUTECT2_TUMOR_NORMAL_SOMATIC_GPU(
        ch_mutect2_input,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_alleles,
        ch_alleles_tbi,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_panel_of_normals,
        ch_panel_of_normals_tbi,
        ch_mutect_interval
    )

    // // Output results
    // MUTECT2_TUMOR_NORMAL_SOMATIC_GPU.out.filtered_vcf.view { meta, vcf ->
    //     "Filtered VCF: ${vcf}"
    // }
}

workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
}
