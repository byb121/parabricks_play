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

nextflow.enable.dsl = 2

include { PARABRICKS_FQ2BAM                           } from './modules/nf-core/parabricks/fq2bam/main'
include { MUTECT2_TUMOR_NORMAL_SOMATIC_GPU           } from './subworkflows/local/mutect2_tumor_normal_somatic_gpu'
include { PARABRICKS_APPLYBQSR                       } from './modules/nf-core/parabricks/applybqsr/main'

workflow {

    // Input channels
    ch_tumor_fastq = Channel.fromFilePairs(params.tumor_fastq_pattern, checkIfExists: true)
        .map { sample_id, fastq_files ->
            def meta = [id: sample_id, sample: sample_id, single_end: false]
            [meta, fastq_files]
        }

    ch_normal_fastq = Channel.fromFilePairs(params.normal_fastq_pattern, checkIfExists: true)
        .map { sample_id, fastq_files ->
            def meta = [id: sample_id, sample: sample_id, single_end: false]
            [meta, fastq_files]
        }

    // Reference genome
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
        .map { fasta ->
            [[id: 'genome'], fasta]
        }

    ch_fasta_index = Channel.fromPath(params.fasta_fai, checkIfExists: true)
        .map { fai ->
            [[id: 'genome'], fai]
        }

    ch_dict = Channel.fromPath(params.dict, checkIfExists: true)
        .map { dict ->
            [[id: 'genome'], dict]
        }

    // Optional reference files for variant calling
    ch_alleles = Channel.fromPath(params.alleles, checkIfExists: true)
        .map { alleles ->
            [[id: 'alleles'], alleles]
        }

    ch_alleles_tbi = Channel.fromPath(params.alleles_tbi, checkIfExists: true)
        .map { tbi ->
            [[id: 'alleles'], tbi]
        }

    ch_germline_resource = Channel.fromPath(params.germline_resource, checkIfExists: true)
        .map { resource ->
            [[id: 'germline'], resource]
        }

    ch_germline_resource_tbi = Channel.fromPath(params.germline_resource_tbi, checkIfExists: true)
        .map { tbi ->
            [[id: 'germline'], tbi]
        }

    ch_panel_of_normals = Channel.fromPath(params.panel_of_normals, checkIfExists: true)
        .map { pon ->
            [[id: 'pon'], pon]
        }

    ch_panel_of_normals_tbi = Channel.fromPath(params.panel_of_normals_tbi, checkIfExists: true)
        .map { tbi ->
            [[id: 'pon'], tbi]
        }

    ch_interval_file = Channel.fromPath(params.interval_file, checkIfExists: true)

    // BWA index
    ch_bwa_index = Channel.fromPath(params.bwa_index, checkIfExists: true)
        .map { index ->
            [[id: 'genome'], index]
        }

    // Known sites for BQSR
    ch_known_sites = Channel.fromPath(params.known_sites, checkIfExists: true)

    // Align tumor sample with fq2bam
    PARABRICKS_FQ2BAM(
        ch_tumor_fastq,
        ch_fasta,
        ch_bwa_index,
        ch_interval_file,
        ch_known_sites,
        'bam'
    )

    // Set output name for tumor BAM
    ch_tumor_bam = PARABRICKS_FQ2BAM.out.bam.map { meta, bam ->
        [[id: "${meta.id}_tumor"], bam]
    }

    ch_tumor_bai = PARABRICKS_FQ2BAM.out.bai.map { meta, bai ->
        [[id: "${meta.id}_tumor"], bai]
    }

    // Align normal sample with fq2bam
    PARABRICKS_FQ2BAM(
        ch_normal_fastq,
        ch_fasta,
        ch_bwa_index,
        ch_interval_file,
        ch_known_sites,
        'bam'
    )

    // Set output name for normal BAM
    ch_normal_bam = PARABRICKS_FQ2BAM.out.bam.map { meta, bam ->
        [[id: "${meta.id}_normal"], bam]
    }

    ch_normal_bai = PARABRICKS_FQ2BAM.out.bai.map { meta, bai ->
        [[id: "${meta.id}_normal"], bai]
    }

    // Combine tumor and normal BAMs for MUTECT2 workflow
    // Input format: [ val(meta), path(input), path(input_index), val(which_norm) ]
    ch_bam_pair = ch_tumor_bam
        .join(ch_tumor_bai)
        .combine(
            ch_normal_bam.join(ch_normal_bai),
            by: 0  // This won't work directly, need a different approach
        )

    // Alternative approach: Create a combined channel with sample linking
    // Cross product of tumor and normal samples
    ch_bam_pairs = ch_tumor_bam
        .join(ch_tumor_bai)
        .map { meta, tumor_bam, tumor_bai ->
            [sample: meta.id.replaceAll('_tumor', ''), meta: meta, tumor_bam: tumor_bam, tumor_bai: tumor_bai]
        }
        .combine(
            ch_normal_bam
                .join(ch_normal_bai)
                .map { meta, normal_bam, normal_bai ->
                    [sample: meta.id.replaceAll('_normal', ''), meta: meta, normal_bam: normal_bam, normal_bai: normal_bai]
                },
            by: 0  // Match by sample ID
        )
        .map { sample, tumor_data, normal_data ->
            def combined_meta = [id: sample, single_end: false]
            [
                combined_meta,
                [tumor_data.tumor_bam, normal_data.normal_bam],
                [tumor_data.tumor_bai, normal_data.normal_bai],
                0  // 0 indicates tumor is reference (first input)
            ]
        }

    // Call somatic variants using MUTECT2 workflow
    MUTECT2_TUMOR_NORMAL_SOMATIC_GPU(
        ch_bam_pairs,
        ch_fasta,
        ch_fasta_index,
        ch_dict,
        ch_alleles,
        ch_alleles_tbi,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_panel_of_normals,
        ch_panel_of_normals_tbi,
        ch_interval_file
    )

    // Output results
    MUTECT2_TUMOR_NORMAL_SOMATIC_GPU.out.filtered_vcf.view { meta, vcf ->
        "Filtered VCF: ${vcf}"
    }
}

workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
}
