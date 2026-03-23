// Run PARABRICKS_MUTECTCALLER in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls

include { GATK4_CALCULATECONTAMINATION                                } from '../../../modules/nf-core/gatk4/calculatecontamination'
include { GATK4_FILTERMUTECTCALLS                                     } from '../../../modules/nf-core/gatk4/filtermutectcalls'
include { GATK4_GETPILEUPSUMMARIES as GATK4_GETPILEUPSUMMARIES_NORMAL } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as GATK4_GETPILEUPSUMMARIES_TUMOR  } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_LEARNREADORIENTATIONMODEL                             } from '../../../modules/nf-core/gatk4/learnreadorientationmodel'
include { PARABRICKS_MUTECTCALLER                                     } from '../../../modules/nf-core/parabricks/mutectcaller'

workflow MUTECT2_TUMOR_NORMAL_SOMATIC_GPU {

    take:
    ch_input // channel: [ val(meta), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]
    ch_dict // channel: [ val(meta), path(dict) ]
    ch_alleles // channel: /path/to/alleles
    ch_alleles_tbi // channel: /path/to/alleles/index
    ch_germline_resource // channel: /path/to/germline/resource
    ch_germline_resource_tbi // channel: /path/to/germline/index
    ch_panel_of_normals // channel: /path/to/panel/of/normals
    ch_panel_of_normals_tbi // channel: /path/to/panel/of/normals/index
    ch_interval_file // channel: /path/to/interval/file

    main:
    // Perform variant calling using PARABRICKS_MUTECTCALLER module in tumor single mode.
    PARABRICKS_MUTECTCALLER(
        ch_input.combine(ch_interval_file.first()), // Combine interval file with input channel to pass intervals to the module
        ch_fasta,
        ch_alleles,
        ch_alleles_tbi,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_panel_of_normals,
        ch_panel_of_normals_tbi,
    )

    // Generate artifactpriors using learnreadorientationmodel on the f1r2 output of PARABRICKS_MUTECTCALLER.
    GATK4_LEARNREADORIENTATIONMODEL(PARABRICKS_MUTECTCALLER.out.f1r2)

    // Generate pileup summary tables using getepileupsummaries
    // Tumor sample should always be passed in as the first input and input list entries of ch_input,
    // to ensure correct file order for calculatecontamination.
    ch_pileup_tumor_input = ch_input
        .combine(ch_interval_file)
        .map { meta, tumor_bam, tumor_bam_index, normal_bam, normal_bam_index, intervals ->
            [meta, tumor_bam, tumor_bam_index, intervals]
        }

    ch_pileup_normal_input = ch_input
        .combine(ch_interval_file)
        .map { meta, tumor_bam, tumor_bam_index, normal_bam, normal_bam_index, intervals ->
            [meta, normal_bam, normal_bam_index, intervals]
        }

    GATK4_GETPILEUPSUMMARIES_TUMOR(
        ch_pileup_tumor_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi,
    )

    GATK4_GETPILEUPSUMMARIES_NORMAL(
        ch_pileup_normal_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi,
    )

    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    ch_pileup_tumor = GATK4_GETPILEUPSUMMARIES_TUMOR.out.table.collect()
    ch_pileup_normal = GATK4_GETPILEUPSUMMARIES_NORMAL.out.table.collect()
    ch_calccon_in = ch_pileup_tumor.join(ch_pileup_normal, failOnDuplicate: true, failOnMismatch: true)

    GATK4_CALCULATECONTAMINATION(ch_calccon_in)

    // PARABRICKS_MUTECTCALLER calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables.
    ch_vcf = PARABRICKS_MUTECTCALLER.out.vcf.collect()
    ch_tbi = PARABRICKS_MUTECTCALLER.out.tbi.collect()
    ch_stats = PARABRICKS_MUTECTCALLER.out.stats.collect()
    ch_orientation = GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior.collect()
    ch_segment = GATK4_CALCULATECONTAMINATION.out.segmentation.collect()

    // [] is used as a placeholder for optional input to specify the contamination estimate as a value, since the contamination table is used, this is not needed.
    ch_contamination = GATK4_CALCULATECONTAMINATION.out.contamination.map { meta, table -> [meta, table, []] }.collect()

    ch_filtermutect_in = ch_vcf
        .join(ch_tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_stats, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_orientation, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_segment, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_contamination, failOnDuplicate: true, failOnMismatch: true)

    GATK4_FILTERMUTECTCALLS(
        ch_filtermutect_in,
        ch_fasta,
        ch_fai,
        ch_dict,
    )

    emit:
    artifact_priors     = GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior // channel: [ val(meta), path(artifactprior) ]
    contamination_table = GATK4_CALCULATECONTAMINATION.out.contamination // channel: [ val(meta), path(table) ]
    filtered_stats      = GATK4_FILTERMUTECTCALLS.out.stats // channel: [ val(meta), path(stats) ]
    filtered_tbi        = GATK4_FILTERMUTECTCALLS.out.tbi // channel: [ val(meta), path(tbi) ]
    filtered_vcf        = GATK4_FILTERMUTECTCALLS.out.vcf // channel: [ val(meta), path(vcf) ]
    mutect2_f1r2        = PARABRICKS_MUTECTCALLER.out.f1r2 // channel: [ val(meta), path(f1r2) ]
    mutect2_stats       = PARABRICKS_MUTECTCALLER.out.stats // channel: [ val(meta), path(stats) ]
    mutect2_tbi         = PARABRICKS_MUTECTCALLER.out.tbi // channel: [ val(meta), path(tbi) ]
    mutect2_vcf         = PARABRICKS_MUTECTCALLER.out.vcf // channel: [ val(meta), path(vcf) ]
    pileup_table_normal = GATK4_GETPILEUPSUMMARIES_NORMAL.out.table // channel: [ val(meta), path(table) ]
    pileup_table_tumor  = GATK4_GETPILEUPSUMMARIES_TUMOR.out.table // channel: [ val(meta), path(table) ]
    segmentation_table  = GATK4_CALCULATECONTAMINATION.out.segmentation // channel: [ val(meta), path(table) ]
}