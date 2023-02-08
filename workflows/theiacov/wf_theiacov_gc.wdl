version 1.0

# Workflows from Theiagen's public_health_viral_genomics
# Source: https://github.com/theiagen/public_health_viral_genomics
import "wf_theiacov_clearlabs.wdl" as clearlabs
import "wf_theiacov_illumina_pe.wdl" as illumina_pe
import "wf_theiacov_illumina_se.wdl" as illumina_se
import "wf_theiacov_ont.wdl" as ont
import "../tasks/task_theiacov_summary.wdl" as summary

struct parseJSON {
  String sample
  String theiacov_wf
  File   r1
  File   r2
  File   primers
}
workflow theiacov_gc {
  meta {
    description: "Incorporates each of the TheiaCoV workflows (clearlabs, illumina_pe, illumina_se, ont) into a single run."
    author: "Robert A. Petit III"
    email:  "robert.petit@theiagen.com"
  }
  input {
    Array[parseJSON] samples
  }

  scatter (sample in samples) {
    if (sample.theiacov_wf == "clearlabs") {
      call clearlabs.theiacov_clearlabs as theiacov_clearlabs {
        input:
          samplename = sample.sample,
          clear_lab_fastq = sample.r1,
          primer_bed = sample.primers
      }
      call summary.theiacov_summary as clearlabs_summary {
        input:
          samplename = sample.sample,
          theiacov_workflow = 'theiacov_clearlabs',
          theiacov_version = theiacov_clearlabs.theiacov_clearlabs_version,
          theiacov_analysis_date = theiacov_clearlabs.theiacov_clearlabs_analysis_date,
          seq_platform = theiacov_clearlabs.seq_platform,
          primer_bed_name = theiacov_clearlabs.primer_bed_name,
          percent_reference_coverage = theiacov_clearlabs.percent_reference_coverage,
          number_N = theiacov_clearlabs.number_N,
          pango_lineage = theiacov_clearlabs.pango_lineage,
          pangolin_conflicts = theiacov_clearlabs.pangolin_conflicts,
          pangolin_notes = theiacov_clearlabs.pangolin_notes,
          pangolin_assignment_version = theiacov_clearlabs.pangolin_assignment_version,
          pangolin_docker = theiacov_clearlabs.pangolin_docker,
          pango_lineage_expanded = theiacov_clearlabs.pango_lineage_expanded,
          pangolin_versions = theiacov_clearlabs.pangolin_versions,
          nextclade_clade = theiacov_clearlabs.nextclade_clade,
          nextclade_aa_subs = theiacov_clearlabs.nextclade_aa_subs,
          nextclade_aa_dels = theiacov_clearlabs.nextclade_aa_dels,
          vadr_num_alerts = theiacov_clearlabs.vadr_num_alerts,
          assembly_length_unambiguous = theiacov_clearlabs.assembly_length_unambiguous,
          assembly_mean_coverage = theiacov_clearlabs.assembly_mean_coverage,
          s_gene_mean_coverage = theiacov_clearlabs.sc2_s_gene_mean_coverage,
          assembly_method = theiacov_clearlabs.assembly_method,
          number_Degenerate = theiacov_clearlabs.number_Degenerate,
          number_Total = theiacov_clearlabs.number_Total,
          meanbaseq_trim = theiacov_clearlabs.meanbaseq_trim,
          meanmapq_trim = theiacov_clearlabs.meanmapq_trim,
          num_reads_clean1 = theiacov_clearlabs.num_reads_clean,
          num_reads_raw1 = theiacov_clearlabs.num_reads_raw,
          fastq_scan_version = theiacov_clearlabs.fastq_scan_version,
          kraken_human = theiacov_clearlabs.kraken_human,
          kraken_human_dehosted = theiacov_clearlabs.kraken_human_dehosted,
          kraken_sc2 = theiacov_clearlabs.kraken_sc2,
          kraken_sc2_dehosted = theiacov_clearlabs.kraken_sc2_dehosted,
          artic_version = theiacov_clearlabs.artic_version,
          artic_docker = theiacov_clearlabs.artic_docker,
          medaka_reference = theiacov_clearlabs.medaka_reference,
          kraken_version = theiacov_clearlabs.kraken_version,
          nextclade_version = theiacov_clearlabs.nextclade_version,
          nextclade_docker = theiacov_clearlabs.nextclade_docker,
          nextclade_ds_tag = theiacov_clearlabs.nextclade_ds_tag,
          samtools_version = theiacov_clearlabs.samtools_version,
          vadr_docker = theiacov_clearlabs.vadr_docker
      }
    }
    if (sample.theiacov_wf == "illumina_pe") {
      call illumina_pe.theiacov_illumina_pe as theiacov_illumina_pe {
        input:
          samplename = sample.sample,
          read1_raw = sample.r1,
          read2_raw = sample.r2,
          primer_bed = sample.primers
      }
      call summary.theiacov_summary as illumina_pe_summary {
        input:
          samplename = sample.sample,
          theiacov_workflow = 'theiacov_illumina_pe',
          theiacov_version = theiacov_illumina_pe.theiacov_illumina_pe_version,
          theiacov_analysis_date = theiacov_illumina_pe.theiacov_illumina_pe_analysis_date,
          seq_platform = theiacov_illumina_pe.seq_platform,
          primer_bed_name = theiacov_illumina_pe.primer_bed_name,
          percent_reference_coverage = theiacov_illumina_pe.percent_reference_coverage,
          consensus_n_variant_min_depth = theiacov_illumina_pe.consensus_n_variant_min_depth,
          number_N = theiacov_illumina_pe.number_N,
          pango_lineage = theiacov_illumina_pe.pango_lineage,
          pangolin_conflicts = theiacov_illumina_pe.pangolin_conflicts,
          pangolin_notes = theiacov_illumina_pe.pangolin_notes,
          pangolin_assignment_version = theiacov_illumina_pe.pangolin_assignment_version,
          pango_lineage_expanded = theiacov_illumina_pe.pango_lineage_expanded,
          pangolin_versions = theiacov_illumina_pe.pangolin_versions,
          pangolin_docker = theiacov_illumina_pe.pangolin_docker,
          nextclade_clade = theiacov_illumina_pe.nextclade_clade,
          nextclade_aa_subs = theiacov_illumina_pe.nextclade_aa_subs,
          nextclade_aa_dels = theiacov_illumina_pe.nextclade_aa_dels,
          vadr_num_alerts = theiacov_illumina_pe.vadr_num_alerts,
          assembly_length_unambiguous = theiacov_illumina_pe.assembly_length_unambiguous,
          assembly_mean_coverage = theiacov_illumina_pe.assembly_mean_coverage,
          s_gene_mean_coverage = theiacov_illumina_pe.sc2_s_gene_mean_coverage,
          s_gene_percent_coverage = theiacov_illumina_pe.sc2_s_gene_percent_coverage,
          assembly_method = theiacov_illumina_pe.assembly_method,
          number_Degenerate = theiacov_illumina_pe.number_Degenerate,
          number_Total = theiacov_illumina_pe.number_Total,
          meanbaseq_trim = theiacov_illumina_pe.meanbaseq_trim,
          meanmapq_trim = theiacov_illumina_pe.meanmapq_trim,
          num_reads_clean1 = theiacov_illumina_pe.num_reads_clean1,
          num_reads_clean2 = theiacov_illumina_pe.num_reads_clean2,
          num_reads_clean_pairs = theiacov_illumina_pe.num_reads_clean_pairs,
          num_reads_raw1 = theiacov_illumina_pe.num_reads_raw1,
          num_reads_raw2 = theiacov_illumina_pe.num_reads_raw2,
          num_reads_raw_pairs = theiacov_illumina_pe.num_reads_raw_pairs,
          kraken_human = theiacov_illumina_pe.kraken_human,
          kraken_human_dehosted = theiacov_illumina_pe.kraken_human_dehosted,
          kraken_sc2 = theiacov_illumina_pe.kraken_sc2,
          kraken_sc2_dehosted = theiacov_illumina_pe.kraken_sc2_dehosted,
          primer_trimmed_read_percent = theiacov_illumina_pe.primer_trimmed_read_percent,
          bbduk_docker = theiacov_illumina_pe.bbduk_docker,
          bwa_version = theiacov_illumina_pe.bwa_version,
          fastq_scan_version = theiacov_illumina_pe.fastq_scan_version,
          ivar_variant_version = theiacov_illumina_pe.ivar_variant_version,
          ivar_version_consensus = theiacov_illumina_pe.ivar_version_consensus,
          ivar_version_primtrim = theiacov_illumina_pe.ivar_version_primtrim,
          kraken_version = theiacov_illumina_pe.kraken_version,
          nextclade_version = theiacov_illumina_pe.nextclade_version,
          nextclade_docker = theiacov_illumina_pe.nextclade_docker,
          nextclade_ds_tag = theiacov_illumina_pe.nextclade_ds_tag,
          samtools_version = theiacov_illumina_pe.samtools_version,
          samtools_version_consensus = theiacov_illumina_pe.samtools_version_consensus,
          samtools_version_primtrim = theiacov_illumina_pe.samtools_version_primtrim,
          samtools_version_stats = theiacov_illumina_pe.samtools_version_stats,
          trimmomatic_version = theiacov_illumina_pe.trimmomatic_version,
          vadr_docker = theiacov_illumina_pe.vadr_docker
      }
    }
    if (sample.theiacov_wf == "illumina_se") {
      call illumina_se.theiacov_illumina_se as theiacov_illumina_se {
        input:
          samplename = sample.sample,
          read1_raw  = sample.r1,
          primer_bed = sample.primers
      }
      call summary.theiacov_summary as illumina_se_summary {
        input:
          samplename = sample.sample,
          theiacov_workflow = 'theiacov_illumina_se',
          theiacov_version = theiacov_illumina_se.theiacov_illumina_se_version,
          theiacov_analysis_date = theiacov_illumina_se.theiacov_illumina_se_analysis_date,
          seq_platform = theiacov_illumina_se.seq_platform,
          primer_bed_name = theiacov_illumina_se.primer_bed_name,
          percent_reference_coverage = theiacov_illumina_se.percent_reference_coverage,
          consensus_n_variant_min_depth = theiacov_illumina_se.consensus_n_variant_min_depth,
          number_N = theiacov_illumina_se.number_N,
          pango_lineage = theiacov_illumina_se.pango_lineage,
          pangolin_conflicts = theiacov_illumina_se.pangolin_conflicts,
          pangolin_notes = theiacov_illumina_se.pangolin_notes,
          pango_lineage_expanded = theiacov_illumina_se.pango_lineage_expanded,
          pangolin_versions = theiacov_illumina_se.pangolin_versions,
          pangolin_assignment_version = theiacov_illumina_se.pangolin_assignment_version,
          pangolin_docker = theiacov_illumina_se.pangolin_docker,
          nextclade_clade = theiacov_illumina_se.nextclade_clade,
          nextclade_aa_subs = theiacov_illumina_se.nextclade_aa_subs,
          nextclade_aa_dels = theiacov_illumina_se.nextclade_aa_dels,
          vadr_num_alerts = theiacov_illumina_se.vadr_num_alerts,
          assembly_length_unambiguous = theiacov_illumina_se.assembly_length_unambiguous,
          assembly_mean_coverage = theiacov_illumina_se.assembly_mean_coverage,
          s_gene_mean_coverage = theiacov_illumina_se.sc2_s_gene_mean_coverage,
          s_gene_percent_coverage = theiacov_illumina_se.sc2_s_gene_percent_coverage,
          assembly_method = theiacov_illumina_se.assembly_method,
          number_Degenerate = theiacov_illumina_se.number_Degenerate,
          number_Total = theiacov_illumina_se.number_Total,
          meanbaseq_trim = theiacov_illumina_se.meanbaseq_trim,
          meanmapq_trim = theiacov_illumina_se.meanmapq_trim,
          num_reads_clean1 = theiacov_illumina_se.num_reads_clean,
          num_reads_raw1 = theiacov_illumina_se.num_reads_raw,
          fastq_scan_version = theiacov_illumina_se.fastq_scan_version,
          kraken_human = theiacov_illumina_se.kraken_human,
          kraken_sc2 = theiacov_illumina_se.kraken_sc2,
          primer_trimmed_read_percent = theiacov_illumina_se.primer_trimmed_read_percent,
          bbduk_docker = theiacov_illumina_se.bbduk_docker,
          bwa_version = theiacov_illumina_se.bwa_version,
          fastq_scan_version = theiacov_illumina_se.fastq_scan_version,
          ivar_variant_version = theiacov_illumina_se.ivar_variant_version,
          ivar_version_consensus = theiacov_illumina_se.ivar_version_consensus,
          ivar_version_primtrim = theiacov_illumina_se.ivar_version_primtrim,
          kraken_version = theiacov_illumina_se.kraken_version,
          nextclade_version = theiacov_illumina_se.nextclade_version,
          nextclade_docker = theiacov_illumina_se.nextclade_docker,
          nextclade_ds_tag = theiacov_illumina_se.nextclade_ds_tag,
          samtools_version = theiacov_illumina_se.samtools_version,
          samtools_version_consensus = theiacov_illumina_se.samtools_version_consensus,
          samtools_version_primtrim = theiacov_illumina_se.samtools_version_primtrim,
          samtools_version_stats = theiacov_illumina_se.samtools_version_stats,
          trimmomatic_version = theiacov_illumina_se.trimmomatic_version,
          vadr_docker = theiacov_illumina_se.vadr_docker
      }
    }
    if (sample.theiacov_wf == "ont") {
      call ont.theiacov_ont as theiacov_ont {
          input:
              samplename = sample.sample,
              demultiplexed_reads = sample.r1,
              primer_bed = sample.primers
      }
      call summary.theiacov_summary as ont_summary {
        input:
          samplename = sample.sample,
          theiacov_workflow = 'theiacov_ont',
          theiacov_version = theiacov_ont.theiacov_ont_version,
          theiacov_analysis_date = theiacov_ont.theiacov_ont_analysis_date,
          seq_platform = theiacov_ont.seq_platform,
          primer_bed_name = theiacov_ont.primer_bed_name,
          percent_reference_coverage = theiacov_ont.percent_reference_coverage,
          number_N = theiacov_ont.number_N,
          pango_lineage = theiacov_ont.pango_lineage,
          pangolin_conflicts = theiacov_ont.pangolin_conflicts,
          pangolin_notes = theiacov_ont.pangolin_notes,
          pango_lineage_expanded = theiacov_ont.pango_lineage_expanded,
          pangolin_versions = theiacov_ont.pangolin_versions,
          pangolin_assignment_version = theiacov_ont.pangolin_assignment_version,
          pangolin_docker = theiacov_ont.pangolin_docker,
          nextclade_clade = theiacov_ont.nextclade_clade,
          nextclade_aa_subs = theiacov_ont.nextclade_aa_subs,
          nextclade_aa_dels = theiacov_ont.nextclade_aa_dels,
          vadr_num_alerts = theiacov_ont.vadr_num_alerts,
          assembly_length_unambiguous = theiacov_ont.assembly_length_unambiguous,
          assembly_mean_coverage = theiacov_ont.assembly_mean_coverage,
          s_gene_mean_coverage = theiacov_ont.sc2_s_gene_mean_coverage,
          assembly_method = theiacov_ont.assembly_method,
          number_Degenerate = theiacov_ont.number_Degenerate,
          number_Total = theiacov_ont.number_Total,
          meanbaseq_trim = theiacov_ont.meanbaseq_trim,
          meanmapq_trim = theiacov_ont.meanmapq_trim,
          num_reads_clean1 = theiacov_ont.num_reads_clean,
          num_reads_raw1 = theiacov_ont.num_reads_raw,
          fastq_scan_version = theiacov_ont.fastq_scan_version,
          kraken_human = theiacov_ont.kraken_human,
          kraken_human_dehosted = theiacov_ont.kraken_human_dehosted,
          kraken_sc2 = theiacov_ont.kraken_sc2,
          kraken_sc2_dehosted = theiacov_ont.kraken_sc2_dehosted,
          artic_version = theiacov_ont.artic_version,
          artic_docker = theiacov_ont.artic_docker,
          medaka_reference = theiacov_ont.medaka_reference,
          kraken_version = theiacov_ont.kraken_version,
          nextclade_version = theiacov_ont.nextclade_version,
          nextclade_docker = theiacov_ont.nextclade_docker,
          nextclade_ds_tag = theiacov_ont.nextclade_ds_tag,
          samtools_version = theiacov_ont.samtools_version,
          vadr_docker = theiacov_ont.vadr_docker
      }
    }
  }
  call summary.merge_theiacov_summary {
    input:
      clearlabs_summaries = clearlabs_summary.summary,
      illumina_pe_summaries = illumina_pe_summary.summary,
      illumina_se_summaries = illumina_se_summary.summary,
      ont_summaries = ont_summary.summary
  }
  output {
    # TheiaCoV outputs
    File summaries_tsv = merge_theiacov_summary.summaries_tsv
    File summaries_json = merge_theiacov_summary.summaries_json
    Array[File] reads_dehosted = flatten([
      select_all(theiacov_clearlabs.reads_dehosted), select_all(theiacov_illumina_pe.read1_dehosted),
      select_all(theiacov_illumina_pe.read2_dehosted), select_all(theiacov_ont.reads_dehosted)
    ])
    Array[File] aligned_bam = flatten([
      select_all(theiacov_clearlabs.aligned_bam), select_all(theiacov_illumina_pe.aligned_bam),
      select_all(theiacov_illumina_se.aligned_bam), select_all(theiacov_ont.aligned_bam)
    ])
    Array[File] aligned_bai = flatten([
      select_all(theiacov_clearlabs.aligned_bai), select_all(theiacov_illumina_pe.aligned_bai),
      select_all(theiacov_illumina_se.aligned_bai), select_all(theiacov_ont.aligned_bai)
    ])
    Array[File] assembly_fasta = flatten([
      select_all(theiacov_clearlabs.assembly_fasta), select_all(theiacov_illumina_pe.assembly_fasta),
      select_all(theiacov_illumina_se.assembly_fasta), select_all(theiacov_ont.assembly_fasta)
    ])
    Array[File] consensus_stats = flatten([
      select_all(theiacov_clearlabs.consensus_stats), select_all(theiacov_illumina_pe.consensus_stats),
      select_all(theiacov_illumina_se.consensus_stats), select_all(theiacov_ont.consensus_stats)
    ])
    Array[File] consensus_flagstat = flatten([
      select_all(theiacov_clearlabs.consensus_flagstat), select_all(theiacov_illumina_pe.consensus_flagstat),
      select_all(theiacov_illumina_se.consensus_flagstat), select_all(theiacov_ont.consensus_flagstat)
    ])
    Array[File] consensus_variants = flatten([
      select_all(theiacov_clearlabs.variants_from_ref_vcf), select_all(theiacov_illumina_pe.ivar_vcf),
      select_all(theiacov_illumina_se.ivar_vcf), select_all(theiacov_ont.variants_from_ref_vcf)
    ])
    Array[File] pango_lineage_report = flatten([
      select_all(theiacov_clearlabs.pango_lineage_report), select_all(theiacov_illumina_pe.pango_lineage_report),
      select_all(theiacov_illumina_se.pango_lineage_report), select_all(theiacov_ont.pango_lineage_report)
    ])
    Array[File] nextclade_json = flatten([
      select_all(theiacov_clearlabs.nextclade_json), select_all(theiacov_illumina_pe.nextclade_json),
      select_all(theiacov_illumina_se.nextclade_json), select_all(theiacov_ont.nextclade_json)
    ])
    Array[File] auspice_json = flatten([
      select_all(theiacov_clearlabs.auspice_json), select_all(theiacov_illumina_pe.auspice_json),
      select_all(theiacov_illumina_se.auspice_json), select_all(theiacov_ont.auspice_json)
    ])
    Array[File] nextclade_tsv = flatten([
      select_all(theiacov_clearlabs.nextclade_tsv), select_all(theiacov_illumina_pe.nextclade_tsv),
      select_all(theiacov_illumina_se.nextclade_tsv), select_all(theiacov_ont.nextclade_tsv)
    ])
    Array[File] vadr_alerts_list = flatten([
      select_all(theiacov_clearlabs.vadr_alerts_list), select_all(theiacov_illumina_pe.vadr_alerts_list),
      select_all(theiacov_illumina_se.vadr_alerts_list), select_all(theiacov_ont.vadr_alerts_list)
    ])
    Array[File] kraken_report = flatten([
      select_all(theiacov_clearlabs.kraken_report), select_all(theiacov_illumina_pe.kraken_report),
      select_all(theiacov_illumina_se.kraken_report), select_all(theiacov_ont.kraken_report)
    ])
    Array[File] kraken_report_dehosted = flatten([
      select_all(theiacov_clearlabs.kraken_report_dehosted), select_all(theiacov_illumina_pe.kraken_report_dehosted),
      select_all(theiacov_ont.kraken_report_dehosted)
    ])
    Array[File] json_summary = flatten([
      select_all(clearlabs_summary.summary), select_all(illumina_pe_summary.summary),
      select_all(illumina_se_summary.summary), select_all(ont_summary.summary)
    ])
  }
}