version 1.0

import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/assembly/task_ivar_consensus.wdl" as consensus_task
import "../../tasks/assembly/task_ivar_primer_trim.wdl" as primer_trim_task
import "../../tasks/assembly/task_ivar_variant_call.wdl" as variant_call_task
import "../../tasks/quality_control/task_assembly_metrics.wdl" as assembly_metrics

workflow ivar_consensus {
  meta {
    description: "Reference-based consensus calling using ivar"
  }
  input {
    String samplename
    File read1
    File? read2
    File? reference_genome
    File? reference_gff
    Boolean trim_primers
    Int min_depth
    Float variant_min_freq 
    Float consensus_min_freq 
    File? primer_bed
  }
  call bwa_task.bwa {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      reference_genome = reference_genome
  }
  if (trim_primers){
    call primer_trim_task.primer_trim {
      input:
        samplename = samplename,
        primer_bed = select_first([primer_bed]),
        bamfile = bwa.sorted_bam
    }
    call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
    input:
      samplename = samplename,
      bamfile = primer_trim.trim_sorted_bam,
    }
  }
  call variant_call_task.variant_call {
    input:
      samplename = samplename,
      bamfile = select_first([primer_trim.trim_sorted_bam, bwa.sorted_bam]),
      reference_gff = reference_gff,
      reference_genome = reference_genome,
      variant_min_depth = min_depth,
      variant_min_freq = variant_min_freq

  }
  call consensus_task.consensus {
    input:
      samplename = samplename,
      bamfile = select_first([primer_trim.trim_sorted_bam, bwa.sorted_bam]),
      reference_genome = reference_genome,
      consensus_min_depth = min_depth,
      consensus_min_freq = consensus_min_freq
  }
  call assembly_metrics.stats_n_coverage {
    input:
      samplename = samplename,
      bamfile = bwa.sorted_bam
  }
  output {
    # bwa outputs
    String bwa_version = bwa.bwa_version
    String samtools_version = bwa.sam_version
    File read1_aligned = bwa.read1_aligned
    File? read2_aligned = bwa.read2_aligned
    String aligned_bam =  select_first([primer_trim.trim_sorted_bam, bwa.sorted_bam, ""]) # set default values for select_first() to avoid workflow failures
    String aligned_bai = select_first([primer_trim.trim_sorted_bai, bwa.sorted_bai, ""])
    
    # primer trimming outputs
    Float? primer_trimmed_read_percent = primer_trim.primer_trimmed_read_percent
    String? ivar_version_primtrim = primer_trim.ivar_version
    String? samtools_version_primtrim = primer_trim.samtools_version
    String? primer_bed_name = primer_trim.primer_bed_name
    
    # variant call outputs
    File ivar_tsv = variant_call.sample_variants_tsv
    File ivar_vcf = variant_call.sample_variants_vcf
    String ivar_variant_version = variant_call.ivar_version

    # assembly outputs
    String assembly_method_nonflu = "~{bwa.bwa_version}; ~{primer_trim.ivar_version}"
    File assembly_fasta = consensus.consensus_seq
    String ivar_version_consensus = consensus.ivar_version
    String samtools_version_consensus = consensus.samtools_version

    # consensus qc outputs
    Int consensus_n_variant_min_depth = min_depth
    File consensus_stats = stats_n_coverage.stats
    File consensus_flagstat = stats_n_coverage.flagstat
    String meanbaseq_trim = select_first([stats_n_coverage_primtrim.meanbaseq, stats_n_coverage.meanbaseq,""])
    String meanmapq_trim = select_first([stats_n_coverage_primtrim.meanmapq, stats_n_coverage.meanmapq,""])
    String assembly_mean_coverage = select_first([stats_n_coverage_primtrim.depth, stats_n_coverage.depth,""])
    String samtools_version_stats = stats_n_coverage.samtools_version
  
  }
}