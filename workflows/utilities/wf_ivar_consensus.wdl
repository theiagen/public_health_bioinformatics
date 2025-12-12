version 1.0

import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/assembly/task_ivar_consensus.wdl" as consensus_task
import "../../tasks/gene_typing/variant_detection/task_ivar_variant_call.wdl" as variant_call_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics
import "../../tasks/quality_control/read_filtering/task_ivar_primer_trim.wdl" as primer_trim_task

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
    Boolean? skip_N
    Int? ivar_bwa_cpu
    Int? ivar_bwa_memory
    Int? ivar_bwa_disk_size
    String? ivar_bwa_docker
    Int? ivar_trim_primers_cpu
    Int? ivar_trim_primers_memory
    Int? ivar_trim_primers_disk_size
    String? ivar_trim_primers_docker
    Int? stats_n_coverage_primtrim_cpu
    Int? stats_n_coverage_primtrim_memory
    Int? stats_n_coverage_primtrim_disk_size
    String? stats_n_coverage_primtrim_docker
    Int? ivar_variant_cpu
    Int? ivar_variant_memory
    Int? ivar_variant_disk_size
    String? ivar_variant_docker
    Int? ivar_consensus_cpu
    Int? ivar_consensus_memory
    Int? ivar_consensus_disk_size
    String? ivar_consensus_docker
    Int? stats_n_coverage_cpu
    Int? stats_n_coverage_memory
    Int? stats_n_coverage_disk_size
    String? stats_n_coverage_docker
  }
  call bwa_task.bwa {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      reference_genome = reference_genome,
      cpu = ivar_bwa_cpu,
      memory = ivar_bwa_memory,
      disk_size = ivar_bwa_disk_size,
      docker = ivar_bwa_docker
  }
  if (trim_primers) {
    call primer_trim_task.primer_trim {
      input:
        samplename = samplename,
        primer_bed = select_first([primer_bed]),
        bamfile = bwa.sorted_bam,
        cpu = ivar_trim_primers_cpu,
        memory = ivar_trim_primers_memory,
        disk_size = ivar_trim_primers_disk_size,
        docker = ivar_trim_primers_docker
    }
    call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
    input:
      samplename = samplename,
      bamfile = primer_trim.trim_sorted_bam,
      read1 = read1,
      read2 = read2,
      cpu = stats_n_coverage_primtrim_cpu,
      memory = stats_n_coverage_primtrim_memory,
      disk_size = stats_n_coverage_primtrim_disk_size,
      docker = stats_n_coverage_primtrim_docker
    }
  }
  call consensus_task.consensus {
    input:
      samplename = samplename,
      bamfile = select_first([primer_trim.trim_sorted_bam, bwa.sorted_bam]),
      reference_genome = reference_genome,
      consensus_min_depth = min_depth,
      consensus_min_freq = consensus_min_freq,
      skip_N = skip_N,
      cpu = ivar_consensus_cpu,
      memory = ivar_consensus_memory,
      disk_size = ivar_consensus_disk_size,
      docker = ivar_consensus_docker
  }
  call variant_call_task.variant_call {
    input:
      samplename = samplename,
      mpileup = consensus.sample_mpileup,
      reference_gff = reference_gff,
      reference_genome = reference_genome,
      variant_min_depth = min_depth,
      variant_min_freq = variant_min_freq,
      cpu = ivar_variant_cpu,
      memory = ivar_variant_memory,
      disk_size = ivar_variant_disk_size,
      docker = ivar_variant_docker
  }
  call assembly_metrics.stats_n_coverage {
    input:
      samplename = samplename,
      bamfile = bwa.sorted_bam,
      cpu = stats_n_coverage_cpu,
      memory = stats_n_coverage_memory,
      disk_size = stats_n_coverage_disk_size,
      docker = stats_n_coverage_docker
  }
  output {
    # bwa outputs
    String bwa_version = bwa.bwa_version
    String samtools_version = bwa.sam_version
    File read1_aligned = bwa.read1_aligned
    File? read2_aligned = bwa.read2_aligned
    File aligned_bam =  select_first([primer_trim.trim_sorted_bam, bwa.sorted_bam, ""]) # set default values for select_first() to avoid workflow failures
    File aligned_bai = select_first([primer_trim.trim_sorted_bai, bwa.sorted_bai, ""])
    File read1_unaligned = bwa.read1_unaligned
    File? read2_unaligned = bwa.read2_unaligned
    File sorted_bam_unaligned = bwa.sorted_bam_unaligned
    File sorted_bam_unaligned_bai = bwa.sorted_bam_unaligned_bai
    # primer trimming outputs
    Float? primer_trimmed_read_percent = primer_trim.primer_trimmed_read_percent
    String? ivar_version_primtrim = primer_trim.ivar_version
    String? samtools_version_primtrim = primer_trim.samtools_version
    String? primer_bed_name = primer_trim.primer_bed_name
    # variant call outputs
    File ivar_tsv = variant_call.sample_variants_tsv
    File ivar_vcf = variant_call.sample_variants_vcf
    String ivar_variant_proportion_intermediate = variant_call.variant_proportion_intermediate
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
    # Assembly metrics
    String percentage_mapped_reads = select_first([stats_n_coverage_primtrim.percentage_mapped_reads, stats_n_coverage.percentage_mapped_reads,""])
  }
}
