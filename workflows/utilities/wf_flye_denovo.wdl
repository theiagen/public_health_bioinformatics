version 1.0

import "../../tasks/quality_control/read_filtering/task_porechop.wdl" as task_porechop
import "../../tasks/assembly/task_flye.wdl" as task_flye
import "../../tasks/assembly/task_bandage_plot.wdl" as task_bandage
import "../../tasks/polishing/task_medaka.wdl" as task_medaka
import "../../tasks/polishing/task_racon.wdl" as task_racon
import "../../tasks/assembly/task_dnaapler.wdl" as task_dnaapler
import "../../tasks/quality_control/read_filtering/task_filter_contigs.wdl" as task_filter_contigs
import "../../tasks/alignment/task_bwa.wdl" as task_bwa_all
import "../../tasks/polishing/task_polypolish.wdl" as task_polypolish

workflow flye_denovo {
  meta {
    description: "This workflow assembles long-read sequencing data using Flye, optionally trims reads with Porechop, and generates an assembly graph visualization with Bandage. It supports consensus polishing with Medaka or Racon for long reads, or Polypolish for hybrid assemblies with Illumina short reads. The workflow concludes by reorienting contigs with Dnaapler for a final assembly."
  }
  input {
    File read1 
    File? illumina_read1                   
    File? illumina_read2                  
    String samplename
    String polisher = "medaka"
    Int polish_rounds = 1                   
    Boolean skip_trim_reads = true # Default: No trimming
    Boolean skip_polishing = false # Default: Polishing enabled
    
    # Porechop inputs
    Int porechop_cpu = 4
    Int porechop_memory = 8
    Int porechop_disk_size = 50
    String? porechop_trimopts # Optional Porechop trimming options

    # Flye inputs
    String flye_read_type = "--nano-hq" 
    Int? flye_genome_length # Requires `asm_coverage`
    Int? flye_asm_coverage # Reduced coverage for initial disjointig assembly
    Int flye_polishing_iterations = 1 # Default polishing iterations
    Int? flye_minimum_overlap # Minimum overlap between reads
    Float? flye_read_error_rate # Maximum expected read error rate
    Boolean flye_uneven_coverage_mode = false
    Boolean flye_keep_haplotypes = false
    Boolean flye_no_alt_contigs = false
    Boolean flye_scaffold = false
    String? flye_additional_parameters # Any extra Flye-specific parameters
    Int flye_cpu = 4
    Int flye_memory = 32
    Int flye_disk_size = 100

    # Bandage inputs
    Int bandage_cpu = 2
    Int bandage_memory = 4
    Int bandage_disk_size = 10

    # Polypolish inputs
    String? polypolish_pair_orientation
    Float? polypolish_low_percentile_threshold
    Float? polypolish_high_percentile_threshold
    Float? polypolish_fraction_invalid
    Float? polypolish_fraction_valid
    Int? polypolish_maximum_errors
    Int? polypolish_minimum_depth
    Boolean polypolish_careful = false
    Int polypolish_cpu = 1
    Int polypolish_memory = 8
    Int polypolish_disk_size = 100

    # Medaka inputs
    Boolean auto_medaka_model = true # Enable automatic Medaka model selection
    String? medaka_model # Optional user-specified Medaka model
    Int medaka_cpu = 4
    Int medaka_memory = 16
    Int medaka_disk_size = 100

     # Racon inputs
    Int racon_cpu = 8
    Int racon_memory = 16
    Int racon_disk_size = 100

    # Contig filter-specific inputs
    Int filter_contigs_min_length = 1000
    Int filter_contigs_cpu = 4
    Int filter_contigs_memory = 16
    Int filter_contigs_disk_size = 100

    # Dnaapler inputs
    String dnaapler_mode = "all"
    Int dnaapler_cpu = 4
    Int dnaapler_memory = 16
    Int dnaapler_disk_size = 100
  }
  # Optional Porechop trimming before Flye
  if (!skip_trim_reads) {
    call task_porechop.porechop {
      input:
        read1 = read1,
        samplename = samplename,
        trimopts = porechop_trimopts,
        cpu = porechop_cpu,
        memory = porechop_memory,
        disk_size = porechop_disk_size
    }
  }
  # Call Flye using either trimmed reads or raw reads
  call task_flye.flye {
    input:
      read1 = select_first([porechop.trimmed_reads, read1]), # Use trimmed reads if available
      samplename = samplename,
      genome_length = flye_genome_length,
      asm_coverage = flye_asm_coverage,
      flye_polishing_iterations = flye_polishing_iterations,
      minimum_overlap = flye_minimum_overlap,
      read_error_rate = flye_read_error_rate,
      uneven_coverage_mode = flye_uneven_coverage_mode,
      keep_haplotypes = flye_keep_haplotypes,
      no_alt_contigs = flye_no_alt_contigs,
      scaffold = flye_scaffold,
      additional_parameters = flye_additional_parameters,
      cpu = flye_cpu,
      memory = flye_memory,
      disk_size = flye_disk_size
  }
  # Bandage plot generation
  call task_bandage.bandage_plot as bandage {
    input:
      assembly_graph_gfa = flye.assembly_graph_gfa,
      samplename = samplename,
      cpu = bandage_cpu,
      memory = bandage_memory,
      disk_size = bandage_disk_size
  }
  # Polypolish for hybrid assembly
  if (defined(illumina_read1) && defined(illumina_read2)) {
   call task_bwa_all.bwa_all as bwa {
     input:
       draft_assembly_fasta = flye.assembly_fasta,
       read1 = select_first([illumina_read1]),
       read2 = select_first([illumina_read2]),
       samplename = samplename
    }
    call task_polypolish.polypolish {
      input:
        assembly_fasta = flye.assembly_fasta,
        read1_sam = bwa.read1_sam,
        read2_sam = bwa.read2_sam,
        samplename = samplename,
        illumina_polishing_rounds = polish_rounds,
        pair_orientation = polypolish_pair_orientation,
        low_percentile_threshold = polypolish_low_percentile_threshold,
        high_percentile_threshold = polypolish_high_percentile_threshold,
        fraction_invalid = polypolish_fraction_invalid,
        fraction_valid = polypolish_fraction_valid,
        maximum_errors = polypolish_maximum_errors,
        minimum_depth = polypolish_minimum_depth,
        careful = polypolish_careful,
        cpu = polypolish_cpu,
        memory = polypolish_memory,
        disk_size = polypolish_disk_size
    }
  }
  # ONT-only Polishing Path: Medaka or Racon
  if (!skip_polishing) {
    if (polisher == "medaka") {
      call task_medaka.medaka {
        input:
          unpolished_fasta = flye.assembly_fasta,
          samplename = samplename,
          read1 = select_first([porechop.trimmed_reads, read1]),
          medaka_model = medaka_model,
          auto_model = auto_medaka_model,
          cpu = medaka_cpu,
          memory = medaka_memory,
          disk_size = medaka_disk_size
      }
    }
    if (polisher == "racon") {
      call task_racon.racon {
        input:
          unpolished_fasta = flye.assembly_fasta,
          read1 = select_first([porechop.trimmed_reads, read1]),
          samplename = samplename,
          polishing_rounds = polish_rounds,
          cpu = racon_cpu,
          memory = racon_memory,
          disk_size = racon_disk_size
      }
    }
  }
  # Contig Filtering and Final Assembly orientation
  call task_filter_contigs.filter_contigs {
    input:
      assembly_fasta = select_first([polypolish.polished_assembly, medaka.medaka_fasta, racon.polished_fasta, flye.assembly_fasta]), # Use Flye assembly if no polishing
      min_length = filter_contigs_min_length,
      cpu = filter_contigs_cpu,
      memory = filter_contigs_memory,
      disk_size = filter_contigs_disk_size
  }
  call task_dnaapler.dnaapler {
    input:
      input_fasta = filter_contigs.filtered_fasta,   
      samplename = samplename,
      dnaapler_mode = dnaapler_mode,
      cpu = dnaapler_cpu,
      memory = dnaapler_memory,
      disk_size = dnaapler_disk_size
  }
  output { 
    File assembly_fasta = dnaapler.reoriented_fasta
    File bandage_plot = bandage.plot
    File contigs_gfa = flye.assembly_graph_gfa
    File filtered_contigs_metrics = filter_contigs.assembly_filtering_metrics
    File flye_assembly_info = flye.assembly_info
    String? medaka_model_used = medaka.resolved_medaka_model
    String? porechop_version = porechop.porechop_version
    String flye_version = flye.flye_version
    String bandage_version = bandage.bandage_version
    String? medaka_version = medaka.medaka_version
    String? racon_version = racon.racon_version
    String? bwa_version = bwa.bwa_version
    String? polypolish_version = polypolish.polypolish_version
    String dnaapler_version = dnaapler.dnaapler_version
  }
}