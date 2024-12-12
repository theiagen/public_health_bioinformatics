version 1.0

import "../../tasks/quality_control/read_filtering/task_porechop.wdl" as task_porechop
import "../../tasks/assembly/task_flye.wdl" as task_flye
import "../../tasks/assembly/task_bandageplot.wdl" as task_bandage
import "../../tasks/polishing/task_medaka.wdl" as task_medaka
import "../../tasks/polishing/task_racon.wdl" as task_racon
import "../../tasks/assembly/task_dnaapler.wdl" as task_dnaapler
import "../../tasks/task_versioning.wdl" as versioning_task
import "../../tasks/quality_control/read_filtering/task_filtercontigs.wdl" as task_filtercontigs
import "../../tasks/alignment/task_bwa.wdl" as task_bwaall
import "../../tasks/polishing/task_polypolish.wdl" as task_polypolish

workflow flye_denovo {
  meta {
    description: "This workflow assembles long-read sequencing data using Flye, optionally trims reads with Porechop, and generates an assembly graph visualization with Bandage. It supports consensus polishing with Medaka or Racon for long reads, or Polypolish for hybrid assemblies with Illumina short reads. The workflow concludes by reorienting contigs with Dnaapler for a final assembly."
  }
  input {
    File read1                    # raw input reads
    File? illumina_read1          # Optional Illumina short-read R1 for hybrid assembly
    File? illumina_read2          # Optional Illumina short-read R2 for hybrid assembly
    String samplename
    String polisher = "medaka"
    Int polish_rounds = 1        # Optional number of polishing rounds
    
    # Task-level inputs for Porechop
    Int porechop_cpu = 4
    Int porechop_memory = 8
    Int porechop_disk_size = 50
    String? porechop_trimopts     # Optional Porechop trimming options

    # Task-level parameters for Flye
    Boolean ont_corrected = false
    Boolean ont_high_quality = false
    Boolean pacbio_raw = false
    Boolean pacbio_corrected = false
    Boolean pacbio_hifi = false
    Int? genome_length              # Requires `asm_coverage`
    Int? asm_coverage               # Reduced coverage for initial disjointig assembly
    Int flye_polishing_iterations = 1    # Default polishing iterations
    Int? minimum_overlap            # Minimum overlap between reads
    Float? read_error_rate          # Maximum expected read error rate
    Boolean uneven_coverage_mode = false
    Boolean keep_haplotypes = false
    Boolean no_alt_contigs = false
    Boolean scaffold = false
    String? additional_parameters   # Any extra Flye-specific parameters
    Int flye_cpu = 4
    Int flye_memory = 32
    Int flye_disk_size = 100

    # Bandage-specific inputs
    Int bandage_cpu = 2
    Int bandage_memory = 4
    Int bandage_disk_size = 10

    #Task-level parameters for Polypolish
    # Polypolish-specific inputs
    String? pair_orientation
    Float? low_percentile_threshold
    Float? high_percentile_threshold
    Float? fraction_invalid
    Float? fraction_valid
    Int? maximum_errors
    Int? minimum_depth
    Boolean careful = false
    Int polypolish_cpu = 1
    Int polypolish_memory = 8
    Int polypolish_disk_size = 100

    # Medaka-specific inputs
    Boolean auto_medaka_model = true  # Enable automatic Medaka model selection
    String? medaka_model_override # Optional user-specified Medaka model
    Int medaka_cpu = 4
    Int medaka_memory = 16
    Int medaka_disk_size = 100

     # Racon-specific inputs
    Int racon_cpu = 8
    Int racon_memory = 16
    Int racon_disk_size = 100

    # Contig filter-specific inputs
    Int contig_filter_min_len = 1000
    Boolean contig_filter_homopolymers = true
    Int contig_filter_cpu = 4
    Int contig_filter_memory = 16
    Int contig_filter_disk_size = 100

    # Dnaapler-specific inputs
    String dnaapler_mode = "all"
    Int dnaapler_cpu = 4
    Int dnaapler_memory = 16
    Int dnaapler_disk_size = 100

    #Other workdlow-level inputs
    Boolean skip_trim_reads = true   # Default: No trimming
    Boolean skip_polishing = false    # Default: Polishing enabled
  }
  call versioning_task.version_capture {
    input:
  }
  # Optional Porechop trimming before Flye
  if (!skip_trim_reads) {
    call task_porechop.porechop as porechop {
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
  call task_flye.flye as flye {
    input:
      read1 = select_first([porechop.trimmed_reads, read1]),  # Use trimmed reads if available
      samplename = samplename,
      ont_corrected = ont_corrected,
      pacbio_corrected = pacbio_corrected,
      pacbio_hifi = pacbio_hifi,
      pacbio_raw = pacbio_raw,
      ont_high_quality = ont_high_quality,
      genome_length = genome_length,
      asm_coverage = asm_coverage,
      flye_polishing_iterations = flye_polishing_iterations,
      minimum_overlap = minimum_overlap,
      read_error_rate = read_error_rate,
      uneven_coverage_mode = uneven_coverage_mode,
      keep_haplotypes = keep_haplotypes,
      no_alt_contigs = no_alt_contigs,
      scaffold = scaffold,
      additional_parameters = additional_parameters,
      cpu = flye_cpu,
      memory = flye_memory,
      disk_size = flye_disk_size
  }
  # Generate Bandage plot
  call task_bandage.bandage_plot as bandage {
    input:
      assembly_graph_gfa = flye.assembly_graph_gfa,
      samplename = samplename,
      cpu = bandage_cpu,
      memory = bandage_memory,
      disk_size = bandage_disk_size
  }
  # Hybrid Assembly Path: Polypolish
  if (defined(illumina_read1) && defined(illumina_read2)) {
   call task_bwaall.bwa_all as bwa {
     input:
       draft_assembly_fasta = flye.assembly_fasta,
       read1 = select_first([illumina_read1]),
       read2 = select_first([illumina_read2]),
       samplename = samplename
    }
    call task_polypolish.polypolish as polypolish {
      input:
        assembly_fasta = flye.assembly_fasta,
        read1_sam = bwa.read1_sam,
        read2_sam = bwa.read2_sam,
        samplename = samplename,
        illumina_polishing_rounds = polish_rounds,
        pair_orientation = pair_orientation,
        low_percentile_threshold = low_percentile_threshold,
        high_percentile_threshold = high_percentile_threshold,
        fraction_invalid = fraction_invalid,
        fraction_valid = fraction_valid,
        maximum_errors = maximum_errors,
        minimum_depth = minimum_depth,
        careful = careful,
        cpu = polypolish_cpu,
        memory = polypolish_memory,
        disk_size = polypolish_disk_size
    }
  }
  # ONT-only Polishing Path: Medaka or Racon
  if (!defined(illumina_read1) || !defined(illumina_read2)) {
    if (!skip_polishing && polisher == "medaka") {
      call task_medaka.medaka_consensus as medaka {
        input:
          unpolished_fasta = flye.assembly_fasta,
          samplename = samplename,
          read1 = select_first([porechop.trimmed_reads, read1]),  # Use trimmed reads if available
          medaka_model_override = medaka_model_override,
          auto_model = auto_medaka_model
      }
    }
    if (!skip_polishing && polisher == "racon") {
      call task_racon.racon as racon {
        input:
          unpolished_fasta = flye.assembly_fasta,
          read1 = select_first([porechop.trimmed_reads, read1]),  # Use trimmed reads if available
          samplename = samplename,
          polishing_rounds = polish_rounds,
          cpu = racon_cpu,
          memory = racon_memory,
          disk_size = racon_disk_size
      }
    }
  }
  # Contig Filtering and Final Assembly orientation
  call task_filtercontigs.contig_filter as contig_filter {
    input:
      assembly_fasta = select_first([polypolish.polished_assembly, medaka.medaka_fasta, racon.polished_fasta, flye.assembly_fasta]), # Use Flye assembly if no polishing
      min_len = contig_filter_min_len,
      filter_homopolymers = contig_filter_homopolymers,
      threads = contig_filter_cpu,
      memory = contig_filter_memory,
      disk_size = contig_filter_disk_size
  }
  call task_dnaapler.dnaapler_all as dnaapler {
    input:
      input_fasta = contig_filter.filtered_fasta,
      samplename = samplename,
      dmaapler_mode = dnaapler_mode,
      cpu = dnaapler_cpu,
      memory = dnaapler_memory,
      disk_size = dnaapler_disk_size
  }
  output {
    File final_assembly = dnaapler.reoriented_fasta
    #add into terra export table
    File bandage_plot = bandage.plot
    File contigs_gfa = flye.assembly_graph_gfa
    String flye_phb_version = version_capture.phb_version
    String flye_analysis_date = version_capture.date
  }
}