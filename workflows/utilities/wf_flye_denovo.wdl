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
    Int? polish_rounds = 1        # Optional number of polishing rounds
    String? medaka_model_override # Optional user-specified Medaka model
    Boolean auto_medaka_model = true  # Enable automatic Medaka model selection
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
        samplename = samplename
    }
  }
  # Call Flye using either trimmed reads or raw reads
  call task_flye.flye as flye {
    input:
      read1 = select_first([porechop.trimmed_reads, read1]),  # Use trimmed reads if available
      samplename = samplename
  }
  # Generate Bandage plot
  call task_bandage.bandage_plot as bandage {
    input:
      assembly_graph_gfa = flye.assembly_graph_gfa,
      samplename = samplename
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
        illumina_polishing_rounds = polish_rounds
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
          polishing_rounds = polish_rounds
      }
    }
  }
  # Contig Filtering and Final Assembly orientation
  call task_filtercontigs.contig_filter as contig_filter {
    input:
      assembly_fasta = select_first([polypolish.polished_assembly, medaka.medaka_fasta, racon.polished_fasta, flye.assembly_fasta]), # Use Flye assembly if no polishing
  }
  call task_dnaapler.dnaapler_all as dnaapler {
    input:
      input_fasta = contig_filter.filtered_fasta,
      samplename = samplename
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