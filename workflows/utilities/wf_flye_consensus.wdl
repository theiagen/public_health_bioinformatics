version 1.0

import "../../tasks/assembly/task_flye.wdl" as task_flye
import "../../tasks/assembly/task_bandageplot.wdl" as task_bandage
import "../../tasks/assembly/task_medaka.wdl" as task_medaka
import "../../tasks/assembly/task_racon.wdl" as task_racon
import "../../tasks/assembly/task_dnaapler.wdl" as task_dnaapler
import "../../tasks/task_versioning.wdl" as versioning_task
import "../../tasks/assembly/task_filtercontigs.wdl" as task_filtercontigs

workflow flye_consensus {
  meta {
    description: "Run Flye assembly, Bandage plot, and consensus polishing with Medaka or Racon from long reads"
  }
  input {
    File read1
    String samplename
    String polisher = "medaka"
    Int? polish_rounds
    String? medaka_model
  }
  call versioning_task.version_capture {
    input:
  }
  call task_flye.flye as flye {
    input:
      read1 = read1,
      samplename = samplename
  }
  call task_bandage.bandage_plot as bandage {
    input:
      assembly_graph_gfa = flye.assembly_graph,
      samplename = samplename
  }
  if (polisher == "medaka") {
    call task_medaka.medaka_consensus as medaka {
      input:
        assembly_fasta = flye.assembly_fasta,
        samplename = samplename,
        read1 = read1,
        medaka_model = medaka_model,
        polish_rounds = polish_rounds
    }
  }
  if (polisher == "racon") {
    call task_racon.racon as racon {
      input:
        unpolished_fasta = flye.assembly_fasta,
        reads = read1,
        samplename = samplename
    }
  }
  call task_filtercontigs.contig_filter as contig_filter {
    input:
      assembly_fasta = select_first([medaka.medaka_fasta, racon.polished_fasta])
  }
  call task_dnaapler.dnaapler_all as dnaapler {
    input:
      input_fasta = contig_filter.filtered_fasta,
      samplename = samplename
  }
  output {
    File final_assembly = dnaapler.reoriented_fasta
    File bandage_plot = bandage.plot
    File filtered_assembly = contig_filter.filtered_fasta
    String flye_phb_version = version_capture.phb_version    
    String flye_analysis_date = version_capture.date
  }
}