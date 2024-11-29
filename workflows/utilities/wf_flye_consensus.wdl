version 1.0

import "../../tasks/assembly/task_flye.wdl" as task_flye
import "../../tasks/assembly/task_bandageplot.wdl" as task_bandage
import "../../tasks/assembly/task_medaka.wdl" as task_medaka
import "../../tasks/assembly/task_racon.wdl" as task_racon
import "../../tasks/assembly/task_dnaapler.wdl" as task_dnaapler

workflow flye_consensus {
  meta {
    description: "Run Flye assembly, Bandage plot, and consensus polishing with Medaka or Racon from long reads"
  }

  input {
    File read1
    String samplename
    String medaka_model
    String polisher = "medaka"
    Int polish_rounds
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
  call task_dnaapler.dnaapler_all as dnaapler {
    input:
      input_fasta = select_first([medaka.medaka_fasta, racon.polished_fasta]),
      prefix = samplename,
      threads = polish_rounds
  }

  output {
    File final_assembly = dnaapler.reoriented_fasta
    File bandage_plot = bandage.plot
    String flye_version = flye.flye_version
    String bandage_version = bandage.bandage_version
    String? medaka_version = medaka.medaka_version
    String? racon_version = racon.racon_version
  }
}
