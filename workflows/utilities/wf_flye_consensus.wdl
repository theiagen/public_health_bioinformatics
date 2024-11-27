version 1.0

import "../../tasks/assembly/task_flye.wdl" as task_flye
import "../../tasks/assembly/task_bandageplot.wdl" as task_bandage
import "../../tasks/assembly/task_medaka.wdl" as task_medaka

workflow flye_consensus {
  # Input parameters
  input {
    File read1
    String samplename
    String medaka_model
    Int polish_rounds
  }
  # Step 1: Run Flye assembly
  call task_flye.flye as flye {
    input:
      read1 = read1,
      samplename = samplename,
  }
  # Step 2: Generate Bandage plot from Flye assembly graph
  call task_bandage.bandage_plot as bandage {
    input:
      assembly_graph_gfa = flye.assembly_graph,
      samplename = samplename
  }
  # Step 3: Perform Medaka polishing on Flye assembly
  call task_medaka.medaka_consensus as medaka {
    input:
      assembly_fasta = flye.assembly_fasta,
      samplename = samplename,
      read1 = read1,
      medaka_model = medaka_model,
      polish_rounds = polish_rounds,
  }
  # Outputs
  output {
    File final_assembly = medaka.medaka_fasta
    File bandage_plot = bandage.plot
    File medaka_vcf = medaka.medaka_vcf
    String flye_version = flye.flye_version
    String bandage_version = bandage.bandage_version
    String medaka_version = medaka.medaka_version
  }
}

