version 1.0

import "../../tasks/assembly/task_racon.wdl" as task_racon
import "../../tasks/assembly/task_minimap2.wdl" as task_minimap2

workflow racon_polishing_workflow {
  input {
    File assembly_fasta          # Initial unpolished assembly
    File read1                   # Input sequencing reads
    String samplename            # Base name for output files
    Int racon_rounds = 1         # Number of Racon polishing rounds
  }

  # Step 1: Run Minimap2 to produce alignments
  call minimap2 {
    input:
      query1 = read1,
      reference = assembly_fasta,
      samplename = samplename + "_alignments",
  }

  # Step 2: Iteratively run Racon for polishing
  scatter (round in range(racon_rounds)) {
    if (round == 0) {
      # First round uses the initial assembly and Minimap2 alignments
      call racon {
        input:
          unpolished_fasta = assembly_fasta,
          reads = read1,
          alignments = minimap2.minimap2_out,
          samplename = samplename + "_racon_round" + round,
      }
    } else {
      # Subsequent rounds use the output from the previous Racon run
      call racon {
        input:
          unpolished_fasta = racon.polished_fasta,
          reads = read1,
          alignments = minimap2.minimap2_out,
          samplename = samplename + "_racon_round" + round,
      }
    }
  }

  # Output the final polished assembly
  output {
    File final_polished_fasta = racon.polished_fasta
  }
}
