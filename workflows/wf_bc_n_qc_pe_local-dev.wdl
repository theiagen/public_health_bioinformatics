version 1.0

import "wf_bc_n_qc_pe.wdl" as assembly

workflow bc_n_qc_local {
  input {
    Array[Pair[Array[String], Pair[File,File]]] inputSamples
  }

  scatter (sample in inputSamples) {
    call assembly.bc_n_qc_pe {
      input:
        samplename = sample.left[0],
        read1_raw  = sample.right.left,
        read2_raw  = sample.right.right
    }
  }

}
