version 1.0

import "../../tasks/gene_typing/drug_resistance/task_chroquetas.wdl" as chroquetas_task
import "../../tasks/task_versioning.wdl" as versioning

workflow chroquetas_workflow {
  meta {
    description: "ChroQueTas fungal AMR detection referencing the FungAMR database"
  }
  input {
    File assembly_fasta
    String species
    String samplename
  }
  call versioning.version_capture {
    input:
  }
  call chroquetas_task.chroquetas {
    input:
      assembly_fasta = assembly_fasta,
      samplename = samplename,
      species = species
  }
  output {
    # ChroQueTas outputs
    File? chroquetas_amr_stats_file = chroquetas.amr_stats_file
    File? chroquetas_amr_summary_file = chroquetas.amr_summary_file
    String chroquetas_fungicide_resistance = chroquetas.chroquetas_fungicide_resistance
    String chroquetas_status = chroquetas.chroquetas_status
    String chroquetas_version = chroquetas.chroquetas_version

    # PHB Versioning
    String chroquetas_wf_analysis_date = version_capture.date
    String chroquetas_wf_version = version_capture.phb_version
  }
}
