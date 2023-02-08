version 1.0

import "../../tasks/species_typing/task_kleborate.wdl" as kleborate
import "../../tasks/task_versioning.wdl" as versioning

workflow kleborate_wf {
  input {
      File assembly
      String samplename
    }
  call kleborate.kleborate_standalone {
    input:
      assembly = assembly,
      samplename = samplename
    }
  call versioning.version_capture{
    input:
  }
  output {
    String kleborate_wf_version = version_capture.phb_version
    String kleborate_wf_analysis_date = version_capture.date
    
    File kleborate_report = kleborate_standalone.kleborate_output_file
    String kleborate_version = kleborate_standalone.version
    String kleborate_mlst_sequence_type = kleborate_standalone.mlst_sequence_type
    String kleborate_virulence_score = kleborate_standalone.virulence_score
    String kleborate_resistance_score = kleborate_standalone.resistance_score
    String kleborate_num_resistance_genes = kleborate_standalone.num_resistance_genes
    String kleborate_bla_resistance_genes = kleborate_standalone.bla_resistance_genes
    String kleborate_esbl_resistance_genes = kleborate_standalone.esbl_resistance_genes
    String kleborate_key_resistance_genes = kleborate_standalone.key_resistance_genes
    String kleborate_resistance_mutations = kleborate_standalone.resistance_mutations
    }
 }
