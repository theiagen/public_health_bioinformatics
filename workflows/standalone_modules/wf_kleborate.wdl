version 1.0

import "../../tasks/species_typing/task_kleborate.wdl" as kleborate_task
import "../../tasks/task_versioning.wdl" as versioning

workflow kleborate_wf {
  input {
      File assembly
      String samplename
    }
  call kleborate_task.kleborate {
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
    
    File kleborate_report = kleborate.kleborate_output_file
    String kleborate_version = kleborate.kleborate_version
    String kleborate_docker = kleborate.kleborate_docker
    String kleborate_mlst_sequence_type = kleborate.kleborate_mlst_sequence_type
    String kleborate_virulence_score = kleborate.kleborate_virulence_score
    String kleborate_resistance_score = kleborate.kleborate_resistance_score
    String kleborate_num_resistance_genes = kleborate.kleborate_num_resistance_genes
    String kleborate_bla_resistance_genes = kleborate.kleborate_bla_resistance_genes
    String kleborate_esbl_resistance_genes = kleborate.kleborate_esbl_resistance_genes
    String kleborate_key_resistance_genes = kleborate.kleborate_key_resistance_genes
    String kleborate_resistance_mutations = kleborate.kleborate_genomic_resistance_mutations
    String kleborate_klocus = kleborate.kleborate_klocus
    String kleborate_ktype = kleborate.kleborate_ktype
    String kleborate_olocus = kleborate.kleborate_olocus
    String kleborate_otype = kleborate.kleborate_otype
    String kleborate_klocus_confidence = kleborate.kleborate_klocus_confidence
    String kleborate_olocus_confidence = kleborate.kleborate_olocus_confidence
    }
 }
