version 1.0


import "../tasks/task_taxon_id.wdl" as taxon
import "../tasks/task_versioning.wdl" as versioning

workflow kleborate_wf {
  input {
      File assembly
      String samplename
    }
  call taxon.kleborate_one_sample {
    input:
      assembly = assembly,
      samplename = samplename
    }
  call versioning.version_capture{
    input:
  }
  output {
    String kleborate_wf_version = version_capture.phbg_version
    String kleborate_wf_analysis_date = version_capture.date
    
    File kleborate_report = kleborate_one_sample.kleborate_output_file
    String kleborate_version = kleborate_one_sample.version
    String kleborate_mlst_sequence_type = kleborate_one_sample.mlst_sequence_type
    String kleborate_virulence_score = kleborate_one_sample.virulence_score
    String kleborate_resistance_score = kleborate_one_sample.resistance_score
    String kleborate_num_resistance_genes = kleborate_one_sample.num_resistance_genes
    String kleborate_bla_resistance_genes = kleborate_one_sample.bla_resistance_genes
    String kleborate_esbl_resistance_genes = kleborate_one_sample.esbl_resistance_genes
    String kleborate_key_resistance_genes = kleborate_one_sample.key_resistance_genes
    String kleborate_resistance_mutations = kleborate_one_sample.resistance_mutations
    }
 }
