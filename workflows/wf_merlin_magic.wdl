version 1.0

import "../tasks/species_typing/task_serotypefinder.wdl" as serotypefinder
import "../tasks/species_typing/task_ectyper.wdl" as ectyper
import "../tasks/species_typing/task_lissero.wdl" as lissero
import "../tasks/species_typing/task_sistr.wdl" as sistr
import "../tasks/species_typing/task_seqsero2.wdl" as seqsero2
import "../tasks/species_typing/task_kleborate.wdl" as kleborate
import "../tasks/species_typing/task_tbprofiler.wdl" as tbprofiler
import "../tasks/species_typing/task_legsta.wdl" as legsta
import "../tasks/species_typing/task_genotyphi.wdl" as genotyphi
import "../tasks/species_typing/task_kaptive.wdl" as kaptive
import "../tasks/species_typing/task_seroba.wdl" as seroba
import "../tasks/species_typing/task_pbptyper.wdl" as pbptyper
import "../tasks/species_typing/task_poppunk_streppneumo.wdl" as poppunk_spneumo
import "../tasks/gene_typing/task_abricate.wdl" as abricate_task

workflow merlin_magic {
  meta {
    description: "Workflow for bacterial species typing; based on the Bactopia subworkflow Merlin (https://bactopia.github.io/bactopia-tools/merlin/)"
  }
  input {
    String samplename
    String merlin_tag
    File assembly
    File read1
    File? read2
    Boolean paired_end = true
    Boolean call_poppunk = true
  }
    if (merlin_tag == "Acinetobacter baumannii") {
    call kaptive.kaptive {
      input:
        assembly = assembly,
        samplename = samplename
    }
    call abricate_task.abricate {
      input:
        assembly = assembly,
        samplename = samplename,
        database = "AcinetobacterPlasmidTyping",
        minid = 95 # strict threshold of 95% identity for typing purposes
    }
  }
  if (merlin_tag == "Escherichia") {
    call serotypefinder.serotypefinder {
      input:
        assembly = assembly,
        samplename = samplename
    }
    call ectyper.ectyper {
      input:
        assembly = assembly,
        samplename = samplename
    }
  }
  if (merlin_tag == "Listeria") {
    call lissero.lissero {
      input:
        assembly = assembly,
        samplename = samplename
    }
  }
  if (merlin_tag == "Salmonella") {
    call sistr.sistr {
      input: 
        assembly = assembly,
        samplename = samplename
    }
    call seqsero2.seqsero2 as seqsero2 {
      input: 
        read1 = read1,
        read2 = read2,
        samplename = samplename,
        paired_end = paired_end
    }
    if( seqsero2.seqsero2_predicted_serotype == "Typhi" || sistr.sistr_predicted_serotype == "Typhi" ) {
      call genotyphi.genotyphi as genotyphi_task {
        input: 
          read1 = read1,
          read2 = read2,
          samplename = samplename
      }
    }
  }
  if (merlin_tag == "Klebsiella") {
    call kleborate.kleborate {
      input:
        assembly = assembly,
        samplename = samplename
    }
  }
  if (merlin_tag == "Mycobacterium tuberculosis") {
    call tbprofiler.tbprofiler {
      input:
        read1 = read1,
        read2 = read2,
        samplename = samplename
    }
  }
  if (merlin_tag == "Legionella pneumophila") {
    call legsta.legsta {
      input:
        assembly = assembly,
        samplename = samplename
    }
  }
  if (merlin_tag == "Streptococcus pneumoniae") {
    if (paired_end) {
      call seroba.seroba as seroba_task {
        input:
          read1 = read1,
          read2 = read2,
          samplename = samplename
      }
    }
    call pbptyper.pbptyper as pbptyper_task {
      input:
        assembly = assembly,
        samplename = samplename
    }      
    if (call_poppunk) {
      call poppunk_spneumo.poppunk as poppunk_task {
        input:
          assembly = assembly,
          samplename = samplename
      }  
    }
  }

  output {
  # Ecoli Typing
  File? serotypefinder_report = serotypefinder.serotypefinder_report
  String? serotypefinder_docker = serotypefinder.serotypefinder_docker
  String? serotypefinder_serotype = serotypefinder.serotypefinder_serotype
  File? ectyper_results = ectyper.ectyper_results
  String? ectyper_version = ectyper.ectyper_version
  String? ectyper_predicted_serotype = ectyper.ectyper_predicted_serotype
  # Listeria Typing
  File? lissero_results = lissero.lissero_results
  String? lissero_version = lissero.lissero_version
  String? lissero_serotype = lissero.lissero_serotype
  # Salmonella Typing
  File? sistr_results = sistr.sistr_results
  File? sistr_allele_json = sistr.sistr_allele_json
  File? sistr_allele_fasta = sistr.sistr_allele_fasta
  File? sistr_cgmlst = sistr.sistr_cgmlst
  String? sistr_version = sistr.sistr_version
  String? sistr_predicted_serotype = sistr.sistr_predicted_serotype
  File? seqsero2_report = seqsero2.seqsero2_report
  String? seqsero2_version = seqsero2.seqsero2_version
  String? seqsero2_predicted_antigenic_profile = seqsero2.seqsero2_predicted_antigenic_profile
  String? seqsero2_predicted_serotype = seqsero2.seqsero2_predicted_serotype
  String? seqsero2_predicted_contamination = seqsero2.seqsero2_predicted_contamination
  # Salmonella serotype Typhi typing
  File? genotyphi_report_tsv = genotyphi_task.genotyphi_report_tsv 
  File? genotyphi_mykrobe_json = genotyphi_task.genotyphi_mykrobe_json
  String? genotyphi_version = genotyphi_task.genotyphi_version
  String? genotyphi_species = genotyphi_task.genotyphi_species
  Float? genotyphi_st_probes_percent_coverage = genotyphi_task.genotyphi_st_probes_percent_coverage
  String? genotyphi_final_genotype = genotyphi_task.genotyphi_final_genotype
  String? genotyphi_genotype_confidence = genotyphi_task.genotyphi_genotype_confidence
  # Klebsiella Typing
  File? kleborate_output_file = kleborate.kleborate_output_file
  String? kleborate_version = kleborate.kleborate_version
  String? kleborate_docker = kleborate.kleborate_docker
  String? kleborate_key_resistance_genes = kleborate.kleborate_key_resistance_genes
  String? kleborate_genomic_resistance_mutations = kleborate.kleborate_genomic_resistance_mutations
  String? kleborate_mlst_sequence_type = kleborate.kleborate_mlst_sequence_type
  String? kleborate_klocus = kleborate.kleborate_klocus
  String? kleborate_ktype = kleborate.kleborate_ktype
  String? kleborate_olocus = kleborate.kleborate_olocus
  String? kleborate_otype = kleborate.kleborate_otype
  String? kleborate_klocus_confidence = kleborate.kleborate_klocus_confidence
  String? kleborate_olocus_confidence = kleborate.kleborate_olocus_confidence
  # Acinetobacter Typing
  File? kaptive_output_file_k = kaptive.kaptive_output_file_k
  File? kaptive_output_file_oc = kaptive.kaptive_output_file_oc
  String? kaptive_version = kaptive.kaptive_version
  String? kaptive_k_match = kaptive.kaptive_k_match
  String? kaptive_k_type = kaptive.kaptive_k_type
  String? kaptive_k_confidence = kaptive.kaptive_k_confidence
  String? kaptive_oc_match = kaptive.kaptive_oc_match
  String? kaptive_oc_type = kaptive.kaptive_oc_type
  String? kaptive_oc_confidence = kaptive.kaptive_oc_confidence
  File? abricate_results = abricate.abricate_results
  String? abricate_genes = abricate.abricate_genes
  String? abricate_database = abricate.abricate_database
  String? abricate_version = abricate.abricate_version
  String? abricate_docker = abricate.abricate_docker
  # Mycobacterium Typing
  File? tbprofiler_output_file = tbprofiler.tbprofiler_output_csv
  File? tbprofiler_output_bam = tbprofiler.tbprofiler_output_bam
  File? tbprofiler_output_bai = tbprofiler.tbprofiler_output_bai
  String? tbprofiler_version = tbprofiler.version
  String? tbprofiler_main_lineage = tbprofiler.tbprofiler_main_lineage
  String? tbprofiler_sub_lineage = tbprofiler.tbprofiler_sub_lineage
  String? tbprofiler_dr_type = tbprofiler.tbprofiler_dr_type
  String? tbprofiler_resistance_genes = tbprofiler.tbprofiler_resistance_genes
  # Legionella pneumophila Typing
  File? legsta_results = legsta.legsta_results
  String? legsta_predicted_sbt = legsta.legsta_predicted_sbt
  String? legsta_version = legsta.legsta_version
  # Streptococcus pneumoniae Typing
  String? pbptyper_predicted_1A_2B_2X = pbptyper_task.pbptyper_predicted_1A_2B_2X
  File? pbptyper_pbptype_predicted_tsv = pbptyper_task.pbptyper_pbptype_predicted_tsv
  File? pbptyper_pbptype_1A_tsv = pbptyper_task.pbptyper_pbptype_1A_tsv
  File? pbptyper_pbptype_2B_tsv = pbptyper_task.pbptyper_pbptype_2B_tsv
  File? pbptyper_pbptype_2X_tsv = pbptyper_task.pbptyper_pbptype_2X_tsv
  String? pbptyper_version = pbptyper_task.pbptyper_version
  String? pbptyper_docker = pbptyper_task.pbptyper_docker
  String? poppunk_gps_cluster = poppunk_task.poppunk_gps_cluster
  File? poppunk_gps_external_cluster_csv = poppunk_task.poppunk_gps_external_cluster_csv
  String? poppunk_GPS_db_version = poppunk_task.poppunk_GPS_db_version
  String? poppunk_version = poppunk_task.poppunk_version
  String? poppunk_docker = poppunk_task.poppunk_docker
  String? seroba_version = seroba_task.seroba_version
  String? seroba_docker = seroba_task.seroba_docker
  String? seroba_serotype = seroba_task.seroba_serotype
  String? seroba_ariba_serotype = seroba_task.seroba_ariba_serotype
  String? seroba_ariba_identity = seroba_task.seroba_ariba_identity
  File? seroba_details = seroba_task.seroba_details
 }
}