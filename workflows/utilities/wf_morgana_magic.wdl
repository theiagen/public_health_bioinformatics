version 1.0

import "../../tasks/quality_control/advanced_metrics/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/species_typing/betacoronavirus/task_pangolin.wdl" as pangolin
import "../../tasks/species_typing/lentivirus/task_quasitools.wdl" as quasitools
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults

workflow morgana_magic {
  input {
    String samplename
    File assembly_fasta
    File read1
    File read2
    String taxon_id
  }
  #### need to add more flu characterization
  call set_organism_defaults.organism_parameters {
    input:
      taxon_id = taxon_id,
      organism = "To Be Determined"
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = assembly_fasta,
      reference_genome = organism_parameters.reference,
      genome_length = organism_parameters.genome_length
  }
  if (organism_parameters.standardized_organism == "sars-cov-2") {
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = assembly_fasta,
        docker = organism_parameters.pangolin_docker
    }
  }
  if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b") { 
    # tasks specific to either MPXV, sars-cov-2, or RSV-A/RSV-B
    call nextclade_task.nextclade_v3 {
      input:
        genome_fasta = assembly_fasta,
        dataset_name = organism_parameters.nextclade_dataset_name,
        dataset_tag = organism_parameters.nextclade_dataset_tag
    }
    call nextclade_task.nextclade_output_parser {
      input:
        nextclade_tsv = nextclade_v3.nextclade_tsv,
        organism = organism_parameters.standardized_organism
    }
  }
  ### add flu
  ##### is running vadr even something we want to do????
  if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "WNV" || organism_parameters.standardized_organism == "flu" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b") { 
    # tasks specific to MPXV, sars-cov-2, WNV, flu, rsv_a, and rsv_b
    call vadr_task.vadr {
      input:
        genome_fasta = assembly_fasta,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        vadr_opts = organism_parameters.vadr_opts,
        max_length = organism_parameters.vadr_maxlength,
        skip_length = organism_parameters.vadr_skiplength,
        memory = organism_parameters.vadr_memory
    }
  }
  ##### is running quasitools even something we want to do????
  if (organism_parameters.standardized_organism == "HIV") {
    call quasitools.quasitools as quasitools_illumina_pe {
      input:
        read1 = read1,
        read2 = read2,
        samplename = samplename
    }
  }
  output {
    String organism = organism_parameters.standardized_organism
    # Pangolin outputs
    String? pango_lineage = pangolin4.pangolin_lineage
    String? pango_lineage_expanded = pangolin4.pangolin_lineage_expanded
    String? pangolin_conflicts = pangolin4.pangolin_conflicts
    String? pangolin_notes = pangolin4.pangolin_notes
    String? pangolin_assignment_version = pangolin4.pangolin_assignment_version
    File? pango_lineage_report = pangolin4.pango_lineage_report
    String? pangolin_docker = pangolin4.pangolin_docker
    String? pangolin_versions = pangolin4.pangolin_versions
    # Nextclade outputs for all organisms
    String nextclade_version = select_first([nextclade_v3.nextclade_version, ""])
    String nextclade_docker = select_first([nextclade_v3.nextclade_docker, ""])
    # Nextclade outputs for non-flu
    File? nextclade_json = nextclade_v3.nextclade_json
    File? auspice_json = nextclade_v3.auspice_json
    File? nextclade_tsv = nextclade_v3.nextclade_tsv
    String nextclade_ds_tag = organism_parameters.nextclade_dataset_tag
    String? nextclade_aa_subs = nextclade_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser.nextclade_aa_dels
    String? nextclade_clade = nextclade_output_parser.nextclade_clade
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    String? nextclade_qc = nextclade_output_parser.nextclade_qc
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    File? vadr_feature_tbl_pass = vadr.feature_tbl_pass
    File? vadr_feature_tbl_fail = vadr.feature_tbl_fail
    File? vadr_classification_summary_file = vadr.classification_summary_file
    File? vadr_all_outputs_tar_gz = vadr.outputs_tgz
    String? vadr_num_alerts = vadr.num_alerts
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    # HIV Outputs
    String? quasitools_version = quasitools_illumina_pe.quasitools_version
    String? quasitools_date = quasitools_illumina_pe.quasitools_date
    File? quasitools_coverage_file = quasitools_illumina_pe.coverage_file
    File? quasitools_dr_report = quasitools_illumina_pe.dr_report
    File? quasitools_hydra_vcf = quasitools_illumina_pe.hydra_vcf
    File? quasitools_mutations_report = quasitools_illumina_pe.mutations_report
  }
}