
version 1.0

import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/gene_typing/task_abricate.wdl" as abricate
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/quality_control/task_qc_check_phb.wdl" as qc_check
import "../../tasks/utilities/task_organism_defaults.wdl" as defaults
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_fasta {
  meta {
    description: "Assessment of the quality of a consensus assembly fasta file for sars-cov-2, MPXV, WNV, flu, or RSV"
  }
  input {
    String samplename
    File assembly_fasta
    String organism = "sars-cov-2" # options: "sars-cov-2" "MPXV" "WNV" "flu" "rsv_a" "rsv_b
    # flu options
    String flu_segment = "HA" # options: HA or NA
    String? flu_subtype # options: "Victoria" "Yamagata" "H3N2" "H1N1"
    # optional reference information
    File? reference_genome
    Int? genome_length
    # nextclade inputs (default SC2)
    String? nextclade_dataset_reference
    String? nextclade_dataset_tag
    String? nextclade_dataset_name
    # sequencing values
    String seq_method
    String input_assembly_method
    # qc check parameters
    File? qc_check_table
    # vadr parameters
    Int? maxlen
    String? vadr_opts
  }
  if (organism == "sars-cov-2") {
    String sc2_reference_genome = "gs://theiagen-public-files-rp/terra/augur-sars-cov-2-references/MN908947.fasta"
    String sc2_nextclade_ds_tag = "2023-08-17T12:00:00Z"
    String sc2_nextclade_ref = "MN908947"
    String sc2_nextclade_ds_name = "sars-cov-2"
    Int sc2_genome_len = 29903
    Int sc2_vadr_max_length = 30000
    String sc2_vadr_options = "--noseqnamemax --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --out_allfasta"
  }
  if (organism == "MPXV") {
    String mpox_reference_genome = "gs://theiagen-public-files/terra/mpxv-files/MPXV.MT903345.reference.fasta"
    String mpox_nextclade_ds_tag = "2023-08-01T12:00:00Z"
    String mpox_nextclade_ref = "pseudo_ON563414"
    String mpox_nextclade_ds_name = "hMPXV_B1"
    # String mpox_target_org = "Monkeypox virus"
    # String mpox_primer_bed_file = "gs://theiagen-public-files/terra/mpxv-files/MPXV.primer.bed"
    # String mpox_reference_gff_file = "gs://theiagen-public-files/terra/mpxv-files/Mpox-MT903345.1.reference.gff3"
    String mpox_vadr_options = "--glsearch -s -r --nomisc --mkey mpxv --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --out_allfasta --minimap2 --s_overhang 150"
    Int mpox_vadr_max_length = 210000
    Int mpox_genome_len = 197200
  }
  if (organism == "WNV") {
    String wnv_reference_genome = "gs://theiagen-public-files/terra/theiacov-files/WNV/NC_009942.1_wnv_L1.fasta"
    # String wnv_target_org = "West Nile virus"
    # String wnv_primer_bed_file = "gs://theiagen-public-files/terra/theiacov-files/WNV/WNV-L1_primer.bed"
    Int wnv_genome_len = 11000
    String wnv_vadr_options = "--mkey flavi --mdir /opt/vadr/vadr-models-flavi/ --nomisc --noprotid --out_allfasta"    
    Int wnv_vadr_max_length = 11000
  }
  if (organism == "flu") {
    if (!defined(flu_subtype)) {
      call abricate.abricate_flu {
        input:
          assembly = assembly_fasta,
          samplename = samplename
      }
      String abricate_subtype = abricate_flu.abricate_flu_subtype
    }
    call defaults.set_organism_defaults_flu {
      input:
        flu_segment = flu_segment,
        flu_subtype = select_first([flu_subtype, abricate_subtype]),
        genome_len = genome_length
    }
  }
  if (organism == "rsv_a") {
    String rsv_a_reference_genome = "gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_a.fasta"
    String rsv_a_nextclade_ds_tag = "2023-02-03T12:00:00Z"
    String rsv_a_nextclade_ref = "EPI_ISL_412866"
    String rsv_a_nextclade_ds_name = "rsv_a"
    Int rsv_a_genome_len = 16000
    String rsv_a_vadr_options = "-r --mkey rsv --xnocomp"
    Int rsv_a_vadr_max_length = 15500
  } 
  if (organism == "rsv_b") {
    String rsv_b_reference_genome = "gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_b.fasta"
    String rsv_b_nextclade_ds_tag = "2023-02-03T12:00:00Z"
    String rsv_b_nextclade_ref = "EPI_ISL_1653999"
    String rsv_b_nextclade_ds_name = "rsv_b"
    Int rsv_b_genome_len = 16000   
    String rsv_b_vadr_options = "-r --mkey rsv --xnocomp"
    Int rsv_b_vadr_max_length = 15500
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = assembly_fasta,
      reference_genome = select_first([reference_genome, sc2_reference_genome, mpox_reference_genome, wnv_reference_genome, set_organism_defaults_flu.reference, rsv_a_reference_genome, rsv_b_reference_genome]),
      genome_length = select_first([genome_length, sc2_genome_len, mpox_genome_len, wnv_genome_len, set_organism_defaults_flu.genome_length, rsv_a_genome_len, rsv_b_genome_len]),
  }
  if (organism == "sars-cov-2") {
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = assembly_fasta
    }
  }
  if (organism == "sars-cov-2" || organism == "MPXV" || organism == "rsv_a" || organism == "rsv_b" || organism == "flu") {
    if (select_first([set_organism_defaults_flu.nextclade_dataset_tag, ""]) != "NA") {
      call nextclade_task.nextclade {
        input:
          genome_fasta = assembly_fasta,
          dataset_name = select_first([nextclade_dataset_name, sc2_nextclade_ds_name, mpox_nextclade_ds_name, rsv_a_nextclade_ds_name, rsv_b_nextclade_ds_name, set_organism_defaults_flu.nextclade_dataset_name]),
          dataset_reference = select_first([nextclade_dataset_reference, sc2_nextclade_ref, mpox_nextclade_ref, rsv_a_nextclade_ref, rsv_b_nextclade_ref, set_organism_defaults_flu.nextclade_reference]),
          dataset_tag = select_first([nextclade_dataset_tag, sc2_nextclade_ds_tag, mpox_nextclade_ds_tag, rsv_a_nextclade_ds_tag, rsv_b_nextclade_ds_tag, set_organism_defaults_flu.nextclade_dataset_tag]),
      }
    }
  }
  # nextclade parser task
  if (organism == "sars-cov-2" || organism == "MPXV" || organism == "rsv_a" || organism == "rsv_b" || organism == "flu") {
    if (defined(nextclade.nextclade_tsv)) {
      call nextclade_task.nextclade_output_parser {
        input:
          nextclade_tsv = select_first([nextclade.nextclade_tsv]),
          organism = organism
      }
    }
  }
  # vadr task
  if (organism == "sars-cov-2" || organism == "MPXV" || organism == "rsv_a" || organism == "rsv_b" || organism == "WNV") {
    call vadr_task.vadr {
      input:
        genome_fasta = assembly_fasta,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        maxlen = select_first([maxlen, sc2_vadr_max_length, mpox_vadr_max_length, wnv_vadr_max_length, rsv_a_vadr_max_length, rsv_b_vadr_max_length]),
        vadr_opts = select_first([vadr_opts, sc2_vadr_options, mpox_vadr_options, wnv_vadr_options, rsv_a_vadr_options, rsv_b_vadr_options]),
    }
  }
  # QC check task
  if(defined(qc_check_table)) {
    call qc_check.qc_check_phb {
      input:
        qc_check_table = qc_check_table,
        expected_taxon = organism,
        number_N = consensus_qc.number_N,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        number_Degenerate = consensus_qc.number_Degenerate,
        percent_reference_coverage =  consensus_qc.percent_reference_coverage,
        vadr_num_alerts = vadr.num_alerts
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String theiacov_fasta_version = version_capture.phb_version
    String theiacov_fasta_analysis_date = version_capture.date
    # Read & Assembly Metadata
    String seq_platform = seq_method
    String assembly_method = input_assembly_method
    # Assembly QC - consensus assembly summary statistics
    Int number_N = consensus_qc.number_N
    Int assembly_length_unambiguous = consensus_qc.number_ATCG
    Int number_Degenerate = consensus_qc.number_Degenerate
    Int number_Total = consensus_qc.number_Total
    Float percent_reference_coverage = consensus_qc.percent_reference_coverage
    # Pangolin outputs
    String? pango_lineage = pangolin4.pangolin_lineage
    String? pango_lineage_expanded = pangolin4.pangolin_lineage_expanded
    String? pangolin_conflicts = pangolin4.pangolin_conflicts
    String? pangolin_notes = pangolin4.pangolin_notes
    String? pangolin_assignment_version = pangolin4.pangolin_assignment_version
    File? pango_lineage_report = pangolin4.pango_lineage_report
    String? pangolin_docker = pangolin4.pangolin_docker
    String? pangolin_versions = pangolin4.pangolin_versions
    # Nextclade outputs
    File? nextclade_json = nextclade.nextclade_json
    File? auspice_json = nextclade.auspice_json
    File? nextclade_tsv = nextclade.nextclade_tsv
    String? nextclade_version = nextclade.nextclade_version
    String? nextclade_docker = nextclade.nextclade_docker
    String nextclade_ds_tag =  select_first([nextclade_dataset_tag, sc2_nextclade_ds_tag, mpox_nextclade_ds_tag, rsv_a_nextclade_ds_tag, rsv_b_nextclade_ds_tag, set_organism_defaults_flu.nextclade_dataset_tag, "NA"])
    String? nextclade_clade = nextclade_output_parser.nextclade_clade
    String? nextclade_aa_subs = nextclade_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser.nextclade_aa_dels
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    # VADR Annotation QC
    File?  vadr_alerts_list = vadr.alerts_list
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    # QC_Check Results
    String? qc_check = qc_check_phb.qc_check
    File? qc_standard = qc_check_phb.qc_standard
    # Flu Outputs
    String? abricate_flu_type = abricate_flu.abricate_flu_type
    String? abricate_flu_subtype =  abricate_flu.abricate_flu_subtype
    File? abricate_flu_results = abricate_flu.abricate_flu_results
    String? abricate_flu_database =  abricate_flu.abricate_flu_database
    String? abricate_flu_version = abricate_flu.abricate_flu_version
  }
}