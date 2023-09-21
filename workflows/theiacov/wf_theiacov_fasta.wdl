
version 1.0

import "../../tasks/utilities/task_theiacov_fasta_utilities.wdl" as fasta_utilities
import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/gene_typing/task_abricate.wdl" as abricate
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/quality_control/task_qc_check_phb.wdl" as qc_check
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_fasta {
  meta {
    description: "Assessment of the quality of a consensus assembly fasta file for sars-cov-2, MPXV, WNV, flu, or RSV"
  }
  input {
    String samplename
    File assembly_fasta
    String organism = "sars-cov-2" # options: "sars-cov-2" "MPXV" "WNV" "flu" "rsv_a" "rsv_b
    File? reference_genome
    # reference genomes for various pathogens
    File ref_sars_cov_2 = "gs://theiagen-public-files-rp/terra/augur-sars-cov-2-references/MN908947.fasta"
    File ref_mpxv = "gs://theiagen-public-files-rp/terra/augur-mpox-references/reconstructed_ancestral_mpox.fasta"
    File ref_wnv = "gs://theiagen-public-files-rp/terra/augur-wnv-references/NC_063383.1.reference.fasta"
    String flu_segment = "HA" # options: HA or NA
    String? flu_subtype # options: "Victoria" "Yamagata" "H3N2" "H1N1"
    File ref_flu_h1n1_ha = "gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_ha.fasta"
    File ref_flu_h1n1_na = "gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_na.fasta"
    File ref_flu_h3n2_ha = "gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_ha.fasta"
    File ref_flu_h3n2_na = "gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_na.fasta"
    File ref_flu_yam_ha = "gs://theiagen-public-files-rp/terra/flu-references/reference_yam_ha.fasta"
    File ref_flu_yam_na = "gs://theiagen-public-files-rp/terra/flu-references/reference_yam_na.fasta"
    File ref_flu_vic_ha = "gs://theiagen-public-files-rp/terra/flu-references/reference_vic_ha.fasta"
    File ref_flu_vic_na = "gs://theiagen-public-files-rp/terra/flu-references/reference_vic_na.fasta"
    File ref_rsv_a = "gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_a.fasta"
    File ref_rsv_b = "gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_b.fasta"
    # nextclade inputs
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2023-08-17T12:00:00Z"
    String nextclade_dataset_name = "sars-cov-2"
    String nextclade_mpxv_tag = "2023-08-01T12:00:00Z"
    String nextclade_rsv_a_tag = "2023-02-03T12:00:00Z"
    String nextclade_rsv_b_tag = "2023-02-03T12:00:00Z"
    # nextclade flu old inputs
    #String nextclade_flu_h1n1_ha_tag = "2023-04-02T12:00:00Z"
    #String nextclade_flu_h1n1_na_tag = "2023-04-02T12:00:00Z"
    #String nextclade_flu_h3n2_ha_tag = "2023-04-02T12:00:00Z"
    #String nextclade_flu_h3n2_na_tag = "2023-04-02T12:00:00Z"
    #String nextclade_flu_vic_ha_tag = "2023-04-02T12:00:00Z"
    #String nextclade_flu_vic_na_tag = "2023-04-02T12:00:00Z"
    #String nextclade_flu_yam_tag = "2022-07-27T12:00:00Z"
    # nextclade flu inputs
    String nextclade_flu_h1n1_ha_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_h1n1_na_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_h3n2_ha_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_h3n2_na_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_vic_ha_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_vic_na_tag = "2023-08-10T12:00:00Z"
    String nextclade_flu_yam_tag = "2022-07-27T12:00:00Z"
    # sequencing values
    String seq_method
    String input_assembly_method
    # qc check parameters
    File? qc_check_table
    # vadr parameters
    Int? maxlen
  }
  call fasta_utilities.create_dummy_file {
    input:
  }
  # sars-cov-2 specific tasks
  if (organism == "sars-cov-2") {
    Int sars_cov_2_maxlen = 30000
    call consensus_qc_task.consensus_qc as sarscov2_consensus_qc {
      input:
        assembly_fasta = assembly_fasta,
        reference_genome = select_first([reference_genome,ref_sars_cov_2])
    }
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = assembly_fasta
    }
    call nextclade_task.nextclade as sarscov2_nextclade {
      input:
        genome_fasta = assembly_fasta,
        dataset_name = nextclade_dataset_name,
        dataset_reference = nextclade_dataset_reference,
        dataset_tag = nextclade_dataset_tag
    }
  }
  # MPXV specific tasks
  if (organism == "MPXV") {
    Int mpxv_maxlen = 200000
    call consensus_qc_task.consensus_qc as mpxv_consensus_qc {
      input:
        assembly_fasta = assembly_fasta,
        reference_genome = select_first([reference_genome,ref_mpxv])
    }
    call nextclade_task.nextclade as mpxv_nextclade {
      input:
        genome_fasta = assembly_fasta,
        dataset_name = "MPXV",
        dataset_reference = "ancestral",
        dataset_tag = nextclade_mpxv_tag
    }
  }
  # RSV specific tasks
  if (organism == "rsv_a") {
    Int rsv_a_maxlen = 16000
    call consensus_qc_task.consensus_qc as rsv_a_consensus_qc {
      input:
        assembly_fasta = assembly_fasta,
        reference_genome = select_first([reference_genome,ref_rsv_a])
    }
    call nextclade_task.nextclade as rsv_a_nextclade {
      input:
        genome_fasta = assembly_fasta,
        dataset_name = "rsv_a",
        dataset_reference = "EPI_ISL_412866",
        dataset_tag = nextclade_rsv_a_tag
    }
  }
  if (organism == "rsv_b") {
    Int rsv_b_maxlen = 16000
    call consensus_qc_task.consensus_qc as rsv_b_consensus_qc {
      input:
        assembly_fasta = assembly_fasta,
        reference_genome = select_first([reference_genome,ref_rsv_b])
    }
    call nextclade_task.nextclade as rsv_b_nextclade {
      input:
        genome_fasta = assembly_fasta,
        dataset_name = "rsv_b",
        dataset_reference = "EPI_ISL_1653999",
        dataset_tag = nextclade_rsv_b_tag
    }
  }
  # WNV specific tasks
  if (organism == "WNV") {
    Int wnv_maxlen = 11000
    call consensus_qc_task.consensus_qc as wnv_consensus_qc {
      input:
        assembly_fasta = assembly_fasta,
        reference_genome = select_first([reference_genome,ref_wnv])
    }
  }
  # flu specific tasks
  if (organism == "flu") {
    # If flu_subtype is not defined, run abricate to determine subtype
    if (!defined(flu_subtype)) {
      call abricate.abricate_flu {
              input:
                assembly = assembly_fasta,
                samplename = samplename,
                nextclade_flu_h1n1_ha_tag = nextclade_flu_h1n1_ha_tag,
                nextclade_flu_h1n1_na_tag = nextclade_flu_h1n1_na_tag,
                nextclade_flu_h3n2_ha_tag = nextclade_flu_h3n2_ha_tag,
                nextclade_flu_h3n2_na_tag = nextclade_flu_h3n2_na_tag,
                nextclade_flu_vic_ha_tag = nextclade_flu_vic_ha_tag,
                nextclade_flu_vic_na_tag = nextclade_flu_vic_na_tag,
                nextclade_flu_yam_tag = nextclade_flu_yam_tag
      }
      String abricate_subtype = abricate_flu.abricate_flu_subtype
    }
    if (flu_segment == "HA") {
        if ((defined(flu_subtype) && flu_subtype == "H1N1") || (defined(abricate_subtype) && (abricate_subtype == "H1N1" || abricate_subtype == "H1"))) {
          call consensus_qc_task.consensus_qc as flu_h1n1_ha_consensus_qc {
            input:
              assembly_fasta = assembly_fasta,
              reference_genome = select_first([reference_genome,ref_flu_h1n1_ha])
          }
          call nextclade_task.nextclade as flu_h1n1_ha_nextclade {
            input:
              genome_fasta = assembly_fasta,
              dataset_name = "flu_h1n1pdm_ha",
              dataset_reference = "MW626062",
              dataset_tag = nextclade_flu_h1n1_ha_tag
          }
        }
        if ((defined(flu_subtype) && flu_subtype == "H3N2") || (defined(abricate_subtype) && (abricate_subtype == "H3N2" || abricate_subtype == "H3"))) {
          call consensus_qc_task.consensus_qc as flu_h3n2_ha_consensus_qc {
            input:
              assembly_fasta = assembly_fasta,
              reference_genome = select_first([reference_genome,ref_flu_h3n2_ha])
          }
          call nextclade_task.nextclade as flu_h3n2_ha_nextclade {
            input:
              genome_fasta = assembly_fasta,
              dataset_name = "flu_h3n2_ha",
              dataset_reference = "EPI1857216",
              dataset_tag = nextclade_flu_h3n2_ha_tag
          }
        }
        if ((defined(flu_subtype) && flu_subtype == "Yamagata") || (defined(abricate_subtype) && abricate_subtype == "Yamagata")) {
          call consensus_qc_task.consensus_qc as flu_yam_ha_consensus_qc {
            input:
              assembly_fasta = assembly_fasta,
              reference_genome = select_first([reference_genome,ref_flu_yam_ha])
          }
          call nextclade_task.nextclade as flu_yam_ha_nextclade {
            input:
              genome_fasta = assembly_fasta,
              dataset_name = "flu_yam_ha",
              dataset_reference = "JN993010",
              dataset_tag = nextclade_flu_yam_tag
          }
        }
        if ((defined(flu_subtype) && flu_subtype == "Victoria") || (defined(abricate_subtype) && abricate_subtype == "Victoria")) {
          call consensus_qc_task.consensus_qc as flu_vic_ha_consensus_qc {
            input:
              assembly_fasta = assembly_fasta,
              reference_genome = select_first([reference_genome,ref_flu_vic_ha])
          }
          call nextclade_task.nextclade as flu_vic_ha_nextclade {
            input:
              genome_fasta = assembly_fasta,
              dataset_name = "flu_vic_ha",
              dataset_reference = "KX058884",
              dataset_tag = nextclade_flu_vic_ha_tag
          }
        }
    }
    if (flu_segment == "NA") {
      if ((defined(flu_subtype) && flu_subtype == "H1N1") || (defined(abricate_subtype) && (abricate_subtype == "H1N1" || abricate_subtype == "N1"))) {
        call consensus_qc_task.consensus_qc as flu_h1n1_na_consensus_qc {
          input:
            assembly_fasta = assembly_fasta,
            reference_genome = select_first([reference_genome,ref_flu_h1n1_na])
        }
        call nextclade_task.nextclade as flu_h1n1_na_nextclade {
          input:
            genome_fasta = assembly_fasta,
            dataset_name = "flu_h1n1pdm_na",
            dataset_reference = "MW626056",
            dataset_tag = nextclade_flu_h1n1_na_tag
        }
      }
      if ((defined(flu_subtype) && flu_subtype == "H3N2") || (defined(abricate_subtype) && (abricate_subtype == "H3N2" || abricate_subtype == "N2"))) {
        call consensus_qc_task.consensus_qc as flu_h3n2_na_consensus_qc {
          input:
            assembly_fasta = assembly_fasta,
            reference_genome = select_first([reference_genome,ref_flu_h3n2_na])
        }
        call nextclade_task.nextclade as flu_h3n2_na_nextclade {
          input:
            genome_fasta = assembly_fasta,
            dataset_name = "flu_h3n2_na",
            dataset_reference = "EPI1857215",
            dataset_tag = nextclade_flu_h3n2_na_tag
        }
      }
      if ((defined(flu_subtype) && flu_subtype == "Yamagata") || (defined(abricate_subtype) && abricate_subtype == "Yamagata")) {
        call consensus_qc_task.consensus_qc as flu_yam_na_consensus_qc {
          input:
            assembly_fasta = assembly_fasta,
            reference_genome = select_first([reference_genome,ref_flu_yam_na])
        }
      }
      if ((defined(flu_subtype) && flu_subtype == "Victoria") || (defined(abricate_subtype) && abricate_subtype == "Victoria")) {
        call consensus_qc_task.consensus_qc as flu_vic_na_consensus_qc {
          input:
            assembly_fasta = assembly_fasta,
            reference_genome = select_first([reference_genome,ref_flu_vic_na])
        }
        call nextclade_task.nextclade as flu_vic_na_nextclade {
          input:
            genome_fasta = assembly_fasta,
            dataset_name = "flu_vic_na",
            dataset_reference = "CY073894",
            dataset_tag = nextclade_flu_vic_na_tag
        }
      }
    }
  }
  # nextclade parser task
  if (organism == "sars-cov-2" || organism == "MPXV" || organism == "rsv_a" || organism == "rsv_b" || organism == "flu") {
    if (defined(sarscov2_nextclade.nextclade_tsv) || defined(mpxv_nextclade.nextclade_tsv) || defined(rsv_a_nextclade.nextclade_tsv) || defined(rsv_b_nextclade.nextclade_tsv) || defined(flu_h1n1_ha_nextclade.nextclade_tsv) || defined(flu_h1n1_na_nextclade.nextclade_tsv) || defined(flu_h3n2_ha_nextclade.nextclade_tsv) || defined(flu_h3n2_na_nextclade.nextclade_tsv) || defined(flu_yam_ha_nextclade.nextclade_tsv) || defined(flu_vic_ha_nextclade.nextclade_tsv) || defined(flu_vic_na_nextclade.nextclade_tsv)) {
      call nextclade_task.nextclade_output_parser {
        input:
          nextclade_tsv = select_first([sarscov2_nextclade.nextclade_tsv, mpxv_nextclade.nextclade_tsv, rsv_a_nextclade.nextclade_tsv, rsv_b_nextclade.nextclade_tsv, flu_h1n1_ha_nextclade.nextclade_tsv, flu_h1n1_na_nextclade.nextclade_tsv, flu_h3n2_ha_nextclade.nextclade_tsv, flu_h3n2_na_nextclade.nextclade_tsv, flu_yam_ha_nextclade.nextclade_tsv, flu_vic_ha_nextclade.nextclade_tsv, flu_vic_na_nextclade.nextclade_tsv]),
          organism = organism
      }
    }
  }
  # vadr task
  if (organism == "sars-cov-2" || organism == "MPXV" || organism == "rsv_a" || organism == "rsv_b" || organism == "WNV") {
    call vadr_task.vadr {
      input:
        genome_fasta = assembly_fasta,
        maxlen = select_first([maxlen, sars_cov_2_maxlen, mpxv_maxlen, rsv_a_maxlen, rsv_b_maxlen, wnv_maxlen]),
        assembly_length_unambiguous = select_first([sarscov2_consensus_qc.number_ATCG, mpxv_consensus_qc.number_ATCG, rsv_a_consensus_qc.number_ATCG, rsv_b_consensus_qc.number_ATCG, wnv_consensus_qc.number_ATCG])
    }
  }
  # QC check task
  if(defined(qc_check_table)) {
    call qc_check.qc_check_phb {
      input:
        qc_check_table = qc_check_table,
        expected_taxon = organism,
        number_N = select_first([sarscov2_consensus_qc.number_N, mpxv_consensus_qc.number_N, rsv_a_consensus_qc.number_N, rsv_b_consensus_qc.number_N, wnv_consensus_qc.number_N, flu_h1n1_ha_consensus_qc.number_N, flu_h1n1_na_consensus_qc.number_N, flu_h3n2_ha_consensus_qc.number_N, flu_h3n2_na_consensus_qc.number_N, flu_yam_ha_consensus_qc.number_N, flu_yam_na_consensus_qc.number_N, flu_vic_ha_consensus_qc.number_N, flu_vic_na_consensus_qc.number_N]),
        assembly_length_unambiguous = select_first([sarscov2_consensus_qc.number_ATCG, mpxv_consensus_qc.number_ATCG, rsv_a_consensus_qc.number_ATCG, rsv_b_consensus_qc.number_ATCG, wnv_consensus_qc.number_ATCG, flu_h1n1_ha_consensus_qc.number_ATCG, flu_h1n1_na_consensus_qc.number_ATCG, flu_h3n2_ha_consensus_qc.number_ATCG, flu_h3n2_na_consensus_qc.number_ATCG, flu_yam_ha_consensus_qc.number_ATCG, flu_yam_na_consensus_qc.number_ATCG, flu_vic_ha_consensus_qc.number_ATCG, flu_vic_na_consensus_qc.number_ATCG]),
        number_Degenerate = select_first([sarscov2_consensus_qc.number_Degenerate, mpxv_consensus_qc.number_Degenerate, rsv_a_consensus_qc.number_Degenerate, rsv_b_consensus_qc.number_Degenerate, wnv_consensus_qc.number_Degenerate, flu_h1n1_ha_consensus_qc.number_Degenerate, flu_h1n1_na_consensus_qc.number_Degenerate, flu_h3n2_ha_consensus_qc.number_Degenerate, flu_h3n2_na_consensus_qc.number_Degenerate, flu_yam_ha_consensus_qc.number_Degenerate, flu_yam_na_consensus_qc.number_Degenerate, flu_vic_ha_consensus_qc.number_Degenerate, flu_vic_na_consensus_qc.number_Degenerate]),
        percent_reference_coverage =  select_first([sarscov2_consensus_qc.percent_reference_coverage, mpxv_consensus_qc.percent_reference_coverage, rsv_a_consensus_qc.percent_reference_coverage, rsv_b_consensus_qc.percent_reference_coverage, wnv_consensus_qc.percent_reference_coverage, flu_h1n1_ha_consensus_qc.percent_reference_coverage, flu_h1n1_na_consensus_qc.percent_reference_coverage, flu_h3n2_ha_consensus_qc.percent_reference_coverage, flu_h3n2_na_consensus_qc.percent_reference_coverage, flu_yam_ha_consensus_qc.percent_reference_coverage, flu_yam_na_consensus_qc.percent_reference_coverage, flu_vic_ha_consensus_qc.percent_reference_coverage, flu_vic_na_consensus_qc.percent_reference_coverage]),
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
    String number_N = select_first([sarscov2_consensus_qc.number_N, mpxv_consensus_qc.number_N, rsv_a_consensus_qc.number_N, rsv_b_consensus_qc.number_N, wnv_consensus_qc.number_N, flu_h1n1_ha_consensus_qc.number_N, flu_h1n1_na_consensus_qc.number_N, flu_h3n2_ha_consensus_qc.number_N, flu_h3n2_na_consensus_qc.number_N, flu_yam_ha_consensus_qc.number_N, flu_yam_na_consensus_qc.number_N, flu_vic_ha_consensus_qc.number_N, flu_vic_na_consensus_qc.number_N, "untyped flu; not calculated"])
    String assembly_length_unambiguous = select_first([sarscov2_consensus_qc.number_ATCG, mpxv_consensus_qc.number_ATCG, rsv_a_consensus_qc.number_ATCG, rsv_b_consensus_qc.number_ATCG, wnv_consensus_qc.number_ATCG, flu_h1n1_ha_consensus_qc.number_ATCG, flu_h1n1_na_consensus_qc.number_ATCG, flu_h3n2_ha_consensus_qc.number_ATCG, flu_h3n2_na_consensus_qc.number_ATCG, flu_yam_ha_consensus_qc.number_ATCG, flu_yam_na_consensus_qc.number_ATCG, flu_vic_ha_consensus_qc.number_ATCG, flu_vic_na_consensus_qc.number_ATCG, "untyped flu; not calculated"])
    String number_Degenerate = select_first([sarscov2_consensus_qc.number_Degenerate, mpxv_consensus_qc.number_Degenerate, rsv_a_consensus_qc.number_Degenerate, rsv_b_consensus_qc.number_Degenerate, wnv_consensus_qc.number_Degenerate, flu_h1n1_ha_consensus_qc.number_Degenerate, flu_h1n1_na_consensus_qc.number_Degenerate, flu_h3n2_ha_consensus_qc.number_Degenerate, flu_h3n2_na_consensus_qc.number_Degenerate, flu_yam_ha_consensus_qc.number_Degenerate, flu_yam_na_consensus_qc.number_Degenerate, flu_vic_ha_consensus_qc.number_Degenerate, flu_vic_na_consensus_qc.number_Degenerate, "untyped flu; not calculated"])
    String number_Total = select_first([sarscov2_consensus_qc.number_Total, mpxv_consensus_qc.number_Total, rsv_a_consensus_qc.number_Total, rsv_b_consensus_qc.number_Total, wnv_consensus_qc.number_Total, flu_h1n1_ha_consensus_qc.number_Total, flu_h1n1_na_consensus_qc.number_Total, flu_h3n2_ha_consensus_qc.number_Total, flu_h3n2_na_consensus_qc.number_Total, flu_yam_ha_consensus_qc.number_Total, flu_yam_na_consensus_qc.number_Total, flu_vic_ha_consensus_qc.number_Total, flu_vic_na_consensus_qc.number_Total, "untyped flu; not calculated"])
    String percent_reference_coverage = select_first([sarscov2_consensus_qc.percent_reference_coverage, mpxv_consensus_qc.percent_reference_coverage, rsv_a_consensus_qc.percent_reference_coverage, rsv_b_consensus_qc.percent_reference_coverage, wnv_consensus_qc.percent_reference_coverage, flu_h1n1_ha_consensus_qc.percent_reference_coverage, flu_h1n1_na_consensus_qc.percent_reference_coverage, flu_h3n2_ha_consensus_qc.percent_reference_coverage, flu_h3n2_na_consensus_qc.percent_reference_coverage, flu_yam_ha_consensus_qc.percent_reference_coverage, flu_yam_na_consensus_qc.percent_reference_coverage, flu_vic_ha_consensus_qc.percent_reference_coverage, flu_vic_na_consensus_qc.percent_reference_coverage, "untyped flu; not calculated"])
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
    File? nextclade_json = select_first([sarscov2_nextclade.nextclade_json, mpxv_nextclade.nextclade_json, rsv_a_nextclade.nextclade_json, rsv_b_nextclade.nextclade_json, flu_h1n1_ha_nextclade.nextclade_json, flu_h1n1_na_nextclade.nextclade_json, flu_h3n2_ha_nextclade.nextclade_json, flu_h3n2_na_nextclade.nextclade_json, flu_yam_ha_nextclade.nextclade_json, flu_vic_ha_nextclade.nextclade_json, flu_vic_na_nextclade.nextclade_json, create_dummy_file.outputFile])
    File? auspice_json = select_first([sarscov2_nextclade.auspice_json, mpxv_nextclade.auspice_json, rsv_a_nextclade.auspice_json, rsv_b_nextclade.auspice_json, flu_h1n1_ha_nextclade.auspice_json, flu_h1n1_na_nextclade.auspice_json, flu_h3n2_ha_nextclade.auspice_json, flu_h3n2_na_nextclade.auspice_json, flu_yam_ha_nextclade.auspice_json, flu_vic_ha_nextclade.auspice_json, flu_vic_na_nextclade.auspice_json, create_dummy_file.outputFile])
    File nextclade_tsv = select_first([sarscov2_nextclade.nextclade_tsv, mpxv_nextclade.nextclade_tsv, rsv_a_nextclade.nextclade_tsv, rsv_b_nextclade.nextclade_tsv, flu_h1n1_ha_nextclade.nextclade_tsv, flu_h1n1_na_nextclade.nextclade_tsv, flu_h3n2_ha_nextclade.nextclade_tsv, flu_h3n2_na_nextclade.nextclade_tsv, flu_yam_ha_nextclade.nextclade_tsv, flu_vic_ha_nextclade.nextclade_tsv, flu_vic_na_nextclade.nextclade_tsv, create_dummy_file.outputFile])
    String? nextclade_version = select_first([sarscov2_nextclade.nextclade_version, mpxv_nextclade.nextclade_version, rsv_a_nextclade.nextclade_version, rsv_b_nextclade.nextclade_version, flu_h1n1_ha_nextclade.nextclade_version, flu_h1n1_na_nextclade.nextclade_version, flu_h3n2_ha_nextclade.nextclade_version, flu_h3n2_na_nextclade.nextclade_version, flu_yam_ha_nextclade.nextclade_version, flu_vic_ha_nextclade.nextclade_version, flu_vic_na_nextclade.nextclade_version, "not applicable"])
    String? nextclade_docker = select_first([sarscov2_nextclade.nextclade_docker, mpxv_nextclade.nextclade_docker, rsv_a_nextclade.nextclade_docker, rsv_b_nextclade.nextclade_docker, flu_h1n1_ha_nextclade.nextclade_docker, flu_h1n1_na_nextclade.nextclade_docker, flu_h3n2_ha_nextclade.nextclade_docker, flu_h3n2_na_nextclade.nextclade_docker, flu_yam_ha_nextclade.nextclade_docker, flu_vic_ha_nextclade.nextclade_docker, flu_vic_na_nextclade.nextclade_docker, "not applicable"])
    String nextclade_ds_tag = select_first([sarscov2_nextclade.nextclade_dataset_tag, mpxv_nextclade.nextclade_dataset_tag, rsv_a_nextclade.nextclade_dataset_tag, rsv_b_nextclade.nextclade_dataset_tag, flu_h1n1_ha_nextclade.nextclade_dataset_tag, flu_h1n1_na_nextclade.nextclade_dataset_tag, flu_h3n2_ha_nextclade.nextclade_dataset_tag, flu_h3n2_na_nextclade.nextclade_dataset_tag, flu_yam_ha_nextclade.nextclade_dataset_tag, flu_vic_ha_nextclade.nextclade_dataset_tag, flu_vic_na_nextclade.nextclade_dataset_tag, "not applicable"])
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