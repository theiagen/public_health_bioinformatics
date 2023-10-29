version 1.0

import "../../tasks/utilities/task_file_handling.wdl" as concatenate
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin_task
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/utilities/task_theiacov_fasta_batch.wdl" as theiacov_fasta_wrangling_task
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_fasta_batch {
  meta {
    description: "TheiaCoV_FASTA for multiple samples"
  }
  input {
    Array[String] samplenames
    Array[File] assembly_fastas
    String organism = "sars-cov-2"
    # sequencing values
    String seq_method
    String input_assembly_method
    # nextclade inputs
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2023-09-21T12:00:00Z"
    String? nextclade_dataset_name
  }
  call concatenate.cat_files {
    input: 
      files_to_cat = assembly_fastas,
      concatenated_file_name = "concatenated_assemblies.fasta"
  }
  if (organism == "sars-cov-2") {
    # sars-cov-2 specific tasks
    call pangolin_task.pangolin4 {
      input:
        samplename = "placeholder", # what should go here?
        fasta = cat_files.concatenated_files
    }
  }
  if (organism == "MPXV" || organism == "sars-cov-2"){
    # tasks specific to either MPXV or sars-cov-2 
    call nextclade_task.nextclade {
      input:
      genome_fasta = cat_files.concatenated_files,
      dataset_name = select_first([nextclade_dataset_name, organism]),
      dataset_reference = nextclade_dataset_reference,
      dataset_tag = nextclade_dataset_tag
    }
  }
  call theiacov_fasta_wrangling_task.theiacov_fasta_batch {
    input:
      samplenames = samplenames,
      organism = organism,
      pangolin_lineage_report = pangolin4.pango_lineage_report,
      nextclade_tsv = nextclade.nextclade_tsv,
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
    # Pangolin outputs
    File? pango_lineage_report = pangolin4.pango_lineage_report
    # Nextclade outputs
    File? nextclade_json = nextclade.nextclade_json
    File? auspice_json = nextclade.auspice_json
    File? nextclade_tsv = nextclade.nextclade_tsv
  }
}