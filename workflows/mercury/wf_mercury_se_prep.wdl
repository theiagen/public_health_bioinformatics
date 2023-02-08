version 1.0

import "../tasks/task_versioning.wdl" as versioning
import "../tasks/task_pub_repo_prep.wdl" as submission_prep

workflow mercury_se_prep {
  input {
    # Required Files
    File assembly_fasta
    File reads_dehosted
    # Required Metadata (TheiaCoV GC Outputs)
    String assembly_method
    Float assembly_mean_coverage
    # Required Metadata (User Inputs)
    String authors
    String bioproject_accession
    String collecting_lab
    String collecting_lab_address
    String collection_date
    String continent
    String country
    String dehosting_method = "NCBI Human Scrubber"
    String gisaid_submitter
    String gisaid_organism = "hCoV-19"
    String filetype = "fastq"
    String host ="Human"
    String host_disease
    String host_sci_name = "Homo sapiens"
    String instrument_model
    String isolation_source
    String library_id
    String library_layout = "single"
    String library_selection
    String library_source
    String library_strategy
    String organism
    Int number_N
    String seq_platform
    String state
    String submission_id
    String submitting_lab
    String submitting_lab_address
    # Optional Metadata
    String? amplicon_primer_scheme
    String? amplicon_size
    String? biosample_accession
    String? gisaid_accession
    String? county
    String? patient_gender
    String? patient_age
    String? purpose_of_sampling
    String? purpose_of_sequencing
    String? submitter_email
    String? treatment
    # Optional User-Defined Thresholds for Generating Submission Files
    Int number_N_threshold = 5000
  }
  if (number_N <= number_N_threshold) {
    call submission_prep.ncbi_prep_one_sample_se {
      input:
        amplicon_primer_scheme = amplicon_primer_scheme,
        amplicon_size = amplicon_size,
        assembly_fasta = assembly_fasta,
        assembly_method = assembly_method,
        bioproject_accession = bioproject_accession,
        biosample_accession = biosample_accession,
        collecting_lab = collecting_lab,
        collection_date = collection_date,
        country = country,
        dehosting_method = dehosting_method,
        design_description = "Whole genome sequencing of ~{organism}",
        filetype = filetype,
        gisaid_accession = gisaid_accession,
        gisaid_organism = gisaid_organism,
        host = host,
        host_disease = host_disease,
        host_sci_name = host_sci_name,
        instrument_model = instrument_model,
        isolation_source = isolation_source,
        library_id = library_id,
        library_layout = library_layout,
        library_selection = library_selection,
        library_source = library_source,
        library_strategy = library_strategy,
        organism = organism,
        patient_age = patient_age,
        patient_gender = patient_gender,
        purpose_of_sampling = purpose_of_sampling,
        purpose_of_sequencing = purpose_of_sequencing,
        reads_dehosted = reads_dehosted,
        seq_platform = seq_platform,
        state = state,
        submission_id = submission_id,
        submitter_email = submitter_email,
        treatment = treatment
    }
    call submission_prep.gisaid_prep_one_sample {
      input:
        assembly_fasta = assembly_fasta,
        authors = authors,
        assembly_method = assembly_method,
        collecting_lab = collecting_lab,
        collecting_lab_address = collecting_lab_address,
        collection_date = collection_date,
        continent = continent,
        assembly_mean_coverage = assembly_mean_coverage,
        country = country,
        gisaid_submitter = gisaid_submitter,
        host = host,
        organism = gisaid_organism,
        seq_platform = seq_platform,
        state = state,
        submission_id = submission_id,
        submitting_lab = submitting_lab,
        submitting_lab_address = submitting_lab_address,
        county = county,
        patient_gender = patient_gender,
        patient_age = patient_age,
        purpose_of_sequencing = purpose_of_sequencing,
        treatment = treatment
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String mercury_pe_prep_version = version_capture.phvg_version
    String mercury_pe_prep_analysis_date = version_capture.date
    # NCBI Submission Files
    File? biosample_attributes = ncbi_prep_one_sample_se.biosample_attributes
    File? sra_metadata = ncbi_prep_one_sample_se.sra_metadata
    File? genbank_assembly = ncbi_prep_one_sample_se.genbank_assembly
    File? genbank_modifier = ncbi_prep_one_sample_se.genbank_modifier
    File? sra_reads = ncbi_prep_one_sample_se.sra_reads
    # GISAID Submission Files
    File? gisaid_assembly = gisaid_prep_one_sample.gisaid_assembly
    File? gisaid_metadata = gisaid_prep_one_sample.gisaid_metadata
  }
}