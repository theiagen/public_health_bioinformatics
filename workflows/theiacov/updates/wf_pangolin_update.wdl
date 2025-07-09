version 1.0

import "../../../tasks/species_typing/betacoronavirus/task_pangolin.wdl" as pangolin
import "../../../tasks/task_versioning.wdl" as versioning
import "../../../workflows/utilities/wf_organism_parameters.wdl" as set_organism_defaults

workflow pangolin_update {
  input {
    String samplename
    File assembly_fasta
    String old_lineage
    String old_pangolin_docker
    String old_pangolin_assignment_version
    String old_pangolin_versions
    String? new_pangolin_docker
    String organism = "sars-cov-2"
    File? lineage_log
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = organism,
      pangolin_docker_image = new_pangolin_docker,
      # including these to block from terra
      flu_segment = "",
      flu_subtype = "",
      reference_gff_file = "gs://theiagen-public-resources-rp/empty_files/empty.gff3",
      reference_genome = "gs://theiagen-public-resources-rp/empty_files/empty.fasta",
      genome_length_input = 0,
      nextclade_dataset_tag_input = "",
      nextclade_dataset_name_input = "",     
      vadr_max_length = 0,
      vadr_skip_length = 0,
      vadr_mem = 0, 
      vadr_options = "",
      primer_bed_file = "gs://theiagen-public-resources-rp/empty_files/empty.bed",
      gene_locations_bed_file = "gs://theiagen-public-resources-rp/empty_files/empty.bed",
      kraken_target_organism_input = "",
      hiv_primer_version = ""
  }
  call pangolin.pangolin4 {
    input:
      samplename = samplename,
      fasta = assembly_fasta,
      docker = organism_parameters.pangolin_docker
  }
  call pangolin.pangolin_update_log {
    input:
      samplename = samplename,
      old_lineage = old_lineage,
      old_pangolin_docker = old_pangolin_docker,
      old_pangolin_assignment_version = old_pangolin_assignment_version,
      old_pangolin_versions = old_pangolin_versions,
      new_lineage = pangolin4.pangolin_lineage,
      new_pangolin_docker = pangolin4.pangolin_docker,
      new_pangolin_assignment_version = pangolin4.pangolin_assignment_version,
      new_pangolin_versions = pangolin4.pangolin_versions,
      lineage_log = lineage_log
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String pangolin_update_version = version_capture.phb_version
    String pangolin_update_analysis_date = version_capture.date
    # Pangolin Assignments
    String pango_lineage = pangolin4.pangolin_lineage
    String pangolin_conflicts = pangolin4.pangolin_conflicts
    String pangolin_notes = pangolin4.pangolin_notes
    String pangolin_assignment_version = pangolin4.pangolin_assignment_version
    String pangolin_versions = pangolin4.pangolin_versions
    File pango_lineage_report = pangolin4.pango_lineage_report
    String pangolin_docker = pangolin4.pangolin_docker
    String pango_lineage_expanded = pangolin4.pangolin_lineage_expanded
    # Update Log
    String pangolin_updates = pangolin_update_log.pangolin_updates
    File pango_lineage_log = pangolin_update_log.pango_lineage_log
  }
}