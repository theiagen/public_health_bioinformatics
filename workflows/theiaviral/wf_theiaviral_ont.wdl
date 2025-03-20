version 1.0
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/contamination/task_metabuli.wdl" as metabuli_task
import "../../tasks/assembly/task_flye.wdl" as flye_task
import "../../tasks/phylogenetic_inference/utilities/task_skani.wdl" as skani_task

workflow theiaviral_ont{
  meta {
    description: "..."
  }
  input {
    File read1
    String taxon_of_interest
    String samplename
    File metabuli_db # delete this later only used for miniwdl testing
    File taxonomy_path # delete this later only used for miniwdl testing
    File skani_db # delete this later only used for miniwdl testing
  }
  call metabuli_task.metabuli as metabuli {
    input:
      read1 = read1,
      samplename = samplename,
      taxon_of_interest = taxon_of_interest,
      metabuli_db = metabuli_db,
      taxonomy_path = taxonomy_path
  }
  call flye_task.flye as flye {
    input:
      read1 = metabuli.metabuli_read1_extract,
  call skani_task.skani as skani {
    input:
      assembly_fasta = flye.assembly_fasta,
      samplename = samplename,
      skani_db = skani_db
  }
      samplename = samplename
  }
  call versioning.version_capture {
    input:
  }
  output {
    # metabuli outputs
    File metabuli_report = metabuli.metabuli_report
    File metabuli_classified = metabuli.metabuli_classified
    File metabuli_read1_extract = metabuli.metabuli_read1_extract
    String metabuli_database = metabuli.metabuli_database
    String metabuli_version = metabuli.metabuli_version
    String metabuli_docker = metabuli.metabuli_docker

    # flye outputs
    File flye_assembly = flye.assembly_fasta
    File flye_assembly_graph = flye.assembly_graph_gfa
    File flye_assembly_info = flye.assembly_info
    String flye_docker = flye.flye_docker
    # skani outputs - ANI-based reference genome selection
    File skani_report = skani.skani_report
    String skani_top_ani_accession = skani.skani_top_ani_accession
    String skani_database = skani.skani_database
    String skani_version = skani.skani_version
    String skani_docker = skani.skani_docker
    # versioning outputs
    String theiaviral_ont_version = version_capture.phb_version
    String theiaviral_ont_date = version_capture.date
  }
}