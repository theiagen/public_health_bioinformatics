version 1.0
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/contamination/task_metabuli.wdl" as metabuli_task
import "../../tasks/assembly/task_flye.wdl" as flye_task

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
    String metabuli_docker = metabuli.metabuli_docker

    # flye outputs
    File flye_assembly = flye.assembly_fasta
    File flye_assembly_graph = flye.assembly_graph_gfa
    File flye_assembly_info = flye.assembly_info
    String flye_docker = flye.flye_docker

    # versioning outputs
    String theiaviral_ont_version = version_capture.phb_version
    String theiaviral_ont_date = version_capture.date
  }
}