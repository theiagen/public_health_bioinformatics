version 1.0

import "../tasks/gene_typing/task_amrfinderplus.wdl" as amrfindertask
import "../tasks/task_versioning.wdl" as versioning

workflow amrfinderplus_wf {
  input {
      File assembly
      String samplename
    }
  call amrfindertask.amrfinderplus_nuc {
    input:
      assembly = assembly,
      samplename = samplename
    }
  call versioning.version_capture{
    input:
  }
  output {
    String amrfinderplus_version = amrfinderplus_nuc.amrfinderplus_version
    String amrfinderplus_db_version = amrfinderplus_nuc.amrfinderplus_db_version
    String amrfinderplus_wf_version = version_capture.phbg_version
    String amrfinderplus_wf_analysis_date = version_capture.date
    File amrfinderplus_all_report = amrfinderplus_nuc.amrfinderplus_all_report
    File amrfinderplus_amr_report = amrfinderplus_nuc.amrfinderplus_amr_report
    File amrfinderplus_stress_report = amrfinderplus_nuc.amrfinderplus_stress_report
    File amrfinderplus_virulence_report = amrfinderplus_nuc.amrfinderplus_virulence_report
    String amrfinderplus_amr_genes = amrfinderplus_nuc.amrfinderplus_amr_genes
    String amrfinderplus_stress_genes = amrfinderplus_nuc.amrfinderplus_stress_genes
    String amrfinderplus_virulence_genes = amrfinderplus_nuc.amrfinderplus_virulence_genes
    String amrfinderplus_amr_classes = amrfinderplus_nuc.amrfinderplus_amr_classes
    String amrfinderplus_amr_subclasses = amrfinderplus_nuc.amrfinderplus_amr_subclasses
    }
 }
