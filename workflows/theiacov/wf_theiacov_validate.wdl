version 1.0

import "../tasks/task_validate.wdl" as validation
import "../tasks/task_versioning.wdl" as versioning

workflow theiacov_validate {
  input {
    String terra_project
    String terra_workspace
    String datatable1
    String datatable2
    String out_dir = "./"
    String out_prefix = "VALIDATION"
  }
  call validation.export_two_tsvs {
    input:
      terra_project = terra_project,
      terra_workspace = terra_workspace,
      datatable1 = datatable1,
      datatable2 = datatable2
  }
  call validation.compare_two_tsvs {
    input:
      datatable1_tsv = export_two_tsvs.datatable1_tsv,
      datatable2_tsv = export_two_tsvs.datatable2_tsv,
      out_dir = out_dir,
      out_prefix = out_prefix
  }
  call versioning.version_capture {
    input:
  }
  output {
    String theiacov_validation_version = version_capture.phvg_version
    String theiacov_validation_date = version_capture.date
    File theiacov_validation_report_pdf = compare_two_tsvs.pdf_report
    File theiacov_validation_report_xl = compare_two_tsvs.xl_report
  }
}