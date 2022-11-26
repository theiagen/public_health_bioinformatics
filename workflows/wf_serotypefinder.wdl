version 1.0

import "../tasks/task_taxon_id.wdl" as taxon_ID
import "../tasks/task_versioning.wdl" as versioning

workflow serotypefinder {
    input {
        String  samplename
        File    ecoli_assembly
    }
    call taxon_ID.serotypefinder_one_sample {
    input:
        samplename = samplename,
        ecoli_assembly = ecoli_assembly
    }
    call versioning.version_capture{
      input:
    }
    output {
      String  serotypefinder_wf_version = version_capture.phbg_version
      String  serotypefinder_wf_analysis_date = version_capture.date
      
      String serotypefinder_report  = serotypefinder_one_sample.serotypefinder_report
      String serotypefinder_docker  = serotypefinder_one_sample.serotypefinder_docker
      String serotypefinder_serotype = serotypefinder_one_sample.serotypefinder_serotype
    }
}
