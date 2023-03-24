version 1.0

import "../../tasks/species_typing/task_serotypefinder.wdl" as serotypefinder_file
import "../../tasks/task_versioning.wdl" as versioning

workflow serotypefinder {
    input {
        String  samplename
        File    ecoli_assembly
    }
    call serotypefinder_file.serotypefinder as serotypefinder_task {
    input:
        samplename = samplename,
        assembly = ecoli_assembly
    }
    call versioning.version_capture{
      input:
    }
    output {
      String  serotypefinder_wf_version = version_capture.phb_version
      String  serotypefinder_wf_analysis_date = version_capture.date
      
      File serotypefinder_report  = serotypefinder_task.serotypefinder_report
      String serotypefinder_docker  = serotypefinder_task.serotypefinder_docker
      String serotypefinder_serotype = serotypefinder_task.serotypefinder_serotype
    }
}
