version 1.0

import "../tasks/tools/task_pmga.wdl" as pmga
import "../tasks/task_versioning.wdl" as versioning

workflow pmga_wf {
    input {
        File assembly
        String samplename
    }
    call pmga.pmga {
        input:
            assembly = assembly,
            samplename = samplename
    }
    call versioning.version_capture{
        input:
    }
    output {
        String pmga_wf_version = version_capture.phbg_version
        String pmga_wf_analysis_date = version_capture.date
        String pmga_version = pmga.version
        String pmga_docker = pmga.docker
        String pmga_speciesdb = pmga.pmga_speciesdb
        String pmga_serotype = pmga.pmga_serotype
        String pmga_genes = pmga.pmga_genes
        String pmga_notes = pmga.pmga_notes
        File pmga_results = pmga.pmga_results
        File pmga_allele_matrix = pmga.pmga_allele_matrix
        File pmga_blast_final = pmga.pmga_blast_final
        File pmga_blast_raw = pmga.pmga_blast_raw
        File pmga_loci_counts = pmga.pmga_loci_counts
        File pmga_gff = pmga.pmga_gff
    }
}
