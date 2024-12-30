version 1.0

task multiqc {
  input {
    Array[File] qc_files
    Int memory = 4
    Int cpu = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/multiqc:1.25"
  }
  command <<<  
    # Run MultiQC on the provided QC files
    multiqc \
      --outdir . \
      --force \
      ~{sep=' ' qc_files}
  >>>
  output {
    File multiqc_html = "multiqc_report.html"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}