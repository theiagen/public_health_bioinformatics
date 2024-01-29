version 1.0

task semibin {
  input {
    File sorted_bam
    File sorted_bai
    String samplename
    File assembly_fasta
    String environment = "global"
    Int cpu = 6
    Int mem = 8
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/semibin:2.0.2--pyhdfd78af_0"
  }
  command <<<
    # date and version control
    date | tee DATE
    SemiBin -v | tee SEMIBIN_VERSION

    # run SemiBin
    SemiBin single_easy_bin \
      -i ~{assembly_fasta} \
      -b ~{sorted_bam} \
      -o ~{samplename} \
      -t ~{cpu} \
      --environment ~{environment}
  >>>
  output {
    String semibin_version = read_string("SEMIBIN_VERSION")
    Array[File] semibin_bins = glob("~{samplename}/output_recluster_bins/bin.*.fa")
  }
  runtime {
    docker: docker
    memory: mem + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}