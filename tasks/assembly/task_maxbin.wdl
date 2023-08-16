version 1.0

task maxbin2 {
  input {
    File assembly
    File read1
    File read2
    String samplename
    Int min_contig_length = 1000
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/maxbin2:2.2.7--hdbdd923_5"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 16
    String? maxbin2_opts
  }
  command <<<
    run_MaxBin.pl -v | head -1 | cut -d ' ' -f 2 | tee VERSION
    run_MaxBin.pl \
      -reads ~{read1} \
      -reads2 ~{read2} \
      -min_contig_length ~{min_contig_length}
      -thread ~{cpu} \
      -contig ~{assembly} \
      -out ~{samplename} \
      ~{maxbin2_opts}

    # compress bins into zip file
    tar -cvf ~{samplename}_bins.tar.gz ~{samplename}.*.fasta
  >>>
  output {
    File maxbin2_bins_fasta = "~{samplename}_bins.tar.gz"
    File maxbin2_summary = "~{samplename}.summary"
    String maxbin2_version = read_string("VERSION")
    String maxbin2_docker = '~{docker}'
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}

