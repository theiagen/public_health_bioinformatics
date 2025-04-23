version 1.0

task megahit_pe {
  input {
    File read1
    File read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/megahit:1.2.9"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 16
    String? kmers
    String? megahit_opts
  }
  command <<<
    # fail hard
    set -euo pipefail

    # get version
    megahit --version | tee VERSION

    # get the memory required, assuming its input is GB
    memory=$(python3 -c "print(~{memory} * 1000000000)")

    # run megathit
    megahit \
      -1 ~{read1} \
      -2 ~{read2} \
      ~{'--k-list ' + kmers} \
      -m ${memory} \
      -t ~{cpu} \
      -o megahit/ \
      ~{megahit_opts}

    mv megahit/final.contigs.fa ~{samplename}_megahit.fasta

  >>>
  output {
    File assembly_fasta = "~{samplename}_megahit.fasta"
    String megahit_version = read_string("VERSION")
    String megahit_docker = '~{docker}'
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