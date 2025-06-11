version 1.0

task megahit {
  input {
    File read1
    File? read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/megahit:1.2.9"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 16
    Int min_contig_length = 1
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

    if [ -n "~{read2}" ]; then
      # if read2 is provided, use paired-end mode
      megahit \
        -1 ~{read1} \
        -2 ~{read2} \
        ~{'--k-list ' + kmers} \
        --min-contig-len ~{min_contig_length} \
        -m ${memory} \
        -t ~{cpu} \
        -o megahit/ \
        ~{megahit_opts}
    else
      # if read2 is not provided, use single-end mode
      megahit \
        -r ~{read1} \
        ~{'--k-list ' + kmers} \
        --min-contig-len ~{min_contig_length} \
        -m ${memory} \
        -t ~{cpu} \
        -o megahit/ \
        ~{megahit_opts}
    fi

    if [ ! -s megahit/final.contigs.fa ]; then
      echo "DEBUG: Megahit assembly failed. No output generated."
      echo "FAIL" > STATUS
    else
      echo "DEBUG: Megahit assembly completed successfully."
      echo "PASS" > STATUS
    fi

    mv megahit/final.contigs.fa ~{samplename}_megahit_contigs.fasta

  >>>
  output {
    File? assembly_fasta = "~{samplename}_megahit_contigs.fasta"
    String megahit_status = read_string("STATUS")
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