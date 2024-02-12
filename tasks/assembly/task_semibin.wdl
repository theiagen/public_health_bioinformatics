version 1.0

task semibin {
  input {
    File sorted_bam
    File sorted_bai
    String samplename
    File assembly_fasta
    String environment = "global"
    Int min_length = 1000
    Float ratio = 0.05
    Int cpu = 6
    Int memory = 8
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/semibin:2.0.2--pyhdfd78af_0"
  }
  command <<<
    # date and version control
    date | tee DATE
    SemiBin -v | tee SEMIBIN_VERSION

    # check the number of contigs are greater than min_len
    count=$(awk '/^>/ {if (length(seq) > ~{min_length}) count++; seq = ""; next} {seq = seq $0} END {if (length(seq) > ~{min_length}) count++; print count}' ~{assembly_fasta})
    echo "Number of contigs greater than ~{min_length} characters: $count"

    # SeminBin2 requires at least two contigs greater than the min_len to run
    if [ $count -gt 1 ]; then
        echo "Running SemiBin2"
        
        # run SemiBin
        SemiBin single_easy_bin \
          -i ~{assembly_fasta} \
          -b ~{sorted_bam} \
          -o ~{samplename} \
          -t ~{cpu} \
          --environment ~{environment} \
          --min-len ~{min_length} \
          --ratio ~{ratio}
    else
        echo "One or fewer contigs found with ~{min_length} in ~{assembly_fasta}."
        echo "Aborting binning with SemiBin2..."
    fi

  >>>
  output {
    String semibin_version = read_string("SEMIBIN_VERSION")
    String semibin_docker = docker
    Array[File] semibin_bins = glob("~{samplename}/output_recluster_bins/bin.*.fa")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}