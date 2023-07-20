version 1.0

task kmc {
  meta {
    description: "Estimate genome size using KMC (https://github.com/refresh-bio/KMC)"
  }
  input {
    File read1
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/kmc:3.2.1--h9ee0642_0"
    Int disk_size = 100
    Int cpu = 8
    Int mem = 32
    Int kmer_length = 21
    Int min_kmer_count = 10
  }
  command <<<
    kmc | head -n 1 | tee VERSION

    # run kmc
    # kmc [options] <input_file> <output_file> <working_dir>
    # -sm - uses strict memory mode (memory from -m<size> switch will not be exceeded)
    # -m<size> - the max amount of RAM in GB (default: 12) 
    # -t<value> - total number of threads to use
    # -k<len> - k-mer length (default: 25)
    # -ci<value> - exclude k-mers occuring less than <value> times (default: 1e9)
    kmc \
      -sm \
      -m"~{mem}" \
      -t"~{cpu}" \
      -k"~{kmer_length}" \
      -ci"~{min_kmer_count}" \
      ~{read1} \
      kmc_outputs \
      . \
      > LOG

    # kmc_outputs is a mess of files that are not human readable
    # however, the stdout does produce some useful stats. 
    #  the no. of unique counted k-mers can be used as an estimate of genome size
    grep "unique counted k" LOG | tr -s ' ' | cut -d ' ' -f8 > UNIQUE_COUNTED

    # extracting only the kmer statistics and writing to file:
    tail -n8 LOG > ~{samplename}_kmer_stats.txt
  >>>
  output {
    String est_genome_size = read_string("UNIQUE_COUNTED") 
    File kmer_stats = "~{samplename}_kmer_stats.txt"
    String kmc_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
