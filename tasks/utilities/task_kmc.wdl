version 1.0

task kmc {
  meta {
    description: "Estimate genome size using KMC (https://github.com/refresh-bio/KMC)"
  }
  input {
    File read1
    String docker = "quay.io/biocontainers/kmc:3.2.1--h9ee0642_0"
    Int disk_size = 100
    Int cpu = 4 
  }
  command <<<
    kmc | head -n 1 | tee VERSION

    # run kmc
    kmc \
      -sm \
      -m8 \
      -t"~{cpu}" \
      -k21 \
      -ci10 \
      ~{read1} \
      kmc.res \ 
      . \
      > LOG

    cat LOG | grep "unique counted" | tr -s ' ' | cut -d ' ' -f8 > UNIQUE_COUNTED

  >>>
  output {
    String genome_size = read_string("UNIQUE_COUNTED") 
    String kmc_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
