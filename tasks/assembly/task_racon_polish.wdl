version 1.0

task racon {
  input {
    File unpolished_fasta     # Assembly (contigs)
    File reads                # Sequencing reads
    File alignments           # Minimap2 alignments (PAF/SAM)
    Int cpu = 4               
    Int memory = 16           
    Int disk_size = 100       
    String samplename         # Output file name prefix
    String docker = "staphb/racon:1.4.20"
  }

  command <<< 
    # Capture version for reproducibility
    racon --version | tee VERSION

    # Run Racon for polishing
    racon \
      -t ~{cpu} \
      ~{reads} \
      ~{alignments} \
      ~{unpolished_fasta} \
      > ~{samplename}.polished.fasta
  >>>

  output {
    File polished_fasta = "~{samplename}.polished.fasta" 
    String racon_version = read_string("VERSION")        
    String racon_docker = "~{docker}"                  
  }

  runtime {
    docker: "~{docker}"         
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
