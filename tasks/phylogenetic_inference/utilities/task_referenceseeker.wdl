version 1.0

task referenceseeker {
  input {
    File assembly_fasta
    String samplename
    File referenceseeker_db = "gs://theiagen-public-resources-rp/reference_data/databases/referenceseeker/referenceseeker-bacteria-refseq-205.v20210406.tar.gz"
    Float referenceseeker_ani_threshold = 0.95
    Float referenceseeker_conserved_dna_threshold = 0.69
    Int disk_size = 200 
    Int cpu = 4
    Int memory = 16
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/referenceseeker:1.8.0--pyhdfd78af_0"
  }
  command <<<
    # get version
    referenceseeker --version | tee VERSION

    # get DB and uncompress
    echo "Uncompressing referenceseeker DB..."
    mkdir db
    tar -C ./db/ -xzf  ~{referenceseeker_db} --strip-components 1

    # run referenceseeker
    echo "Running referenceseeker..."
    referenceseeker \
      --threads ~{cpu} \
      --conserved-dna ~{referenceseeker_conserved_dna_threshold} \
      --ani ~{referenceseeker_ani_threshold} \
      ./db/ ~{assembly_fasta} > ~{samplename}_referenceseeker.tsv

    # grab accession from top hit: remove header line, grab the first line, cut for first column 
    # containing the accessions and remove the version number 
    tail -n +2 ~{samplename}_referenceseeker.tsv | head -n 1 | cut -f 1 | cut -d '.' -f 1 | tee REFERENCE
  >>>
  output {
    String referenceseeker_top_hit_ncbi_accession = read_string("REFERENCE")
    String referenceseeker_version = read_string("VERSION")
    File referenceseeker_tsv = "~{samplename}_referenceseeker.tsv"
    String referenceseeker_docker = docker
    String referenceseeker_database = referenceseeker_db
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: memory + " GB"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}