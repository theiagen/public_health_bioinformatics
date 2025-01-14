version 1.0

task hicap {
  meta {
    description: "cap locus typing for H. influenzae assembly into serotypes a-f"
  }
  input {
    File assembly
    String samplename

    Boolean output_full_sequence = false
    Float gene_coverage = 0.80
    Float gene_identity = 0.70
    Float broken_gene_identity = 0.80
    Int broken_gene_length = 60

    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/hicap:1.0.3--py_0"
    Int cpu = 2
    Int memory = 8
    Int disk_size = 50

    #Parameters
    #-q QUERY_FP, --query_fp QUERY_FP              Input FASTA query
    #-o OUTPUT_DIR, --output_dir OUTPUT_DIR        Output directory
    #-d DATABASE_DIR, --database_dir DATABASE_DIR  Directory containing locus database. [default: /usr/local/lib/python3.6/site-
    #                                              packages/hicap/database]
    #-m MODEL_FP, --model_fp MODEL_FP              Path to prodigal model. [default: /usr/local/lib/python3.6/site-
    #                                              packages/hicap/model/prodigal_hi.bin]
    #-s, --full_sequence                           Write the full input sequence out to the genbank file rather than just the region
    #                                              surrounding and including the locus
    #--gene_coverage GENE_COVERAGE                 Minimum percentage coverage to consider a single gene complete. [default: 0.80]
    #--gene_identity GENE_IDENTITY                 Minimum percentage identity to consider a single gene complete. [default: 0.70]
    #--broken_gene_length BROKEN_GENE_LENGTH       Minimum length to consider a broken gene. [default: 60]
    #--broken_gene_identity BROKEN_GENE_IDENTITY   Minimum percentage identity to consider a broken gene. [default: 0.80]
    #--threads THREADS                             Threads to use for BLAST+. [default: 1]
    #--log_fp LOG_FP                               Record logging messages to file
    #--debug                                       Print debug messages
  }
  command <<<
    echo $(hicap --version 2>&1) | sed 's/^hicap //' | tee VERSION

    mkdir output_dir

    hicap \
      -q ~{assembly} \
      -o output_dir \
      --gene_coverage ~{gene_coverage} \
      --gene_identity ~{gene_identity} \
      --broken_gene_length ~{broken_gene_length} \
      --broken_gene_identity ~{broken_gene_identity} \
      ~{true='--full_sequence' false='' output_full_sequence} \
      --threads ~{cpu}
      
    filename=$(basename ~{assembly})

    # if there are no hits for a cap locus, no file is produced
    if [ ! -f ./output_dir/${filename%.*}.tsv ]; then
      echo "No_hits" > hicap_serotype
      echo "No_hits" > hicap_genes
      touch ~{samplename}.hicap.tsv
    else
      tail -n1 ./output_dir/${filename%.*}.tsv | cut -f 2 > hicap_serotype
      tail -n1 ./output_dir/${filename%.*}.tsv | cut -f 4 > hicap_genes
      mv ./output_dir/${filename%.*}.tsv ~{samplename}.hicap.tsv
    fi
  >>>
  output {
    String hicap_serotype = read_string("hicap_serotype")
    String hicap_genes = read_string("hicap_genes")
    File hicap_results_tsv = "~{samplename}.hicap.tsv"
    String hicap_version = read_string("VERSION")
    String hicap_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu    
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
