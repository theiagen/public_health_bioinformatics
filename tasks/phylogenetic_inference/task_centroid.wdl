version 1.0

task centroid {
  input {
    Array[File] assembly_fasta
    Int disk_size = 50
    Int cpu = 1
    Int memory = 4
    String docker = "quay.io/theiagen/centroid:0.1.0"
  }
  command <<<
    centroid.py --version | tee CENTROID_VER
    
    mkdir INPUT_DIR
    # copy all asms to a single directory
    cp -v ~{sep=' ' assembly_fasta} INPUT_DIR
    
    # centroid.py expects a positional argument with a path to a directory of FASTA files
    centroid.py INPUT_DIR/
    
    # set bash variable which ONLY has the filename, e.g. "SAMN19774644_contigs.fasta"
    CENTROID_FASTA_FILENAME=$(cat centroid_out.txt)

    # rename the centroid genome so it can be accessed outside this task
    mv -v INPUT_DIR/"$CENTROID_FASTA_FILENAME" centroid.fasta

    # capture samplename of centroid genome for use downstream
    cut -d '_' -f 1 centroid_out.txt | tee CENTROID_GENOME_SAMPLENAME
  >>>
  output {
    String centroid_genome_samplename = read_string("CENTROID_GENOME_SAMPLENAME")
    File centroid_genome_fasta = "centroid.fasta"
    File centroid_mash_tsv = "mash-results.tsv"
    String centroid_docker = docker
    String centroid_version = read_string("CENTROID_VER")
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