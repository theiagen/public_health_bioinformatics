version 1.0

task centroid {
  input {
    Array[File] assembly_fasta
    Int disk_size = 50
    Int cpu = 1
    Int memory = 4
  }
  command <<<
    mkdir INPUT_DIR
    ln -s ~{sep=' ' assembly_fasta} INPUT_DIR
    
    # centroid.py expects a positional argument with a path to a directory of FASTA files
    centroid.py INPUT_DIR/

    # rename the centroid genome so it can be accessed outside this task
    mv -v $(ls */* | grep $(cat centroid_out.txt)) centroid.fasta

  >>>
  output {
    String centroid_genome_fasta_filename = read_string("centroid_out.txt")
    File centroid_genome_fasta_file = "centroid.fasta"
    File centroid_mash_tsv = "mash-results.tsv"
  }
  runtime {
    # hardcoding docker as I do not expect updates to occur to Centroid in the future
    docker: "quay.io/theiagen/centroid:latest"
    cpu: cpu
    memory: memory + " GB"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}