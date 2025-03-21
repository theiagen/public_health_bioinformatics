version 1.0

task hicap {
  meta {
    description: "cap locus typing for H. influenzae assembly into serotypes a-f"
  }
  input {
    File assembly
    String samplename
    Float min_gene_percent_coverage = 0.80 #Minimum percentage coverage to consider a single gene complete. [default: 0.80]
    Float min_gene_percent_identity = 0.70 #Minimum percentage identity to consider a single gene complete. [default: 0.70]
    Float min_broken_gene_percent_identity = 0.80 #Minimum percentage identity to consider a broken gene. [default: 0.80]
    Int broken_gene_length = 60 #Minimum length to consider a broken gene. [default: 60]
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/hicap:1.0.3--py_0"
    Int cpu = 2
    Int memory = 8
    Int disk_size = 50
  }
  command <<<
    echo $(hicap --version 2>&1) | sed 's/^hicap //' | tee VERSION

    mkdir output_dir

    hicap \
      -q ~{assembly} \
      -o output_dir \
      --gene_coverage ~{min_gene_percent_coverage} \
      --gene_identity ~{min_gene_percent_identity} \
      --broken_gene_length ~{broken_gene_length} \
      --broken_gene_identity ~{min_broken_gene_percent_identity} \
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
