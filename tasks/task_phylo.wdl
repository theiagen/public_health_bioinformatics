version 1.0

task ksnp3 {
  input {
    Array[File] assembly_fasta
    Array[String] samplename
    String cluster_name
    Int kmer_size = 19
    String docker_image = "quay.io/staphb/ksnp3:3.1"
    Int mem_size_gb = 8
    Int CPUs = 4
  }
  command <<<
  assembly_array=(~{sep=' ' assembly_fasta})
  assembly_array_len=$(echo "${#assembly_array[@]}")
  samplename_array=(~{sep=' ' samplename})
  samplename_array_len=$(echo "${#samplename_array[@]}")
  
  # Ensure assembly, and samplename arrays are of equal length
  if [ "$assembly_array_len" -ne "$samplename_array_len" ]; then
    echo "Assembly array (length: $assembly_array_len) and samplename array (length: $samplename_array_len) are of unequal length." >&2
    exit 1
  fi

  # create file of filenames for kSNP3 input
  touch ksnp3_input.tsv
  for index in ${!assembly_array[@]}; do
    assembly=${assembly_array[$index]}
    samplename=${samplename_array[$index]}
    
    echo -e "${assembly}\t${samplename}" >> ksnp3_input.tsv
  done
  # run ksnp3 on input assemblies
  kSNP3 -in ksnp3_input.tsv -outdir ksnp3 -k ~{kmer_size} -core 
  
  # rename ksnp3 outputs with cluster name 
  mv ksnp3/core_SNPs_matrix.fasta ~{cluster_name}_core_SNPs_matrix.fasta
  mv ksnp3/tree.core.tre ~{cluster_name}_core.tree

  >>>
  output {
    File ksnp3_matrix = "${cluster_name}_core_SNPs_matrix.fasta"
    File ksnp3_tree = "${cluster_name}_core.tree"
    String ksnp3_docker_image = docker_image
  }
  runtime {
    docker: docker_image
    memory: "~{mem_size_gb} GB"
    cpu: CPUs
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
task snp_dists {
  input {
    File alignment
    String cluster_name
  }
  command{
    # date and version control
    date | tee DATE
    snp-dists -v | tee VERSION

    snp-dists ${alignment} > ${cluster_name}_snp_distance_matrix.tsv
}
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File snp_matrix = "${cluster_name}_snp_distance_matrix.tsv"
  }
  runtime {
    docker: "quay.io/staphb/snp-dists:0.6.2"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
