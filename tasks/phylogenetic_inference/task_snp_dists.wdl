version 1.0

task snp_dists {
  input {
    File alignment
    String cluster_name
    Int disk_size = 50
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/snp-dists:0.8.2"
  }
  command <<<
    # date and version control
    date | tee DATE
    snp-dists -v | tee VERSION

    # create snp-dists matrix file
    snp-dists ~{alignment} > ~{cluster_name}_snp_distance_matrix.tsv
  >>>
  output {
    String date = read_string("DATE")
    String snp_dists_version = read_string("VERSION")
    String snp_dists_docker = docker
    File snp_matrix = "~{cluster_name}_snp_distance_matrix.tsv"
  }
  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}