version 1.0

task gubbins {
  input {
    File alignment
    String cluster_name
    String docker = "sangerpathogens/gubbins"
  }
  command <<<
    # date and version control
    date | tee DATE
    run_gubbins.py --version | tee VERSION

    run_gubbins.py \
    ~{alignment} \
    --prefix ~{cluster_name} \
    --first-tree-builder fasttree
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File gubbins_final_tree = "~{cluster_name}.final_tree.tre"
    File gubbins_final_labelled_tree = "~{cluster_name}.node_labelled.final_tree.tre"
    File gubbins_polymorphic_fasta = "~{cluster_name}.filtered_polymorphic_sites.fasta"
    File gubbins_recombination_gff = "~{cluster_name}.recombination_predictions.gff"
    File gubbins_branch_stats = "~{cluster_name}.per_branch_statistics.csv"
  }
  runtime {
    docker: "~{docker}"
    memory: "32 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}
