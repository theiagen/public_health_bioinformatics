version 1.0

task gubbins {
  input {
    File alignment
    String cluster_name
    String docker = "sangerpathogens/gubbins"
    Int? filter_percent = 25 #default is 25%
    Int? iterations = 5
    String? first_tree_nuc_subst_model = "GTRGAMMA"
    Array[String]? first_tree_args
    String? tree_builder = "raxml"
    Array[String]? tree_args
    String? nuc_subst_model = "GTRGAMMA"
    Array[String]? nuc_subst_model_args = ["4"]
    String? best_nuc_subst_model
    Array[String]? outgroup
    Int? bootstrap = 0
    File? dates_file
  }
  command <<<
    # date and version control
    date | tee DATE
    run_gubbins.py --version | tee VERSION

    run_gubbins.py \
    ~{alignment} \
    --prefix ~{cluster_name} \
    --first-tree-builder fasttree \
    --filter-percentage ~{filter_percent} \
    --iterations ~{iterations} \
    --first-model ~{first_tree_nuc_subst_model} \
    ~{'--first-tree-args ' + first_tree_args} \
    --tree-builder ~{tree_builder} \
    ~{'--tree-args ' + tree_args} \
    ~{'--best-model ' + best_nuc_subst_model} \
    --model ~{nuc_subst_model} \
    --model-args ~{sep="," nuc_subst_model_args} \
    --bootstrap ~{bootstrap} \
    ~{'--outgroup ' + outgroup} \
    ~{'--date ' + dates_file} \
    --threads 2
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File gubbins_final_tree = "~{cluster_name}.final_tree.tre"
    File gubbins_final_labelled_tree = "~{cluster_name}.node_labelled.final_tree.tre"
    File gubbins_polymorphic_fasta = "~{cluster_name}.filtered_polymorphic_sites.fasta"
    File gubbins_recombination_gff = "~{cluster_name}.recombination_predictions.gff"
    File gubbins_branch_stats = "~{cluster_name}.per_branch_statistics.csv"
    File? gubbins_timetree = "~{cluster_name}.final_tree.timetree.tre"
    File? gubbins_timetree_stats = "~{cluster_name}.lsd.out"
  }
  runtime {
    docker: "~{docker}"
    memory: "32 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}
