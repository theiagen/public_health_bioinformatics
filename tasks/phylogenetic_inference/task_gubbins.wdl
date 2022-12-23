version 1.0

task gubbins {
  input {
    File alignment
    String cluster_name
    String docker = "sangerpathogens/gubbins"
    Int? filter_percent = 25 #default is 25%
    Int? iterations = 5
    String? tree_builder = "raxml"
    String? tree_args
    String? nuc_subst_model = "GTRCAT"
    Boolean? best_nuc_subst_model = false
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
    --filter-percentage ~{filter_percent} \
    --iterations ~{iterations} \
    --tree-builder ~{tree_builder} \
    ~{'--tree-args ' + tree_args} \
    ~{true="--best-model" false="" best_nuc_subst_model} \
    ~{'--model ' + nuc_subst_model} \
    --bootstrap ~{bootstrap} \
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
