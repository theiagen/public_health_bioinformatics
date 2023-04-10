version 1.0

task gubbins {
  input {
    File alignment
    String cluster_name
    String docker = "quay.io/biocontainers/gubbins:3.3--py310pl5321h8472f5a_0"
    Int filter_percent = 25 #default is 25%
    Int iterations = 5
    String tree_builder = "raxml"
    String? tree_args
    String nuc_subst_model = "GTRGAMMA"
    Int bootstrap = 0
    String? outgroup
    File? dates_file
    Int cpu = 4
    Int disk_size = 100
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
    ~{'--model ' + nuc_subst_model} \
    --bootstrap ~{bootstrap} \
    ~{'--outgroup ' + outgroup} \
    ~{'--date ' + dates_file} \
    --threads ~{cpu}

    # rename newick files, TSV file, and text file to have matching and appropriate file endings
    mv -v ~{cluster_name}.node_labelled.final_tree.tre ~{cluster_name}.node_labelled.final_tree.nwk
    mv -v ~{cluster_name}.per_branch_statistics.csv ~{cluster_name}.per_branch_statistics.tsv
    mv -v ~{cluster_name}.final_tree.tre ~{cluster_name}.final_tree.nwk
    
    # rename optional output files to have appropriate file endings
    # echo if file not found to prevent throwing error codes
    mv -v ~{cluster_name}.final_tree.timetree.tre ~{cluster_name}.final_tree.timetree.nwk || echo "No timetree .nwk file created, skipping renaming"
    mv -v ~{cluster_name}.lsd.out ~{cluster_name}.lsd.out.txt || echo "No treetime .lsd.out file created, skipping renaming"

  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    String gubbins_docker = docker
    File gubbins_final_tree = "~{cluster_name}.final_tree.nwk"
    File gubbins_final_labelled_tree = "~{cluster_name}.node_labelled.final_tree.nwk"
    File gubbins_polymorphic_fasta = "~{cluster_name}.filtered_polymorphic_sites.fasta"
    File gubbins_recombination_gff = "~{cluster_name}.recombination_predictions.gff"
    File gubbins_branch_stats = "~{cluster_name}.per_branch_statistics.tsv"
    File? gubbins_timetree = "~{cluster_name}.final_tree.timetree.nwk"
    File? gubbins_timetree_stats = "~{cluster_name}.lsd.out.txt"
  }
  runtime {
    docker: "~{docker}"
    memory: "32 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
