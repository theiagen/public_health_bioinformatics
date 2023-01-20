version 1.0

task snp_dists {
  input {
    File alignment
    String cluster_name
    Int disk_size = 100
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
    String version = read_string("VERSION")
    File snp_matrix = "~{cluster_name}_snp_distance_matrix.tsv"
  }
  runtime {
    docker: "quay.io/staphb/snp-dists:0.8.2"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}

task reorder_matrix {
  input {
    File tree
    File matrix
    String cluster_name
    Int disk_size = 100
  }
  command <<<
    python3 <<CODE
    from Bio import Phylo
    import pandas as pd

    # read in newick tree
    tree = Phylo.read("~{tree}", "newick")
    
    # extract ordered terminal ends
    term_names = [term.name for term in tree.get_terminals()]

    # read in matrix into pandas data frame
    snps = pd.read_csv("~{matrix}", header=0, index_col=0, delimiter="\t")

    # ensure all header and index values are strings for proper reindexing
    # this is because if sample_name is entirely composed of integers, pandas 
    # auto-casts them as integers; get_terminals() interprets those as strings. 
    # this incompatibility leads to failure and an empty ordered SNP matrix
    snps.columns = snps.columns.astype(str)
    snps.index = snps.index.astype(str)

    # reorder matrix according to terminal ends
    snps = snps.reindex(index=term_names, columns=term_names)

    # add phandango suffix to ensure continuous coloring
    snps_out1 = snps.add_suffix(":c1")

    # write out reordered matrix to a file
    snps_out1.to_csv("~{cluster_name}_ordered_snp_distance_matrix.csv", sep=",")

    # reroot tree with midpoint
    tree.root_at_midpoint()

    # re-extract ordered terminal ends of rerooted tree
    term_names = [term.name for term in tree.get_terminals()]

    # reorder matrix with re-ordered terminal ends
    snps = snps.reindex(index=term_names, columns=term_names)

    # add phandango suffix to ensure continuous coloring
    snps_out2 = snps.add_suffix(":c1")

    # write out reordered matrix of rerooted tree to a file
    snps_out2.to_csv("~{cluster_name}_midpoint_snp_distance_matrix.csv", sep=",")

    # write rerooted tree to a file
    Phylo.write(tree, "~{cluster_name}_midpoint_tree.nwk", "newick")

    CODE
  >>>
  output{
    File ordered_matrix = "~{cluster_name}_ordered_snp_distance_matrix.csv"
    File ordered_midpoint_matrix = "~{cluster_name}_midpoint_snp_distance_matrix.csv"
    File midpoint_rooted_tree = "~{cluster_name}_midpoint_tree.nwk"
  }
  runtime {
    docker: "staphb/mykrobe:0.12.1" # used because it contains both biopython and pandas
    memory: "2 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
   # maxRetries: 3
    preemptible: 0
  }

}