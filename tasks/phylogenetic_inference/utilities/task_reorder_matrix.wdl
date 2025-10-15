version 1.0

task reorder_matrix {
  input {
    File input_tree
    File matrix
    String cluster_name

    String? outgroup_root # will preferentially root with outgroup if defined
    Boolean? midpoint_root_tree
    
    Int disk_size = 100
    Int cpu = 1
    Int memory = 2
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/mykrobe:0.12.1" # used because it contains both biopython and pandas
    Boolean phandango_coloring = false
  }
  command <<<
    # removing any "_contigs" suffixes from the tree and matrix
    sed 's/_contigs//g' ~{input_tree} > temporary_tree.nwk
    sed 's/_contigs//g' ~{matrix} > temporary_matrix.tsv

    python3 <<CODE
    from Bio import Phylo
    import pandas as pd
    import os

    # read in newick tree
    tree = Phylo.read("temporary_tree.nwk", "newick")
    
    # read in matrix into pandas data frame
    snps = pd.read_csv("temporary_matrix.tsv", header=0, index_col=0, delimiter="\t")

    # ensure all header and index values are strings for proper reindexing
    # this is because if sample_name is entirely composed of integers, pandas 
    # auto-casts them as integers; get_terminals() interprets those as strings. 
    # this incompatibility leads to failure and an empty ordered SNP matrix
    snps.columns = snps.columns.astype(str)
    snps.index = snps.index.astype(str)

    # reroot tree with midpoint, if midpoint_root_tree is set to true
    if "~{outgroup_root}":
        tree.root_with_outgroup(["~{outgroup_root}"])
    elif "~{midpoint_root_tree}" == "true":
        tree.root_at_midpoint()

    # extract ordered terminal ends of tree (could be midpoint rooted or not, depending on midpoint_root optional input)
    term_names = [term.name for term in tree.get_terminals()]

    # reorder matrix according to the order of tree terminal ends
    snps = snps.reindex(index=term_names, columns=term_names)

    # add phandango suffix to ensure continuous coloring
    if ("~{phandango_coloring}" == "true"):
        snps_out2 = snps.add_suffix(":c1")
    else:
        snps_out2 = snps

    # remove header from index column
    snps_out2.index.name = ""

    # write out the reordered matrix to a file
    snps_out2.to_csv("~{cluster_name}_distance_matrix.csv", sep=",")

    # write tree to a file (same as input tree if not midpoint rooted)
    Phylo.write(tree, "~{cluster_name}_tree.nwk", "newick")

    CODE
  >>>
  output {
    File ordered_matrix = "~{cluster_name}_distance_matrix.csv"
    File tree = "~{cluster_name}_tree.nwk"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
