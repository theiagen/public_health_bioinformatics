version 1.0

task mycosnptree {
  input {
    Array[File] assembly_fasta
    Array[String] samplename
    String docker="quay.io/theiagen/mycosnp:dev"
    String strain="B11205"
    String accession="GCA_016772135"
  }
  command <<<
    date | tee DATE
    echo $(nextflow pull rpetit3/mycosnp-nf 2>&1) | sed 's/^.*revision: //;' | tee MYCOSNPTREE_VERSION

    assembly_array=(~{sep=' ' assembly_fasta})
    assembly_array_len=$(echo "${#assembly_array[@]}")
    samplename_array=(~{sep=' ' samplename})
    samplename_array_len=$(echo "${#samplename_array[@]}")

    # Ensure assembly, and samplename arrays are of equal length
    if [ "$assembly_array_len" -ne "$samplename_array_len" ]; then
      echo "Assembly array (length: $assembly_array_len) and samplename array (length: $samplename_array_len) are of unequal length." >&2
      exit 1
    fi

    # Make sample FOFN
    echo "sample,fasta" > samples.csv
    for index in ${!assembly_array[@]}; do
      assembly=${assembly_array[$index]}
      samplename=${samplename_array[$index]}
      echo -e "${samplename},${assembly}" >> samples.csv
    done

    # Run MycoSNP
    mkdir mycosnptree
    cd mycosnptree
    if nextflow run rpetit3/mycosnp-nf -entry NFCORE_MYCOSNPTREE --input ../samples.csv --fasta /reference/~{accession}/masked/reference-consensus.fa --publish_dir_mode copy --rapidnj False --fasttree False --iqtree; then
      # Everything finished, pack up the results and clean up
      find work/ -name "*.iqtree" | xargs -I {} cp {} ./
      rm -rf .nextflow/ work/
      cd ..
      tar -cf - mycosnptree/ | gzip -n --best  > mycosnptree.tar.gz
    else
      # Run failed
      exit 1
    fi
  >>>
  output {
    String mycosnptree_version = read_string("MYCOSNPTREE_VERSION")
    String mycosnptree_docker = docker
    String analysis_date = read_string("DATE")
    String reference_strain = strain
    String reference_accession = accession
    File mycosnptree_tree = "mycosnptree/results/combined/phylogeny/iqtree/alignment.fasta.treefile"
    File mycosnptree_iqtree_log = "mycosnptree/alignment.fasta.iqtree"
    File mycosnptree_full_results = "mycosnptree.tar.gz"
  }
  runtime {
    docker: "~{docker}"
    memory: "32 GB"
    cpu: 4
    disks:  "local-disk 50 SSD"
    maxRetries: 3
    preemptible: 0
  }
}
