version 1.0

task ksnp4 {
  input {
    Array[File] assembly_fasta
    Array[String] samplename
    String cluster_name
    Int kmer_size = 19
    String ksnp4_args = "" # add -ML to calculate a maximum likelihood tree or -NJ to calculate a neighbor-joining tree
    String docker_image = "us-docker.pkg.dev/general-theiagen/staphb/ksnp4:4.1"
    File? previous_ksnp4_snps
    Int memory = 16
    Int cpu = 2
    Int disk_size = 100
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
    # ensure kSNP file naming convention is met
    assembly_renamed_array=()
    for index in ${!assembly_array[@]}; do
        assembly=${assembly_array[$index]}
        # ensure kSNP file naming convention is met by removing non-id, dot-separated info,
        # e.g. sample01.ivar.consensus.fasta will be renamed to sample01.fasta
        assembly_renamed=$(echo $assembly | sed 's/\.\(.*\)\././')
        echo "ASSEMBLY: $assembly"
        echo "ASSEMBLY_RENAMED: $assembly_renamed"
        mv $assembly $assembly_renamed
        assembly_renamed_array+=($assembly_renamed)
    done

    # create file of filenames for ksnp4 input
    touch ksnp4_input.tsv
    for index in ${!assembly_renamed_array[@]}; do
      assembly=${assembly_renamed_array[$index]}
      samplename=${samplename_array[$index]}
      echo -e "${assembly}\t${samplename}" >> ksnp4_input.tsv
    done

    echo "ksnp4_input.tsv:: "
    cat ksnp4_input.tsv

    # run ksnp4 on input assemblies
    kSNP4 \
      -in ksnp4_input.tsv \
      -outdir ksnp4 \
      -k ~{kmer_size} \
      -core \
      -vcf \
      ~{'-SNPs_all ' + previous_ksnp4_snps} \
      -CPU ~{cpu} \
      ~{ksnp4_args}

    # rename ksnp4 outputs with cluster name
    # sometimes the core nwk and fasta outputs do not have content
    mv -v ksnp4/core_SNPs_matrix.fasta ksnp4/~{cluster_name}_core_SNPs_matrix.fasta
    mv -v ksnp4/tree.core_SNPs.parsimony.tre ksnp4/~{cluster_name}_core.nwk # note the name of this file has changed with kSNP4 (used to be "tree.core.tre")

    # to allow appending new genomes to an existing tree, save the SNPs_all file
    mv -v ksnp4/SNPs_all ksnp4/~{cluster_name}_SNPs_all

    if [ -s ksnp4/~{cluster_name}_core_SNPs_matrix.fasta ]; then # is the file not-empty?
      echo "The core SNP matrix was produced" | tee SKIP_SNP_DIST # then do NOT skip
    else
      echo "The core SNP matrix could not be produced" | tee SKIP_SNP_DIST # otherwise, skip
    fi

    # capture sample name of genome used as reference
    ls ksnp4/*.vcf | cut -d '.' -f 2 | tee ksnp4_VCF_REF_SAMPLENAME.txt

    # capture number of SNPs
    cut -f 2 -d ':' ksnp4/COUNT_SNPs | tr -d ' ' > ksnp4/NUMBER_SNPS
    # capture number of core SNPs
    grep "Number core SNPs: " ksnp4/COUNT_coreSNPs | cut -f 2 -d ':' | tr -d ' ' > ksnp4/NUMBER_CORE_SNPS

    # rename the 2 vcf files by appending ~{cluster_name} and removing the ref genome name to make final filenames predictable
    mv -v ksnp4/VCF.*.vcf ksnp4/~{cluster_name}_VCF.reference_genome.vcf
    mv -v ksnp4/VCF.SNPsNotinRef.* ksnp4/~{cluster_name}_VCF_.SNPsNotinRef.tsv

    mv -v ksnp4/SNPs_all_matrix.fasta ksnp4/~{cluster_name}_pan_SNPs_matrix.fasta
    mv -v ksnp4/tree.parsimony.tre ksnp4/~{cluster_name}_pan_parsimony.nwk

    if [ -f ksnp4/tree.ML.tre ]; then
      mv -v ksnp4/tree.ML.tre ksnp4/~{cluster_name}_ML.nwk
    fi
    if [ -f ksnp4/tree.NJ.tre ]; then
      mv -v ksnp4/tree.NJ.tre ksnp4/~{cluster_name}_NJ.nwk
    fi

  >>>
  output {
    File ksnp4_core_matrix = "ksnp4/~{cluster_name}_core_SNPs_matrix.fasta"
    File ksnp4_core_tree = "ksnp4/~{cluster_name}_core.nwk"
    File ksnp4_vcf_ref_genome = "ksnp4/~{cluster_name}_VCF.reference_genome.vcf"
    File ksnp4_vcf_snps_not_in_ref = "ksnp4/~{cluster_name}_VCF_.SNPsNotinRef.tsv"
    String ksnp4_vcf_ref_samplename = read_string("ksnp4_VCF_REF_SAMPLENAME.txt")
    File ksnp4_pan_matrix = "ksnp4/~{cluster_name}_pan_SNPs_matrix.fasta"
    File ksnp4_pan_parsimony_tree = "ksnp4/~{cluster_name}_pan_parsimony.nwk"
    File? ksnp4_ml_tree = "ksnp4/~{cluster_name}_ML.nwk"
    File? ksnp4_nj_tree = "ksnp4/~{cluster_name}_NJ.nwk"
    String ksnp4_number_snps = read_string("ksnp4/NUMBER_SNPS")
    String ksnp4_number_core_snps = read_string("ksnp4/NUMBER_CORE_SNPS")
    File ksnp4_snps_all = "ksnp4/~{cluster_name}_SNPs_all"
    File ksnp4_input = "ksnp4_input.tsv"
    String skip_core_snp_dists = read_string("SKIP_SNP_DIST")
    Array[File] ksnp_outs = glob("ksnp4/*")
    String ksnp4_docker_image = docker_image
  }
  runtime {
    docker: docker_image
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 0
    maxRetries: 0
  }
}
