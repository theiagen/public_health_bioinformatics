version 1.0

task lyvset {
  input {
    Array[File] read1
    Array[File] read2
    File reference_genome
    String dataset_name
    String docker_image = "quay.io/staphb/lyveset:1.1.4f"
    Int memory = 16
    Int cpu = 4
    Int disk_size = 100

    # Lyve-SET Parameters
    ##COMMON OPTIONS
    ##--allowedFlanking  0              allowed flanking distance in bp.
    ##                                  Nucleotides this close together cannot be
    ##                                  considered as high-quality.
    ##--min_alt_frac     0.75           The percent consensus that needs
    ##                                  to be reached before a SNP is called.
    ##                                  Otherwise, 'N'
    ##--min_coverage     10             Minimum coverage needed before a
    ##                                  SNP is called. Otherwise, 'N'
    ##--presets          ""             See presets.conf for more information
    ##--numcpus          1              number of cpus
    ##PERFORM CERTAIN STEPS
    ##--mask-phages                  Search for and mask phages in the reference genome
    ##--mask-cliffs                  Search for and mask 'Cliffs' in pileups
    ##
    ##SKIP CERTAIN STEPS
    ##--nomatrix                     Do not create an hqSNP matrix
    ##--nomsa                        Do not make a multiple sequence alignment
    ##--notrees                      Do not make phylogenies
    ##--singleend                    Treat everything like single-end. Useful
    ##                               for when you think there is a single-
    ##                               end/paired-end bias.
    ##OTHER SHORTCUTS
    ##--fast                         Shorthand for --downsample --mapper snap --nomask-phages
    ##                                             --nomask-cliffs --sample-sites
    ##--downsample                   Downsample all reads to 50x. Approximated according
    ##                               to the ref genome assembly
    ##--sample-sites                 Randomly choose a genome and find SNPs in a quick
    ##                               and dirty way. Then on the SNP-calling stage,
    ##                               only interrogate those sites for SNPs for each
    ##                               genome (including the randomly-sampled genome).
    ##
    ##MODULES
    ##--read_cleaner none            Which read cleaner? Choices: none, CGP, BayesHammer
    ##--mapper       smalt           Which mapper? Choices: smalt, snap
    ##--snpcaller    varscan         Which SNP caller? Choices: varscan, vcftools
    
    Int allowedFlanking = 0
    Float min_alt_frac = 0.75
    Int min_coverage = 10
    String? presets
    Boolean mask_phages = false
    Boolean mask_cliffs = false
    Boolean nomatrix = false
    Boolean nomsa = false
    Boolean notrees = false
    Boolean singleend = false
    Boolean fast = false
    Boolean downsample = false
    Boolean sample_sites = false
    String read_cleaner = "none"
    String mapper = "smalt"
    String snpcaller = "varscan"

  }
  command <<<
    date | tee DATE

    read1_array=(~{sep=' ' read1})
    read1_array_len=$(echo "${#read1[@]}")
    read2_array=(~{sep=' ' read2})
    read2_array_len=$(echo "${#read2[@]}")

    # Ensure read arrays are of equal length
    if [ "$read1_array_len" -ne "$read2_index_array_len" ]; then
      echo "read1 array (length: $read1_array_len) and read2 index array (length: $read2_array_len) are of unequal length." >&2
      exit 1
    fi

    echo "before symlinks:"
    ls ./



    echo "after symlinks:"
    ls ./

    # create lyvset project
    set_manage.pl --create ~{dataset_name}
    #shuffle paired end reads
    echo "YOUUO"
    for index in ${!read1_array[@]}; do
      shuffleSplitReads.pl --numcpus ~{cpu} -o ./interleaved ${read1_array[$index]} ${read2_array[$index]}
    done
    echo "YUP!"
    # then moved into your project dir
    mv ./interleaved/*.fastq.gz ~{dataset_name}/reads/
    # cleanup
    rmdir interleaved
    mkdir ~{dataset_name}/ref/
    cp ~{reference_genome} ~{dataset_name}/ref/reference.fasta
    launch_set.pl --numcpus 8 -ref ~{dataset_name}/ref/reference.fasta ~{dataset_name}

  >>>
  output {
    
    String lyveset_docker_image = docker_image
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
