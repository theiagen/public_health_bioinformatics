version 1.0

task snippy_variants {
  input {
    File reference_genome_file
    File? assembly_fasta
    File? read1
    File? read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/snippy:4.6.0"
    Int cpu = 8
    Int memory = 32
    Int disk_size = 100
    # Paramters 
    # --map_qual: Minimum read mapping quality to consider (default '60')
    # --base_quality: Minimum base quality to consider (default '13')
    # --min_coverage: Minimum site depth to for calling alleles (default '10') 
    # --min_frac: Minumum proportion for variant evidence (0=AUTO) (default '0')
    # --min_quality: Minumum QUALITY in VCF column 6 (default '100')
    # --maxsoft: Maximum soft clipping to allow (default '10')
    Int? map_qual
    Int? base_quality
    Int min_coverage = 10
    Float? min_frac
    Int? min_quality
    Int? maxsoft
  }
  command <<<
    snippy --version | head -1 | tee VERSION

    # set input variable
    if [ -f "~{assembly_fasta}" ]; then
      reads="--ctgs ~{assembly_fasta}"
    elif [ -f "~{read1}" ] && [ -z "~{read2}" ]; then
      reads="--se ~{read1}"
    elif [ -f "~{read1}" ] && [ -f "~{read2}" ]; then
      reads="--R1 ~{read1} --R2 ~{read2}"
    else
      echo "ERROR: No reads or assembly provided"
      exit 1
    fi
    
    # call snippy
    snippy \
      --reference ~{reference_genome_file} \
      --outdir ~{samplename} \
      ${reads} \
      --cpus ~{cpu} \
      --ram ~{memory} \
      --prefix ~{samplename} \
      ~{'--mapqual ' + map_qual} \
      ~{'--basequal ' + base_quality} \
      ~{'--mincov ' + min_coverage} \
      ~{'--minfrac ' + min_frac} \
      ~{'--minqual ' + min_quality} \
      ~{'--maxsoft ' + maxsoft}

    # Compress output dir
    tar -cvzf "./~{samplename}_snippy_variants_outdir.tar" "./~{samplename}"

    # compute number of reads aligned to reference
    samtools view -c "~{samplename}/~{samplename}.bam" > READS_ALIGNED_TO_REFERENCE

    # create coverage stats file
    samtools coverage "~{samplename}/~{samplename}.bam" -o "~{samplename}/~{samplename}_coverage.tsv"

    # capture number of variants from summary file
    grep "VariantTotal" ~{samplename}/~{samplename}.txt | cut -f 2 > VARIANTS_TOTAL

    # # compute proportion of reference genome with depth >= min_coverage
    # compute read depth at every position in the genome
    samtools depth -a "~{samplename}/~{samplename}.bam" -o "~{samplename}/~{samplename}_depth.tsv"

    # compute reference genome length
    reference_length=$(cat "~{samplename}/~{samplename}_depth.tsv" | wc -l)
    echo $reference_length | tee REFERENCE_LENGTH

    # filter depth file to only include positions with depth >= min_coverage
    awk -F "\t" -v cov_var=~{min_coverage} '{ if ($3 >= cov_var) print;}' "~{samplename}/~{samplename}_depth.tsv" > "~{samplename}/~{samplename}_depth_~{min_coverage}.tsv"
    
    # compute proportion of genome with depth >= min_coverage
    reference_length_passed_depth=$(cat "~{samplename}/~{samplename}_depth_~{min_coverage}.tsv" | wc -l)
    echo $reference_length_passed_depth | tee REFERENCE_LENGTH_PASSED_DEPTH

    # check if reference_length is equal to 0, if so, output a warning
    if [ "$reference_length" -eq 0 ]; then
      echo "Could not compute percent reference coverage: reference length is 0" > PERCENT_REF_COVERAGE
    else
      # compute percent reference coverage
      echo $reference_length_passed_depth $reference_length | awk '{ print ($1/$2)*100 }' > PERCENT_REF_COVERAGE
    fi

  >>>
  output {
    String snippy_variants_version = read_string("VERSION")
    String snippy_variants_docker = docker
    File snippy_variants_reference_genome = "~{reference_genome_file}"
    File snippy_variants_outdir_tarball = "./~{samplename}_snippy_variants_outdir.tar"
    Array[File] snippy_variants_outputs = glob("~{samplename}/~{samplename}*")
    File snippy_variants_results = "~{samplename}/~{samplename}.csv"
    File snippy_variants_bam = "~{samplename}/~{samplename}.bam"
    File snippy_variants_bai ="~{samplename}/~{samplename}.bam.bai"
    File snippy_variants_summary = "~{samplename}/~{samplename}.txt"
    String snippy_variants_num_reads_aligned = read_string("READS_ALIGNED_TO_REFERENCE")
    File snippy_variants_coverage_tsv = "~{samplename}/~{samplename}_coverage.tsv"
    String snippy_variants_num_variants = read_string("VARIANTS_TOTAL")
    File snippy_variants_depth = "~{samplename}/~{samplename}_depth.tsv"
    String snippy_variants_ref_length = read_string("REFERENCE_LENGTH")
    String snippy_variants_ref_length_passed_depth = read_string("REFERENCE_LENGTH_PASSED_DEPTH")
    String snippy_variants_percent_ref_coverage = read_string("PERCENT_REF_COVERAGE")
  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: "~{cpu}"
      disks: "local-disk " + disk_size + " SSD"
      preemptible: 0
      maxRetries: 3
  }
}