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
    Int map_qual = 60
    Int base_quality = 13
    Int min_coverage = 10
    Float min_frac = 0
    Int min_quality = 100
    Int maxsoft = 10
  }
  command <<<
    # set -euo pipefail to avoid silent failure
    set -euo pipefail

    snippy --version | head -1 | tee VERSION

    # set input variable
    if [ -f "~{assembly_fasta}" ]; then
      echo "DEBUG: Using assembly fasta file"
      reads="--ctgs ~{assembly_fasta}"
    elif [ -f "~{read1}" ] && [ -z "~{read2}" ]; then
      echo "DEBUG: Using single-end read file"
      reads="--se ~{read1}"
    elif [ -f "~{read1}" ] && [ -f "~{read2}" ]; then
      echo "DEBUG: Using paired-end read files"
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

    # compute number of reads aligned to reference (excluding unmapped reads)
    samtools view -c -F 4 "~{samplename}/~{samplename}.bam" > READS_ALIGNED_TO_REFERENCE

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
      echo $reference_length_passed_depth $reference_length | awk '{ printf("%.2f", ($1/$2)*100) }' > PERCENT_REF_COVERAGE
    fi

    # Compute percentage of reads aligned
    reads_aligned=$(cat READS_ALIGNED_TO_REFERENCE)
    total_reads=$(samtools view -c "~{samplename}/~{samplename}.bam")
    echo $total_reads > TOTAL_READS
    if [ "$total_reads" -eq 0 ]; then
      echo "Could not compute percent reads aligned: total reads is 0" > PERCENT_READS_ALIGNED
    else
      echo $reads_aligned $total_reads | awk '{ printf("%.2f", ($1/$2)*100) }' > PERCENT_READS_ALIGNED
    fi

    # Create QC metrics file
    line_count=$(wc -l < "~{samplename}/~{samplename}_coverage.tsv")
    # Check the number of lines in the coverage file, to consider scenarios e.g. for V. cholerae that has two chromosomes and therefore coverage metrics per chromosome
    if [ "$line_count" -eq 2 ]; then
      head -n 1 "~{samplename}/~{samplename}_coverage.tsv" | tr ' ' '\t' > COVERAGE_HEADER
      sed -n '2p' "~{samplename}/~{samplename}_coverage.tsv" | tr ' ' '\t' > COVERAGE_VALUES
    elif [ "$line_count" -gt 2 ]; then
      # Multiple chromosomes (header + multiple data lines)
      header=$(head -n 1 "~{samplename}/~{samplename}_coverage.tsv")
      output_header=""
      output_values=""
      # while loop to iterate over each line in the coverage file
      while read -r line; do
        if [ -z "$output_header" ]; then
          output_header="$header"
          output_values="$line"
        else
          output_header="$output_header\t$header"
          output_values="$output_values\t$line"
        fi
      done < <(tail -n +2 "~{samplename}/~{samplename}_coverage.tsv")
      echo "$output_header" | tr ' ' '\t' > COVERAGE_HEADER
      echo "$output_values" | tr ' ' '\t' > COVERAGE_VALUES
    else
      # Coverage file has insufficient data
      echo "Coverage file has insufficient data." > COVERAGE_HEADER
      echo "" > COVERAGE_VALUES
    fi

    # Build the QC metrics file
    echo -e "samplename\treads_aligned_to_reference\ttotal_reads\tpercent_reads_aligned\tvariants_total\tpercent_ref_coverage\t$(cat COVERAGE_HEADER)" > "~{samplename}/~{samplename}_qc_metrics.tsv"
    echo -e "~{samplename}\t$reads_aligned\t$total_reads\t$(cat PERCENT_READS_ALIGNED)\t$(cat VARIANTS_TOTAL)\t$(cat PERCENT_REF_COVERAGE)\t$(cat COVERAGE_VALUES)" >> "~{samplename}/~{samplename}_qc_metrics.tsv"

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
    File snippy_variants_qc_metrics = "~{samplename}/~{samplename}_qc_metrics.tsv"
    String snippy_variants_percent_reads_aligned = read_string("PERCENT_READS_ALIGNED")
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