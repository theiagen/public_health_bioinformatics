version 1.0

task retrieve_aligned_contig_paf {
  meta {
    description: "Parse minimap2 PAF file and return aligned contigs in FASTA format"
  }
  input {
    File paf
    File assembly
    String samplename
    String docker = "staphb/seqkit:2.3.1"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    # retrieve contig name
    cut -f1 "~{paf}" > contig_names.txt
    
    # extract mapped contigs in FASTA format
    seqkit grep -f contig_names.txt "~{assembly}" > "~{samplename}".fasta

    echo "Original contig number:"
    cat ~{assembly} | grep ">" | wc -l 
    
    echo "Filtered contig number:"
    cat "~{samplename}".fasta | grep ">" | wc -l 

    if [ ! -s "~{samplename}".fasta ]; then
      echo "Filtered assembly file is empty! Removing..."
      rm "~{samplename}".fasta
    fi
  >>>
  output {
    File final_assembly = "~{samplename}.fasta"
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}

task calculate_coverage_paf {
  meta {
    description: "Parse minimap2 PAF file and return the breadth of coverage"
  }
  input {
    File paf
    String docker = "quay.io/quay/ubuntu"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    # get reference length from paf file
    ref_len=$(head -n1 ~{paf} | cut -f7 | tee REFERENCE_LENGTH)

    # get sum of unique contigs that align to reference
    contig_len=$(cat ~{paf} | cut -f2 | uniq | awk '{ sum += $1 } END { print sum }' | tee CONTIG_LENGTH)

    result=$(awk "BEGIN { printf \"%.2f\", ($contig_len / $ref_len) * 100 }")
    echo $result | tee PERCENT_COVERAGE
  >>>
  output {
    String reference_length = read_string("REFERENCE_LENGTH")
    String contig_length = read_string("CONTIG_LENGTH")
    String percent_coverage = read_string("PERCENT_COVERAGE")
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task sam_to_sorted_bam {
  meta {
    description: "Converts SAM file to sorted BAM file"
  }
  input {
    File sam
    String samplename
    String docker = "quay.io/staphb/samtools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    # Convert SAM to BAM, and sort it based on read name
    samtools view -Sb ~{sam} > "~{samplename}".bam
    samtools sort "~{samplename}".bam -o "~{samplename}".sorted.bam

    # index sorted BAM
    samtools index "~{samplename}".sorted.bam > "~{samplename}".sorted.bam.bai
  >>>
  output {
    File bam = "~{samplename}.sorted.bam"
    File bai = "~{samplename}.sorted.bam.bai"
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task retrieve_pe_reads_bam {
  meta {
    description: "Parse minimap2 SAM file and return paired-end reads in FASTQ format"
  }
  input {
    File bam
    String samplename
    Int sam_flag = "4" # unmapped reads (SAM flag 4)
    String docker = "quay.io/staphb/samtools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    # Convert to fastq.gz, discarding the singleton reads
    samtools fastq -f "~{sam_flag}" -1 "~{samplename}"_1.fq.gz -2 "~{samplename}"_2.fq.gz \
    -s singleton.fq.gz "~{bam}"
  >>>
  output {
    File read1 = "~{samplename}_1.fq.gz"
    File read2 = "~{samplename}_2.fq.gz"
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task calculate_coverage {
  meta {
    description: "Use bedtools to calculate average depth of coverage from BAM file and a reference file"
  }
  input {
    File bam
    File bai
    String docker = "quay.io/staphb/bedtools:2.31.0"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<

    bedtools genomecov -d -ibam ~{bam} > coverage_stats.txt

    # get reference length from paf file
    ref_len=$(wc -l coverage_stats.txt | cut -f1 -d" " | tee REFERENCE_LENGTH)

    # get sum of coverage per position
    sum_coverage=$(cat coverage_stats.txt | cut -f3 | awk '{ sum += $1 } END { print sum }' | tee SUM_COVERAGE)

    # calculate the mean depth coverage
    result=$(awk "BEGIN { printf \"%.2f\", $sum_coverage / $ref_len }")
    echo $result | tee AVG_DEPTH_COVERAGE
  >>>
  output {
    File stats = "coverage_stats.txt"
    String ref_len = read_string("REFERENCE_LENGTH")
    String sum_coverage = read_string("SUM_COVERAGE")
    String mean_depth_coverage = read_string("AVG_DEPTH_COVERAGE")
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

