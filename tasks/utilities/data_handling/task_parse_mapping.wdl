version 1.0

task retrieve_aligned_contig_paf {
  meta {
    description: "Parse minimap2 PAF file and return aligned contigs in FASTA format"
  }
  input {
    File paf
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/seqkit:2.4.0--h9ee0642_0"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
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
    memory: memory + " GB"
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
    String docker = "us-docker.pkg.dev/general-theiagen/quay/ubuntu:latest"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
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
    memory: memory + " GB"
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
    Int? min_qual
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    # Samtools verion capture
    samtools --version | head -n1 | cut -d' ' -f2 | tee VERSION

    # Convert SAM to BAM, and sort it based on read name
    samtools view ~{"-q " + min_qual} -Sb ~{sam} > "~{samplename}".bam
    samtools sort "~{samplename}".bam -o "~{samplename}".sorted.bam

    # index sorted BAM
    samtools index "~{samplename}".sorted.bam > "~{samplename}".sorted.bam.bai
  >>>
  output {
    File bam = "~{samplename}.sorted.bam"
    File bai = "~{samplename}.sorted.bam.bai"
    String samtools_version = read_string("VERSION")
    String samtools_docker = "~{docker}"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task bam_to_fastq {
  meta {
    description: "Convert BAM file to FASTQ equivalent to BWA task"
  }
  input {
    File bam
    String samplename
    Boolean paired = false
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8   
  }
  command <<<
    # Samtools verion capture
    samtools --version | head -n1 | cut -d' ' -f2 | tee VERSION

    samtools view \
      -@ ~{cpu} \
      -F 0x904 \
      -b \
      -o ~{samplename}.sorted.aligned-reads.bam \
      ~{bam}

    # convert SAM to BAM that only includes unaligned reads
    samtools view \
      -@ ~{cpu} \
      -f 4 \
      -b \
      -o ~{samplename}.sorted.unaligned-reads.bam \
      ~{bam}

    # Get aligned and unaligned bams
    if ~{paired}; then
      echo -e "\nGenerating FASTQs for aligned paired-end reads"
      samtools fastq \
        -@ ~{cpu} \
        -F 4 \
        -1 ~{samplename}_R1.fastq.gz \
        -2 ~{samplename}_R2.fastq.gz \
        ~{samplename}.sorted.aligned-reads.bam
      echo "Generating FASTQs for unaligned paired-end reads"
      # note the lowercase 'f' here is imporant
      samtools fastq \
        -@ ~{cpu} \
        -f 4 \
        -1 ~{samplename}_unaligned_R1.fastq.gz \
        -2 ~{samplename}_unaligned_R2.fastq.gz \
        ~{samplename}.sorted.unaligned-reads.bam
    else
      echo -e "\nGenerating FASTQs for aligned single-end reads\n"
      samtools fastq \
        -@ ~{cpu} \
        -F 4 \
        -0 ~{samplename}_R1.fastq.gz \
        ~{samplename}.sorted.aligned-reads.bam
      echo -e "Generating FASTQs for unaligned single-end reads\n"
      # again, lowercase 'f' is important for getting all unaligned reads
      samtools fastq \
        -@ ~{cpu} \
        -f 4 \
        -0 ~{samplename}_unaligned_R1.fastq.gz \
        ~{samplename}.sorted.unaligned-reads.bam
    fi
  >>>
  output {
    String sam_version = read_string("VERSION")
    File read1_aligned = "~{samplename}_R1.fastq.gz"
    File? read2_aligned = "~{samplename}_R2.fastq.gz"
    File read1_unaligned = "~{samplename}_unaligned_R1.fastq.gz"
    File? read2_unaligned = "~{samplename}_unaligned_R2.fastq.gz"
    File sorted_bam_aligned = "~{samplename}.sorted.aligned-reads.bam"
    File sorted_bam_unaligned = "~{samplename}.sorted.unaligned-reads.bam"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
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
    String prefix = ""
    Int sam_flag = "4" # unmapped reads (SAM flag 4)
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    # Convert to fastq.gz, discarding the singleton reads
    samtools fastq -f "~{sam_flag}" -1 "~{prefix}"_"~{samplename}"_1.fq.gz -2 "~{prefix}"_"~{samplename}"_2.fq.gz \
    -s singleton.fq.gz "~{bam}"
  >>>
  output {
    File read1 = "~{prefix}_~{samplename}_1.fq.gz"
    File read2 = "~{prefix}_~{samplename}_2.fq.gz"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
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
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bedtools:2.31.0"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    # get version
    bedtools --version | cut -d' ' -f2 | tee VERSION

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
    String bedtools_version = read_string("VERSION")
    String bedtools_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task mask_low_coverage {
  meta {
    description: "Use bedtools to create a coverage mask for low coverage regions in a BAM file"
  }
  input {
    File bam
    File bai
    File reference_fasta
    Int min_depth
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bedtools:2.31.0"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    set -euo pipefail

    # get version
    bedtools --version | cut -d' ' -f2 | tee VERSION

    # make a temporary reference fasta and make sure the header is ONLY the ">" + reference accession number.
    # bedtools maskfasta requires the header to be the same as the reference name in the bam file
    awk '{if (/^>/) print $1; else print}' ~{reference_fasta} > mod_reference.fasta

    # report depth at all regions of a genome including regions with 0 coverage
    bedtools genomecov -bga -ibam ~{bam} > all_coverage_regions.bed

    # filter out regions that have a coverage greater than {min_depth}. Only need low coverage regions here.
    awk -v min_depth=~{min_depth} '$4 < min_depth' all_coverage_regions.bed > low_coverage_regions.bed

    # mask the low coverage regions in the reference fasta file
    bedtools maskfasta -fi mod_reference.fasta -bed low_coverage_regions.bed \
      -fo masked_reference.fasta
  >>>
  output {
    File low_coverage_regions_bed = "low_coverage_regions.bed"
    File all_coverage_regions_bed = "all_coverage_regions.bed"
    File mask_reference_fasta = "masked_reference.fasta"
    String bedtools_version = read_string("VERSION")
    String bedtools_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task assembled_reads_percent {
  input {
    File bam
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    # Count the total number of reads
    total_reads=$(samtools view -c ~{bam})

    # Count the number of mapped reads
    mapped_reads=$(samtools view -F 4 -c ~{bam})

    # Calculate the percentage of mapped reads using awk
    percentage_mapped=$(awk -v mapped=$mapped_reads -v total=$total_reads 'BEGIN { print (mapped / total) * 100 }')

    echo $percentage_mapped | tee PERCENTAGE_MAPPED
  >>>
  output {
    String percentage_mapped = read_string("PERCENTAGE_MAPPED")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}
