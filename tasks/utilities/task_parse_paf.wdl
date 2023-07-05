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

    apt-get install bc

    # get reference length from paf file
    ref_len=$(head -n1 ~{paf} | cut -f7 | tee REFERENCE_LENGTH)

    # get sum of contig lengths that align to reference
    contig_len=$(cat ~{paf} | cut -f2 | awk '{ sum += $1 } END { print sum }' | tee CONTIG_LENGTH)

    result=$(echo "scale=2; ($contig_len / $ref_len)*100" | bc)
    
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

task retrieve_pe_reads_sam {
  meta {
    description: "Parse minimap2 SAM file and return paired-end reads in FASTQ format"
  }
  input {
    File sam
    String samplename
    Int sam_flag = "4" # unmapped reads (SAM flag 4)
    String docker = "staphb/samtools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    # Convert SAM to BAM, and sort it based on read name
    samtools sort -n "~{sam}" -O BAM > mapped_name_sorted_reads.bam

    # Convert to fastq.gz, discarding the singleton reads
    samtools fastq -f "~{sam_flag}" -1 "~{samplename}"_1.fq.gz -2 "~{samplename}"_2.fq.gz \
    -s singleton.fq.gz mapped_name_sorted_reads.bam
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

