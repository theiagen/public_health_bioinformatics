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
    seqkit grep -f contig_names.txt "~{assembly}" > "~{samplename}"_mapped_contigs.fasta

  >>>
  output {
    File parse_paf_contigs = "~{samplename}_mapped_contigs.fasta"
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

task retrieve_unaligned_pe_reads_sam {
  meta {
    description: "Parse minimap2 SAM file and return unaligned paired-end reads in FASTQ format"
  }
  input {
    File sam
    String samplename
    String docker = "staphb/samtools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    # Convert SAM to BAM
    samtools view -Sb "~{sam}" > mapped_reads.bam

    # Sort the BAM file
    samtools sort -o mapped_reads.sorted.bam mapped_reads.bam
    
    # Index the sorted BAM file
    samtools index mapped_reads.sorted.bam
    
    # Extract the unmapped reads
    samtools view -b -f 4 -F 264 mapped_reads.sorted.bam > unmapped_reads.bam

  >>>
  output {
    File unmapped_bam = "unmapped_reads.bam"
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

task bam2fastq {
  meta {
    description: "Convert bam to FASTQ file with bedtools "
  }
  input {
    File bam
    String samplename
    String docker = "staphb/bedtools:2.30.0"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<

    # Convert unmapped reads to FASTQ
    bedtools bamtofastq -i "~{bam}" -fq "~{samplename}"_1.fq -fq2 "~{samplename}"_2.fq

    # Compress fq files
    gzip "~{samplename}"_1.fq "~{samplename}"_2.fq
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