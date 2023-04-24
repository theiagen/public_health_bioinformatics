version 1.0

task bowtie_retrieve_mapped_pe {
  input {
    File read1
    File read2
    String samplename
    File reference
    String docker = "staphb/bowtie2:2.5.1"
    Int disk_size = 100
    Int cpu = 4
  }
  command <<<
    # date and version control
    date | tee DATE

    # Build a Bowtie2 index for the reference genome
    bowtie2-build "~{reference}" "~{reference}"

    # Map paired-end reads to reference genome using bowtie2
    #bowtie2 -x "~{reference}" -1 "~{read1}" -2 "~{read2}" -S mapped_reads.sam --no-unal
    bowtie2 -x "~{reference}" -1 "~{read1}" -2 "~{read2}" p 4 --un-conc-gz "unmapped.fastq.gz" -S "mapped_reads.sam" --no-mixed --no-discordant

    # Convert SAM file to BAM file using samtools
    samtools view -Sb mapped_reads.sam > mapped_reads.bam

    # Sort BAM file by coordinates using samtools
    samtools sort mapped_reads.bam -o mapped_reads.sorted.bam

    # Index sorted BAM file using samtools
    samtools index mapped_reads.sorted.bam

    # Retrieve all mapped reads using samtools
    samtools view -h -b -F 4 mapped_reads.sorted.bam | samtools fastq -1 "~{samplename}_mapped_1.fastq.gz" -2 "~{samplename}_mapped_2.fastq.gz" - 
    
  >>>
  output {
    File read1_mapped = "~{samplename}_mapped_1.fastq.gz"
    File read2_mapped = "~{samplename}_mapped_2.fastq.gz"
  }
  runtime {
      docker: "~{docker}"
      memory: "8 GB"
      cpu: cpu
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB" # TES
      preemptible: 0
      maxRetries: 0
  }
}