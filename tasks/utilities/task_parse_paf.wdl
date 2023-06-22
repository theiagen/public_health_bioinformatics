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

    #if [ ! -s "~{samplename}".fasta ]; then
    #  echo "Filtered assembly file is empty! Removing..."
    #  rm "~{samplename}".fasta
    #fi
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
    # Convert SAM to BAM, and sort it based on read name
    samtools sort -n "~{sam}" -O BAM > mapped_name_sorted_reads.bam

    # Convert unmapped reads (SAM flag 4) to fastq.gz, discarding the singleton reads
    samtools fastq -f 4 -1 "~{samplename}"_1.fq.gz -2 "~{samplename}"_2.fq.gz \
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