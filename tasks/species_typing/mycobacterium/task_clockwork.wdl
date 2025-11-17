version 1.0

task clockwork_decon_reads {
  input {
    File read1
    File? read2 # only optional to not fail in merlin_magic
    String samplename
    Int disk_size = 200
    Int cpu = 16
    Int memory = 64
    String docker = "ashedpotatoes/clockwork-plus:v0.12.5.3-CRyPTIC"
  }
  command <<<
    # Print and save version
    clockwork version > VERSION

    # Map reads to the clockwork reference
    clockwork map_reads \
      --unsorted_sam ~{samplename} /ref/Ref.remove_contam/ref.fa \
      "~{samplename}.sam" \
      ~{read1} \
      ~{read2}

    samtools sort -n  "~{samplename}.sam" > ~{samplename}.sorted.sam
    
    # Remove contaminants (reads that map with high identity to non-MTB sequences)
    clockwork remove_contam \
      /ref/Ref.remove_contam/remove_contam_metadata.tsv \
      "~{samplename}.sorted.sam" \
      "~{samplename}_outfile_read_counts" \
      "clockwork_cleaned_~{samplename}_R1.fastq.gz" \
      "clockwork_cleaned_~{samplename}_R2.fastq.gz"

    # Clean up files
    rm "~{samplename}.sam"
  >>>
  output {
    File clockwork_cleaned_read1 = "clockwork_cleaned_~{samplename}_R1.fastq.gz"
    File clockwork_cleaned_read2 = "clockwork_cleaned_~{samplename}_R2.fastq.gz"
    String clockwork_version = read_string("VERSION")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 1
  }
}