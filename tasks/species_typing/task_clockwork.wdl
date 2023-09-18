version 1.0

task clockwork_decon_reads {
  # Inputs
  input {
    File read1
    File? read2
    String samplename
    Int disk_size = 200
    Int cpu = 16
    Int mem = 64
    }

  command <<<
    # Print and save date
    date | tee DATE

    # Print and save version
    clockwork version > VERSION

    # Map reads to the clockwork reference
    clockwork map_reads \
    --unsorted_sam ~{samplename} /varpipe_wgs/tools/clockwork-0.11.3/OUT/ref.fa \
    "~{samplename}.sam" \
    ~{read1} \
    ~{read2}

    # Remove contaminants (reads that map with high identity to non-MTB sequences)
    clockwork remove_contam \
    /varpipe_wgs/tools/clockwork-0.11.3/OUT/remove_contam_metadata.tsv \
    "~{samplename}.sam" \
    "~{samplename}_outfile_read_counts" \
    "./clockwork_cleaned_~{samplename}_R1.fastq.gz" \
    "./clockwork_cleaned_~{samplename}_R2.fastq.gz"

  >>>
  output {
    File clockwork_cleaned_read1 = "./clockwork_cleaned_~{samplename}_R1.fastq.gz"
    File clockwork_cleaned_read2 = "./clockwork_cleaned_~{samplename}_R2.fastq.gz"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/cdcgov/varpipe_wgs_with_refs:2bc7234074bd53d9e92a1048b0485763cd9bbf6f4d12d5a1cc82bfec8ca7d75e"
    memory: "~{mem} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3 
  }
}