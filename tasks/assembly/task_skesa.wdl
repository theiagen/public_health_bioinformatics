version 1.0

task skesa {
  input {
    File read1
    File? read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/skesa:2.4.0"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 16
    Int min_contig_length = 1
    String? skesa_opts
  }
  command <<<
    # Get skesa version
    skesa --version 2>&1 | grep "SKESA" | sed -E 's/^.*SKESA ([0-9.]+).*/\1/' | tee VERSION
    
    # Use paired-end mode if read2 is provided, otherwise single-end
    if [ -n "~{read2}" ]; then
        skesa --gz \
        --fastq ~{read1},~{read2} \
        --use_paired_ends \
        --contigs_out ~{samplename}_skesa_contigs.fasta \
        --min_contig ~{min_contig_length} \
        --memory ~{memory} \
        --cores ~{cpu} \
        --vector_percent 1 \
        ~{skesa_opts}
    else
      skesa --fastq ~{read1} \
        --contigs_out ~{samplename}_skesa_contigs.fasta \
        --min_contig ~{min_contig_length} \
        --memory ~{memory} \
        --cores ~{cpu} \
        --vector_percent 1 \
        ~{skesa_opts}
    fi

    if [ -f "~{samplename}_skesa_contigs.fasta" ]; then
      echo "PASS" > SKESA_STATUS
      if [ ! -s "~{samplename}_skesa_contigs.fasta" ]; then
        echo "ERROR: No contigs were assembled"
        echo "FAIL" > SKESA_STATUS
      fi
    else
      echo "ERROR: Skesa assembly failed"
      echo "FAIL" > SKESA_STATUS
    fi
  >>>
  output {
    File? assembly_fasta = "~{samplename}_skesa_contigs.fasta"
    String skesa_version = read_string("VERSION")
    String skesa_docker = "~{docker}"
    String skesa_status = read_string("SKESA_STATUS")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}