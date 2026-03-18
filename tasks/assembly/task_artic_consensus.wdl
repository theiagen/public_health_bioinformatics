version 1.0

task consensus {
  input {
    String samplename
    String? organism
    File read1
    File? primer_bed
    File? reference_genome
    Int normalise = 20000
    Int cpu = 8
    Int memory = 16
    Int disk_size = 100
    String clair3_model = "r1041_e82_400bps_sup_v500"
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/artic:1.9.0"
  }
  command <<<
    set -euo pipefail

    # Newer version of Artic use a different command to launch the pipeline
    # when using user provided files for bed and reference genome are provided, 
    # we can now provide the bed and reference genome files to the pipeline directly 
    if [[ -z "~{primer_bed}" || -z "~{reference_genome}" ]]; then
      echo "Error: Both primer bed and reference genome must be provided" >&2
      exit 1
    fi

    # Copy reference files to working directory, this is necessary for fiadx to work
    # which is a requirement of clair3, we run into similar issues in the clair3_variants task
    cp "~{reference_genome}" reference.fasta
    # Run ARTIC with user provided files
    artic minion --model ~{clair3_model} \
      --normalise ~{normalise} \
      --threads ~{cpu} \
      --bed ~{primer_bed} \
      --ref reference.fasta \
      --read-file ~{read1} \
      ~{samplename}

    # Set output file contents for local scheme
    head -n1 "~{reference_genome}" | sed 's/>//' > REFERENCE_GENOME
    basename "~{primer_bed}" > PRIMER_NAME
  fi

  # Capture ARTIC version
  artic -v > VERSION

  # Grab reads from alignment - cdph wants this
  # 0x904 means we are now filtering out unaligned, secondary, and supplemental alignments - thanks Curtis
  samtools fastq -F0x904 ~{samplename}.primertrimmed.rg.sorted.bam | gzip > ~{samplename}.fastq.gz  
  >>>
  output {
    File consensus_seq = "~{samplename}.consensus.fasta"
    File artic_clair3_pass_vcf = "~{samplename}.pass.vcf"
    String artic_clair3_model = clair3_model
    String artic_pipeline_reference = read_string("REFERENCE_GENOME")
    String artic_pipeline_version = read_string("VERSION")
    String artic_pipeline_docker = docker
    File sorted_bam = "~{samplename}.trimmed.rg.sorted.bam"
    File trim_sorted_bam = "~{samplename}.primertrimmed.rg.sorted.bam"
    File trim_sorted_bai = "~{samplename}.primertrimmed.rg.sorted.bam.bai"
    File? reads_aligned = "~{samplename}.fastq.gz"
    String primer_bed_name = read_string("PRIMER_NAME")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}