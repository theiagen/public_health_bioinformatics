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

  if [[ -z "~{primer_bed}" || -z "~{reference_genome}" ]]; then
    echo "Error: Both primer bed and reference genome must be provided" >&2
    exit 1
  fi

  # Copy reference files to working directory, this is necessary for fiadx to work
  # which is a requirement of clair3, we run into similar issues in the clair3_variants task
  cp "~{reference_genome}" reference.fasta

  # ARTIC claims to expect 6 column BED files, but BEDs w/o the sequence column fail primalbedtools
  # SPOOF the sequence column if a 6 column BED file is encountered 
  python3 <<CODE
  with open("~{primer_bed}", "r") as infile, open("primer.bed", "w") as outfile:
    for line in infile:
      parts = line.strip().split("\t")
      if len(parts) < 7:
        parts.append("SPOOF")
      outfile.write("\t".join(parts) + "\n")
  CODE

  # Run ARTIC with user-provided files
  artic minion --model ~{clair3_model} \
    --normalise ~{normalise} \
    --threads ~{cpu} \
    --bed primer.bed \
    --ref reference.fasta \
    --read-file ~{read1} \
    ~{samplename}

  # Set output file contents for local scheme
  head -n1 "~{reference_genome}" | sed 's/>//' > REFERENCE_GENOME
  basename "~{primer_bed}" > PRIMER_NAME

  # Capture ARTIC version
  artic -v > VERSION

  # Grab reads from alignment
  # 0x904 means we are now filtering out unaligned, secondary, and supplemental alignments - thanks Curtis
  samtools fastq -F0x904 ~{samplename}.primertrimmed.rg.sorted.bam | gzip > ~{samplename}.fastq.gz  

  # Calculate the primer trimming stats:
    # primer_trimmed_read_percent
    # per primer read counts (?)
  # Based on comparing *primertrimmed.rg.sorted.bam to *trimmed.rg.sorted.bam

  # calculate percent of reads that were primer trimmed
  primers_trimmed=$(tail -n+2 ~{samplename}.alignreport.tsv | wc -l)
  total_reads=$(samtools view -c ~{samplename}.sorted.bam)
  >>>
  output {
    File consensus_seq = "~{samplename}.consensus.fasta"
    File artic_clair3_pass_vcf = "~{samplename}.pass.vcf"
    File? artic_amplicon_depths = "~{samplename}.amplicon_depths.tsv"
    String artic_clair3_model = clair3_model
    String artic_pipeline_reference = read_string("REFERENCE_GENOME")
    String artic_pipeline_version = read_string("VERSION")
    String artic_pipeline_docker = docker
    File sorted_bam = "~{samplename}.sorted.bam"
    File sorted_bai = "~{samplename}.sorted.bam.bai"
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