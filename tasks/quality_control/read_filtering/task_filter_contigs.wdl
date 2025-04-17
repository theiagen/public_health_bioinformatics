version 1.0

task filter_contigs {
  input {
    File assembly_fasta
    Int min_length = 1000
    Int disk_size = 50
    Int memory = 8
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/kaptive:2.0.3" # Staph-b docker image with Biopython
  }
  command <<< 
    set -euo pipefail
    
    echo "Filtering contigs from ~{assembly_fasta}" >&2


    # Run Biopython script for filtering
    python3 <<EOF
import sys
from Bio import SeqIO

# Input parameters
assembly_fasta = "~{assembly_fasta}"
min_length = ~{min_length}

# Read sequences
records = list(SeqIO.parse(assembly_fasta, "fasta"))
total_contigs = len(records)
total_bases = sum(len(record) for record in records)

# Filter by length
length_filtered = [record for record in records if len(record) >= min_length]

# Filter out homopolymers
def has_homopolymers(seq):
    return any(base * len(seq) == seq for base in "ACGT")

filtered_records = [record for record in length_filtered if not has_homopolymers(str(record.seq))]

# Output results
SeqIO.write(filtered_records, "filtered_contigs.fasta", "fasta")

# Metrics output
retained_contigs = len(filtered_records)
retained_bases = sum(len(record) for record in filtered_records)

with open("filtering_metrics.txt", "w") as metrics_file:
    metrics_file.write("Contig Filtering Metrics\\n")
    metrics_file.write("========================\\n")
    metrics_file.write(f"Total contigs before filtering: {total_contigs}\\n")
    metrics_file.write(f"Total sequence length before filtering: {total_bases} bases\\n")
    metrics_file.write(f"Total contigs after filtering: {retained_contigs}\\n")
    metrics_file.write(f"Total sequence length after filtering: {retained_bases} bases\\n")
    metrics_file.write(f"Contigs removed (short length): {total_contigs - len(length_filtered)}\\n")
    metrics_file.write(f"Contigs removed (homopolymers): {len(length_filtered) - retained_contigs}\\n")
    metrics_file.write("========================\\n")
    metrics_file.write("Filtering completed successfully.\\n")
EOF
  >>>
  output {
    File filtered_fasta = "filtered_contigs.fasta"
    File assembly_filtering_metrics = "filtering_metrics.txt"
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 3
    preemptible: 0
  }
}
