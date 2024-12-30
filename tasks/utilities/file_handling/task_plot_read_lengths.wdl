version 1.0

task plot_read_length_distribution {
  input {
    File filtered_fastq
    String sample_name
    Int cpu = 1
    Int memory = 4
    Int disk_size = 10
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.2"
  }
  command <<<
    set -eu pipefail
    # Ensure matplotlib is installed
    python3 -m pip install --no-cache-dir matplotlib

    # Python script to plot read length distribution
    python3 <<EOF
import matplotlib.pyplot as plt

sample_name = "~{sample_name}"
fastq_file = "~{filtered_fastq}"
output_plot = f"{sample_name}_read_length_distribution.png"

read_lengths = []
with open(fastq_file, "rt") as f:
    for i, line in enumerate(f):
        if i % 4 == 1:  # Sequence lines in FASTQ
            read_lengths.append(len(line.strip()))

plt.figure(figsize=(10, 6))
plt.hist(read_lengths, bins=range(min(read_lengths), max(read_lengths) + 1, 10), edgecolor="black")
plt.title(f"Read Length Distribution for {sample_name}")
plt.xlabel("Read Length (bp)")
plt.ylabel("Frequency")
plt.savefig(output_plot)
plt.close()
EOF
  >>>
  output {
    File read_length_plot = "~{sample_name}_read_length_distribution.png"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}