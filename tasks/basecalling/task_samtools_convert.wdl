version 1.0

task samtools_convert {
  input {
    Array[File] sam_files
    Int cpu = 4
    Int memory = 16
    Int disk_size = 50
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15"
  }

  command <<< 
    set -euo pipefail

    # samtools version
    samtools --version | tee SAMTOOLS_VERSION

    bam_output="output/bam/"
    mkdir -p "$bam_output"

    for sam_file in ~{sep=" " sam_files}; do
      base_name=$(basename "$sam_file" .sam)
      bam_file="${bam_output}/${base_name}.bam"

      echo "Converting SAM to BAM: $sam_file -> $bam_file"
      samtools view -@ ~{cpu} -bS "$sam_file" > "$bam_file" || { echo "ERROR: samtools failed to convert $sam_file to BAM" >&2; exit 1; }

      # Count reads in the BAM file
      read_count=$(samtools view -c "$bam_file")
      
      # Log BAM file name and read count to stdout
      echo "BAM file: $bam_file, Read count: $read_count"
    done

    echo "### Listing BAM files in $bam_output ###"
    ls -lh "$bam_output"
  >>>

  output {
    Array[File] bam_files = glob("output/bam/*.bam")
    String samtools_version = read_string("SAMTOOLS_VERSION")
    String samtools_docker = docker
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 0
    maxRetries: 1
  }
}
