version 1.0

task samtools_convert {
  input {
    Array[File] sam_files
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15"
  }

  command <<< 
    set -e

    bam_output="output/bam/"
    mkdir -p $bam_output

    for sam_file in ~{sep=" " sam_files}; do
      base_name=$(basename $sam_file .sam)
      bam_file="${bam_output}/${base_name}.bam"

      echo "Converting SAM to BAM: $sam_file -> $bam_file"
      samtools view -bS "$sam_file" > "$bam_file" || { echo "ERROR: samtools failed to convert $sam_file to BAM" >&2; exit 1; }

      echo "Conversion successful: $bam_file"
    done

    echo "### Listing BAM files in $bam_output ###"
    ls -lh $bam_output
  >>>

  output {
    Array[File] bam_files = glob("output/bam/*.bam")
  }

  runtime {
    docker: docker
    cpu: 4
    memory: "16GB"
    maxRetries: 0
  }
}
