version 1.0

task samtools_process {
  input {
    File sam_file            # Input SAM file
    String sample_name    
    Int cpu = 4             
    Int memory = 16          
    Int disk_size = 50      
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15" 
  }
  command <<< 
    set -euo pipefail

    echo "Samtools version:" 
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # Define output file names based on sample name
    bam_file="~{sample_name}.bam"
    sorted_bam_file="~{sample_name}.sorted.bam"
    coverage_file="~{sample_name}.coverage"

    # Step 1: Convert SAM to BAM
    echo "Converting SAM to BAM: ~{sam_file} -> $bam_file"
    samtools view -@ ~{cpu} -bS "~{sam_file}" > "$bam_file"

    # Step 2: Sort BAM
    echo "Sorting BAM file: $bam_file -> $sorted_bam_file"
    samtools sort -@ ~{cpu} "$bam_file" -o "$sorted_bam_file"

    # Step 3: Index BAM
    echo "Indexing BAM file: $sorted_bam_file"
    samtools index "$sorted_bam_file"

    # Step 4: Calculate coverage
    echo "Calculating coverage: $sorted_bam_file -> $coverage_file"
    samtools depth "$sorted_bam_file" > "$coverage_file"

    echo "### BAM and Coverage files created successfully ###"
  >>>
  output {
    File sorted_bam_file = "~{sample_name}.sorted.bam"         # Sorted BAM file
    File bam_index = "~{sample_name}.sorted.bam.bai"           # BAM index file
    File coverage_file = "~{sample_name}.coverage"             # Coverage file
    String samtools_version = read_string("SAMTOOLS_VERSION") 
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk ~{disk_size} SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
