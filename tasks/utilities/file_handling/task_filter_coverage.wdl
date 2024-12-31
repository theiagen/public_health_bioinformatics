version 1.0

task filter_coverage {
  input {
    File bam_file
    File reference
    Int min_coverage = 10
    Int memory = 8
    Int cpu = 4
    Int disk_size = 50
    String sample_name
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bedtools:2.31.1"
  }
  command <<< 
    set -euo pipefail

    # Generate coverage file from BAM file
    echo "Calculating coverage from BAM file: ~{bam_file}"
    bedtools genomecov -ibam ~{bam_file} -d > ~{sample_name}_coverage.txt

    # Step 2: Filter regions based on minimum coverage
    echo "Filtering regions with minimum coverage of ~{min_coverage}"
    # AWK explanation:
    #  - `total_positions`: Total positions processed in the reference.
    #  - `filtered_positions`: Count of positions meeting the minimum coverage threshold.
    #  - Regions are identified and grouped into continuous stretches based on:
    #    * The same reference sequence (`prev_ref`).
    #    * Consecutive positions (`prev_pos`).
    awk -v min_cov=~{min_coverage} '
        BEGIN { OFS="\t"; total_positions=0; filtered_positions=0 }
        {
            total_positions++
            if ($3 >= min_cov) { # Position meets the coverage threshold
                filtered_positions++
                if (prev_ref != $1 || $2 != prev_pos + 1) { # Start a new region
                    if (region_start) { # Print the previous region
                        print prev_ref, region_start, prev_pos
                    }
                    region_start = $2 # Set the new region start
                }
            } else if (region_start) { # Position below threshold; close the region
                print prev_ref, region_start, prev_pos
                region_start = "" # Reset region
            }
            prev_ref = $1 # Track reference name
            prev_pos = $2 # Track position
        }
        END {
            if (region_start) { # Print the last region if it exists
                print prev_ref, region_start, prev_pos
            }
            # Print statistics to STDERR for logging
            print "Total positions processed:", total_positions > "/dev/stderr"
            print "Filtered positions (coverage >= " min_cov "):", filtered_positions > "/dev/stderr"
        }
    ' ~{sample_name}_coverage.txt > ~{sample_name}_filtered_regions.txt

    # Step 3: Create BED file
    echo "Creating BED file"
    # Converts filtered regions into BED format (0-based start, 1-based end)
    awk '{print $1 "\t" $2-1 "\t" $3}' ~{sample_name}_filtered_regions.txt > ~{sample_name}_regions.bed

    # Step 4: Extract scaffolds using bedtools
    echo "Extracting scaffolds from reference genome"
    # Extracts sequences from the reference based on BED file
    bedtools getfasta -fi ~{reference} -bed ~{sample_name}_regions.bed -fo ~{sample_name}_scaffolds_temp.fasta
    # Remove extra information from scaffold headers (optional for cleaner output and downstream processing)
    sed 's/:.*//' ~{sample_name}_scaffolds_temp.fasta > ~{sample_name}_scaffolds.fasta

    echo "Filtering complete. Check stderr for statistics."
  >>>
  output {
    File scaffolds_fasta = "~{sample_name}_scaffolds.fasta"      # Extracted scaffolds
    File coverage_file = "~{sample_name}_coverage.txt"          # Raw coverage file
    File filtered_regions = "~{sample_name}_filtered_regions.txt" # Filtered regions
    File regions_bed = "~{sample_name}_regions.bed"             # BED file for filtered regions
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
