version 1.0

task filter_coverage {
  input {
    File bam_file
    File reference
    Int min_coverage = 10
    Int memory = 8
    Int cpu = 4
    Int disk_size = 20
    String sample_name
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bedtools:2.31.1"
  }
  command <<< 
    set -euo pipefail

    # Generate coverage file from BAM file
    echo "Calculating coverage from BAM file: ~{bam_file}"
    bedtools genomecov -ibam ~{bam_file} -d > coverage.txt

    # Step 2: Filter regions based on minimum coverage
    echo "Filtering regions with minimum coverage of ~{min_coverage}"
    awk -v min_cov=~{min_coverage} '
        BEGIN { OFS="\t"; total_positions=0; filtered_positions=0 }
        {
            total_positions++
            if ($3 >= min_cov) {
                filtered_positions++
                if (prev_ref != $1 || $2 != prev_pos + 1) {
                    if (region_start) {
                        print prev_ref, region_start, prev_pos
                    }
                    region_start = $2
                }
            } else if (region_start) {
                print prev_ref, region_start, prev_pos
                region_start = ""
            }
            prev_ref = $1
            prev_pos = $2
        }
        END {
            if (region_start) {
                print prev_ref, region_start, prev_pos
            }
            print "Total positions processed:", total_positions > "/dev/stderr"
            print "Filtered positions (coverage >= " min_cov "):", filtered_positions > "/dev/stderr"
        }
    ' coverage.txt > filtered_regions.txt

    # Step 3: Create BED file
    echo "Creating BED file"
    awk '{print $1 "\t" $2-1 "\t" $3}' filtered_regions.txt > regions.bed

    # Step 4: Extract scaffolds using bedtools
    echo "Extracting scaffolds from reference genome"
    bedtools getfasta -fi ~{reference} -bed regions.bed -fo scaffolds_temp.fasta
    # Remove extra information from scaffold headers
    sed 's/:.*//' scaffolds_temp.fasta > scaffolds.fasta

    echo "Filtering complete. Check stderr for statistics."
  >>>
  output {
    File scaffolds_fasta = "scaffolds.fasta"      # Extracted scaffolds
    File coverage_file = "coverage.txt"          # Raw coverage file
    File filtered_regions = "filtered_regions.txt" # Filtered regions
    File regions_bed = "regions.bed"             # BED file for filtered regions
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
