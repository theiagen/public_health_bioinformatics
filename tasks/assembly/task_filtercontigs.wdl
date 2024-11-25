version 1.0

task contig_filter {
  input {
    File assembly_fasta
    Int min_len
    Float min_cov = 0.0
    Boolean filter_homopolymers = true
    Int cpu = 4
    Int memory = 8
    Int disk_size = 50
    String docker = "biotools"
  }
  command <<< 
    set -e
    echo "Starting contig filtering"

    # Define output file paths
    output_filtered="filtered_contigs.fasta"

    # Perform filtering using seqkit and awk
    seqkit fx2tab -l -i -H ~{assembly_fasta} | \
    awk -v min_len=~{min_len} -v min_cov=~{min_cov} -v filter_homopolymers=~{filter_homopolymers} '
    BEGIN { OFS="\t" }
    {
        # Parse length and optional coverage from headers
        length = $2
        coverage = (match($1, /cov=([0-9.]+)/, arr) ? arr[1] : 0)

        # Filter criteria
        if (length >= min_len && coverage >= min_cov) {
            # Remove homopolymers if enabled
            if (filter_homopolymers) {
                if ($3 !~ /^(A+|T+|C+|G+)$/) {
                    print $1, $2, $3
                }
            } else {
                print $1, $2, $3
            }
        }
    }' | seqkit tab2fx -o $output_filtered

    # Ensure the output is not empty
    if [ ! -s $output_filtered ]; then
        echo "Error: No contigs passed filtering criteria!" >&2
        exit 1
    fi

    echo "Filtering complete. Filtered contigs written to $output_filtered."
  >>>
  output {
    File filtered_fasta = "filtered_contigs.fasta"
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
