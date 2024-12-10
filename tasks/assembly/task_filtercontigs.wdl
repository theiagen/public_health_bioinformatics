version 1.0

task contig_filter {
  input {
    File assembly_fasta
    Int min_len = 1000
    Boolean filter_homopolymers = true
    Int disk_size = 100
    Int memory = 16
    Int threads = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/seqkit:2.8.2"
  }
  command <<< 
    set -euo pipefail
    echo "Filtering contigs from ~{assembly_fasta}" >&2

    # Calculate initial metrics
    total_contigs=$(grep -c "^>" ~{assembly_fasta})
    echo "Total contigs: $total_contigs" >&2

    # Calculate total base pairs in all contigs
    total_bases=$(awk '/^>/ {if (seqlen){total+=seqlen}; seqlen=0; next} {seqlen += length($0)} END {total+=seqlen; print total}' ~{assembly_fasta})
    echo "Total bases: $total_bases" >&2

    # Filter contigs by length
    # Retain only contigs that are >= minimum length and write to a new FASTA file
    seqkit seq -m ~{min_len} ~{assembly_fasta} > length_filtered.fasta
    echo "Length filtering complete. File size:" >&2
    ls -lh length_filtered.fasta >&2

    # Optionally filter out contigs with homopolymers
    if [ "~{filter_homopolymers}" = "true" ]; then
      awk '
      BEGIN {RS=">"; ORS=""} 
      NR > 1 {
        # Split each contig entry into header and sequence
        split($0, lines, "\n");
        header = lines[1];
        seq = "";
        for (i=2; i<=length(lines); i++) {
          seq = seq lines[i];
        }
         # Exclude sequences composed entirely of a single repeated base
        if (seq !~ /^(.)\1+$/) {
          print ">" header "\n" seq "\n"
        }
      }' length_filtered.fasta > filtered_contigs.fasta
    else
      # If homopolymer filtering is not enabled, use length-filtered file as final output
      mv length_filtered.fasta filtered_contigs.fasta
    fi

    # Validate the final file is not empty
    if [ ! -s filtered_contigs.fasta ]; then
      echo "Error: No contigs passed filtering criteria!" >&2
      exit 1
    fi

    # Calculate final metrics for the filtered contigs
    retained_contigs=$(grep -c "^>" filtered_contigs.fasta)
    retained_bases=$(awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length($0)} END {if (seqlen) print seqlen}' filtered_contigs.fasta)
    contigs_removed_homopolymers=$((total_contigs - retained_contigs))

    # Write filtering metrics to a summary file
    metrics_file="filtering_metrics.txt"
    echo "Total contigs: $total_contigs" > $metrics_file
    echo "Total bases: $total_bases" >> $metrics_file
    echo "Contigs retained: $retained_contigs" >> $metrics_file
    echo "Bases retained: $retained_bases" >> $metrics_file
    echo "Contigs removed (short length): $((total_contigs - retained_contigs))" >> $metrics_file
    echo "Contigs removed (homopolymers): $contigs_removed_homopolymers" >> $metrics_file

    cat $metrics_file >&2
  >>>
  output {
    File filtered_fasta = "filtered_contigs.fasta"
    File assembly_filtering_metrics = "filtering_metrics.txt"
  }
  runtime {
    docker: "~{docker}"
    cpu: threads
    memory: "~{memory}G"
    disks: "local-disk " + disk_size + " SSD"
  }
}
