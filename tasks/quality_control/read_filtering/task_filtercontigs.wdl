version 1.0

task contig_filter {
  input {
    File assembly_fasta
    Int min_len = 1000
    Boolean filter_homopolymers = true
    Int disk_size = 100
    Int memory = 16
    Int cpu = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/seqkit:2.8.2"
  }
  command <<< 
    set -euo pipefail
    echo "Filtering contigs from ~{assembly_fasta}" >&2

    # Calculate initial metrics
    total_contigs=$(grep -c "^>" ~{assembly_fasta})
    echo "Total contigs: $total_contigs" >&2

    # Calculate total base pairs in all contigs and output individual lengths
    awk '/^>/ {if (seqlen) {print seqlen; total+=seqlen}; seqlen=0; next} {seqlen+=length($0)} END {if (seqlen) {print seqlen; total+=seqlen}; print "Total:" total}' ~{assembly_fasta} > contig_lengths_before.txt
    total_bases=$(tail -n 1 contig_lengths_before.txt | cut -d ':' -f 2)
    echo "Total bases: $total_bases" >&2

    # Filter contigs by length
    seqkit seq -m ~{min_len} ~{assembly_fasta} > length_filtered.fasta
    echo "Length filtering complete. File size:" >&2
    ls -lh length_filtered.fasta >&2

    # Optionally filter out contigs with homopolymers
    if [ "~{filter_homopolymers}" = "true" ]; then
      awk '
      BEGIN {RS=">"; ORS=""} 
      NR > 1 {
        split($0, lines, "\n");
        header = lines[1];
        seq = "";
        for (i=2; i<=length(lines); i++) {
          seq = seq lines[i];
        }
        if (seq !~ /^(.)\1+$/) {
          print ">" header "\n" seq "\n"
        }
      }' length_filtered.fasta > filtered_contigs.fasta
    else
      mv length_filtered.fasta filtered_contigs.fasta
    fi

    # Validate the final file is not empty
    if [ ! -s filtered_contigs.fasta ]; then
      echo "Error: No contigs passed filtering criteria!" >&2
      exit 1
    fi

    # Calculate final metrics for the filtered contigs
    awk '/^>/ {if (seqlen) {print seqlen; total+=seqlen}; seqlen=0; next} {seqlen+=length($0)} END {if (seqlen) {print seqlen; total+=seqlen}; print "Total:" total}' filtered_contigs.fasta > contig_lengths_after.txt
    retained_contigs=$(grep -c "^>" filtered_contigs.fasta)
    retained_bases=$(tail -n 1 contig_lengths_after.txt | cut -d ':' -f 2)

    # Calculate removed contigs
    contigs_removed_short=$((total_contigs - $(grep -c "^>" length_filtered.fasta)))
    contigs_removed_homopolymers=$((total_contigs - retained_contigs - contigs_removed_short))

    # Write filtering metrics to a summary file
    metrics_file="filtering_metrics.txt"
    {
      echo "Total contigs: $total_contigs"
      echo "Total bases: $total_bases"
      echo "Contigs retained: $retained_contigs"
      echo "Bases retained: $retained_bases"
      echo "Contigs removed (short length): $contigs_removed_short"
      echo "Contigs removed (homopolymers): $contigs_removed_homopolymers"
    } > $metrics_file

    # Output individual contig lengths before and after filtering
    echo "Contig lengths before filtering:" >> $metrics_file
    cat contig_lengths_before.txt >> $metrics_file
    echo "Contig lengths after filtering:" >> $metrics_file
    cat contig_lengths_after.txt >> $metrics_file

    cat $metrics_file >&2
  >>>
  output {
    File filtered_fasta = "filtered_contigs.fasta"
    File assembly_filtering_metrics = "filtering_metrics.txt"
    File contig_lengths_before = "contig_lengths_before.txt"
    File contig_lengths_after = "contig_lengths_after.txt"
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory}G"
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 3
    preemptible: 0
  }
}
