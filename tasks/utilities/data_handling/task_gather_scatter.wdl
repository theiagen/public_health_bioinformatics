version 1.0

task gather_scatter {
  input {
    Array[Int]? taxon_ids
    Array[String?] organism
    Array[File?] extracted_read1
    Array[File?] extracted_read2
    Array[Int?] fastq_scan_num_reads_binned1
    Array[Int?] fastq_scan_num_reads_binned2
    Array[String?] fastq_scan_num_reads_binned_pairs
    Array[File?] pilon_assembly_fasta ### maybe?????
    Array[Int?] quast_genome_length
    Array[Int?] quast_number_contigs
    Array[Int?] quast_n50
    Array[Float?] quast_gc_percent
    Array[String?] pango_lineage
    Array[String?] pango_lineage_expanded
    Array[String?] pangolin_conflicts
    Array[String?] pangolin_notes
    Array[String?] pangolin_assignment_version
    Array[String?] pangolin_versions
    Array[String?] pangolin_docker
    Array[String?] nextclade_version
    Array[String?] nextclade_docker
    Array[String?] nextclade_ds_tag
    Array[String?] nextclade_aa_subs
    Array[String?] nextclade_aa_dels
    Array[String?] nextclade_clade
    Array[String?] nextclade_lineage
    Array[String?] nextclade_qc


  }
  command <<<
    echo "taxon_ids: ~{sep="," taxon_ids}"
    echo "organism: ~{sep="," organism}"
    echo "extracted_read1: ~{sep="," extracted_read1}"
    echo "extracted_read2: ~{sep="," extracted_read2}"
    echo "fastq_scan_num_reads_binned1: ~{sep="," fastq_scan_num_reads_binned1}"
    echo "fastq_scan_num_reads_binned2: ~{sep="," fastq_scan_num_reads_binned2}"
    echo "fastq_scan_num_reads_binned_pairs: ~{sep="," fastq_scan_num_reads_binned_pairs}"
    echo "pilon_assembly_fasta: ~{sep="," pilon_assembly_fasta}"
    echo "quast_genome_length: ~{sep="," quast_genome_length}"
    echo "quast_number_contigs: ~{sep="," quast_number_contigs}"
    echo "quast_n50: ~{sep="," quast_n50}"
    echo "quast_gc_percent: ~{sep="," quast_gc_percent}"
    echo "pango_lineage: ~{sep="," pango_lineage}"
    echo "pango_lineage_expanded: ~{sep="," pango_lineage_expanded}"
    echo "pangolin_conflicts: ~{sep="," pangolin_conflicts}"
    echo "pangolin_notes: ~{sep="," pangolin_notes}"
    echo "pangolin_assignment_version: ~{sep="," pangolin_assignment_version}"
    echo "pangolin_versions: ~{sep="," pangolin_versions}"
    echo "pangolin_docker: ~{sep="," pangolin_docker}"
    echo "nextclade_version: ~{sep="," nextclade_version}"
    echo "nextclade_docker: ~{sep="," nextclade_docker}"
    echo "nextclade_ds_tag: ~{sep="," nextclade_ds_tag}"
    echo "nextclade_aa_subs: ~{sep="," nextclade_aa_subs}"
    echo "nextclade_aa_dels: ~{sep="," nextclade_aa_dels}"
    echo "nextclade_clade: ~{sep="," nextclade_clade}"
    echo "nextclade_lineage: ~{sep="," nextclade_lineage}"
    echo "nextclade_qc: ~{sep="," nextclade_qc}"


    # turn into tsv?
    # output to file

  >>>
  output {
    Array[Int]? taxon_ids_out = taxon_ids
  }
}