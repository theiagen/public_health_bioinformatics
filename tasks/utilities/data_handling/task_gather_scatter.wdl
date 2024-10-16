version 1.0

task gather_scatter {
  input {
    String samplename
    Array[Int]? taxon_ids
    # krakentools outputs
    Array[String?]+ organism
    Array[File?]+ extracted_read1
    Array[File?]+ extracted_read2
    Array[String]? krakentools_docker
    # fastq_scan outputs
    Array[Int?]+ fastq_scan_num_reads_binned1
    Array[Int?]+ fastq_scan_num_reads_binned2
    Array[String?]+ fastq_scan_num_reads_binned_pairs
    Array[String?]+ fastq_scan_docker
    Array[String?]+ fastq_scan_version
    # Assembly 
    Array[File?]+ pilon_assembly_fasta ### maybe?????
    # quast outputs
    Array[Int?]+ quast_genome_length
    Array[Int?]+ quast_number_contigs
    Array[Int?]+ quast_n50
    Array[Float?]+ quast_gc_percent
    # consensus qc outputs
    Array[Int?]+ number_N
    Array[Int?]+ number_ATCG
    Array[Int?]+ number_Degenerate
    Array[Int?]+ number_Total
    Array[Float?]+ percent_reference_coverage
    # pangolin outputs
    Array[String?]+ pango_lineage
    Array[String?]+ pango_lineage_expanded
    Array[String?]+ pangolin_conflicts
    Array[String?]+ pangolin_notes
    Array[String?]+ pangolin_assignment_version
    Array[String?]+ pangolin_versions
    Array[String?]+ pangolin_docker
    # Nextclade outputs for non-flu
    Array[String?]+ nextclade_version
    Array[String?]+ nextclade_docker
    Array[String?]+ nextclade_ds_tag
    Array[String?]+ nextclade_aa_subs
    Array[String?]+ nextclade_aa_dels
    Array[String?]+ nextclade_clade
    Array[String?]+ nextclade_lineage
    Array[String?]+ nextclade_qc
    # Nextclade outputs for flu HA
    Array[String?]+ nextclade_ds_tag_flu_ha
    Array[String?]+ nextclade_aa_subs_flu_ha
    Array[String?]+ nextclade_aa_dels_flu_ha
    Array[String?]+ nextclade_clade_flu_ha
    Array[String?]+ nextclade_qc_flu_ha
    # Nextclade outputs for flu NA
    Array[String?]+ nextclade_ds_tag_flu_na
    Array[String?]+ nextclade_aa_subs_flu_na
    Array[String?]+ nextclade_aa_dels_flu_na
    Array[String?]+ nextclade_clade_flu_na
    Array[String?]+ nextclade_qc_flu_na
    String docker = "us-docker.pkg.dev/general-theiagen/quay/ubuntu:latest"
    Int disk_size = 50
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    (
      echo -e "taxon_ids\torganism\textracted_read1\textracted_read2\tkrakentools_docker\tfastq_scan_num_reads_binned1\tfastq_scan_num_reads_binned2\tfastq_scan_num_reads_binned_pairs\tfastq_scan_docker\tfastq_scan_version\tpilon_assembly_fasta\tquast_genome_length\tquast_number_contigs\tquast_n50\tquast_gc_percent\tnumber_N\tnumber_ATCG\tnumber_Degenerate\tnumber_Total\tpercent_reference_coverage\tpango_lineage\tpango_lineage_expanded\tpangolin_conflicts\tpangolin_notes\tpangolin_assignment_version\tpangolin_versions\tpangolin_docker\tnextclade_version\tnextclade_docker\tnextclade_ds_tag\tnextclade_aa_subs\tnextclade_aa_dels\tnextclade_clade\tnextclade_lineage\tnextclade_qc\tnextclade_ds_tag_flu_ha\tnextclade_aa_subs_flu_ha\tnextclade_aa_dels_flu_ha\tnextclade_clade_flu_ha\tnextclade_qc_flu_ha\tnextclade_ds_tag_flu_na\tnextclade_aa_subs_flu_na\tnextclade_aa_dels_flu_na\tnextclade_clade_flu_na\tnextclade_qc_flu_na"
      paste <(echo "~{sep="\n" taxon_ids}") \
        <(echo "~{sep="\n" organism}") \
        <(echo "~{sep="\n" extracted_read1}") \
        <(echo "~{sep="\n" extracted_read2}") \
        <(echo "~{sep="\n" krakentools_docker}") \
        <(echo "~{sep="\n" fastq_scan_num_reads_binned1}") \
        <(echo "~{sep="\n" fastq_scan_num_reads_binned2}") \
        <(echo "~{sep="\n" fastq_scan_num_reads_binned_pairs}") \
        <(echo "~{sep="\n" fastq_scan_docker}") \
        <(echo "~{sep="\n" fastq_scan_version}") \
        <(echo "~{sep="\n" pilon_assembly_fasta}") \
        <(echo "~{sep="\n" quast_genome_length}") \
        <(echo "~{sep="\n" quast_number_contigs}") \
        <(echo "~{sep="\n" quast_n50}") \
        <(echo "~{sep="\n" quast_gc_percent}") \
        <(echo "~{sep="\n" number_N}") \
        <(echo "~{sep="\n" number_ATCG}") \
        <(echo "~{sep="\n" number_Degenerate}") \
        <(echo "~{sep="\n" number_Total}") \
        <(echo "~{sep="\n" percent_reference_coverage}") \
        <(echo "~{sep="\n" pango_lineage}") \
        <(echo "~{sep="\n" pango_lineage_expanded}") \
        <(echo "~{sep="\n" pangolin_conflicts}") \
        <(echo "~{sep="\n" pangolin_notes}") \
        <(echo "~{sep="\n" pangolin_assignment_version}") \
        <(echo "~{sep="\n" pangolin_versions}") \
        <(echo "~{sep="\n" pangolin_docker}") \
        <(echo "~{sep="\n" nextclade_version}") \
        <(echo "~{sep="\n" nextclade_docker}") \
        <(echo "~{sep="\n" nextclade_ds_tag}") \
        <(echo "~{sep="\n" nextclade_aa_subs}") \
        <(echo "~{sep="\n" nextclade_aa_dels}") \
        <(echo "~{sep="\n" nextclade_clade}") \
        <(echo "~{sep="\n" nextclade_lineage}") \
        <(echo "~{sep="\n" nextclade_qc}") \
        <(echo "~{sep="\n" nextclade_ds_tag_flu_ha}") \
        <(echo "~{sep="\n" nextclade_aa_subs_flu_ha}") \
        <(echo "~{sep="\n" nextclade_aa_dels_flu_ha}") \
        <(echo "~{sep="\n" nextclade_clade_flu_ha}") \
        <(echo "~{sep="\n" nextclade_qc_flu_ha}") \
        <(echo "~{sep="\n" nextclade_ds_tag_flu_na}") \
        <(echo "~{sep="\n" nextclade_aa_subs_flu_na}") \
        <(echo "~{sep="\n" nextclade_aa_dels_flu_na}") \
        <(echo "~{sep="\n" nextclade_clade_flu_na}") \
        <(echo "~{sep="\n" nextclade_qc_flu_na}")
    ) > ~{samplename}.results.tsv
  >>>
  output {
    File gathered_results = "~{samplename}.results.tsv"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}