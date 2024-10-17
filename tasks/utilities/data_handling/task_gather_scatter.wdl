version 1.0

task gather_scatter {
  input {
    String samplename
    File? taxon_ids
    # krakentools outputs
    File? organism
    File? extracted_read1
    File? extracted_read2
    File? krakentools_docker
    # fastq_scan outputs
    File? fastq_scan_num_reads_binned1
    File? fastq_scan_num_reads_binned2
    File? fastq_scan_num_reads_binned_pairs
    File? fastq_scan_docker
    File? fastq_scan_version
    # Assembly 
    File? pilon_assembly_fasta### maybe?????
    # quast outputs
    File? quast_genome_length
    File? quast_number_contigs
    File? quast_n50
    File? quast_gc_percent
    # consensus qc outputs
    File? number_N
    File? number_ATCG
    File? number_Degenerate
    File? number_Total
    File? percent_reference_coverage
    # pangolin outputs
    File? pango_lineage
    File? pango_lineage_expanded
    File? pangolin_conflicts
    File? pangolin_notes
    File? pangolin_assignment_version
    File? pangolin_versions
    File? pangolin_docker
    # Nextclade outputs for non-flu
    File? nextclade_version
    File? nextclade_docker
    File? nextclade_ds_tag
    File? nextclade_aa_subs
    File? nextclade_aa_dels
    File? nextclade_clade
    File? nextclade_lineage
    File? nextclade_qc
    # Nextclade outputs for flu HA
    File? nextclade_ds_tag_flu_ha
    File? nextclade_aa_subs_flu_ha
    File? nextclade_aa_dels_flu_ha
    File? nextclade_clade_flu_ha
    File? nextclade_qc_flu_ha
    # Nextclade outputs for flu NA
    File? nextclade_ds_tag_flu_na
    File? nextclade_aa_subs_flu_na
    File? nextclade_aa_dels_flu_na
    File? nextclade_clade_flu_na
    File? nextclade_qc_flu_na
    # change to be a docker with pandas
    String docker = "us-docker.pkg.dev/general-theiagen/quay/ubuntu:latest"
    Int disk_size = 50
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    python3<<CODE
    import pandas as pd
    
    taxon_ids = pd.read_csv("~{taxon_ids}", header=None, names=["taxon_id"])
    organism = pd.read_csv("~{organism}", header=None, names=["organism"])
    extracted_read1 = pd.read_csv("~{extracted_read1}", header=None, names=["extracted_read1"])
    extracted_read2 = pd.read_csv("~{extracted_read2}", header=None, names=["extracted_read2"])
    krakentools_docker = pd.read_csv("~{krakentools_docker}", header=None, names=["krakentools_docker"])
    fastq_scan_num_reads_binned1 = pd.read_csv("~{fastq_scan_num_reads_binned1}", header=None, names=["fastq_scan_num_reads_binned1"])
    fastq_scan_num_reads_binned2 = pd.read_csv("~{fastq_scan_num_reads_binned2}", header=None, names=["fastq_scan_num_reads_binned2"])
    fastq_scan_num_reads_binned_pairs = pd.read_csv("~{fastq_scan_num_reads_binned_pairs}", header=None, names=["fastq_scan_num_reads_binned_pairs"])
    fastq_scan_docker = pd.read_csv("~{fastq_scan_docker}", header=None, names=["fastq_scan_docker"])
    fastq_scan_version = pd.read_csv("~{fastq_scan_version}", header=None, names=["fastq_scan_version"])
    pilon_assembly_fasta = pd.read_csv("~{pilon_assembly_fasta}", header=None, names=["pilon_assembly_fasta"])
    quast_genome_length = pd.read_csv("~{quast_genome_length}", header=None, names=["quast_genome_length"])
    quast_number_contigs = pd.read_csv("~{quast_number_contigs}", header=None, names=["quast_number_contigs"])
    quast_n50 = pd.read_csv("~{quast_n50}", header=None, names=["quast_n50"])
    quast_gc_percent = pd.read_csv("~{quast_gc_percent}", header=None, names=["quast_gc_percent"])
    number_N = pd.read_csv("~{number_N}", header=None, names=["number_N"])
    number_ATCG = pd.read_csv("~{number_ATCG}", header=None, names=["number_ATCG"])
    number_Degenerate = pd.read_csv("~{number_Degenerate}", header=None, names=["number_Degenerate"])
    number_Total = pd.read_csv("~{number_Total}", header=None, names=["number_Total"])
    percent_reference_coverage = pd.read_csv("~{percent_reference_coverage}", header=None, names=["percent_reference_coverage"])
    pango_lineage = pd.read_csv("~{pango_lineage}", header=None, names=["pango_lineage"])
    pango_lineage_expanded = pd.read_csv("~{pango_lineage_expanded}", header=None, names=["pango_lineage_expanded"])
    pangolin_conflicts = pd.read_csv("~{pangolin_conflicts}", header=None, names=["pangolin_conflicts"])
    pangolin_notes = pd.read_csv("~{pangolin_notes}", header=None, names=["pangolin_notes"])
    pangolin_assignment_version = pd.read_csv("~{pangolin_assignment_version}", header=None, names=["pangolin_assignment_version"])
    pangolin_versions = pd.read_csv("~{pangolin_versions}", header=None, names=["pangolin_versions"])
    pangolin_docker = pd.read_csv("~{pangolin_docker}", header=None, names=["pangolin_docker"])
    nextclade_version = pd.read_csv("~{nextclade_version}", header=None, names=["nextclade_version"])
    nextclade_docker = pd.read_csv("~{nextclade_docker}", header=None, names=["nextclade_docker"])
    nextclade_ds_tag = pd.read_csv("~{nextclade_ds_tag}", header=None, names=["nextclade_ds_tag"])
    nextclade_aa_subs = pd.read_csv("~{nextclade_aa_subs}", header=None, names=["nextclade_aa_subs"])
    nextclade_aa_dels = pd.read_csv("~{nextclade_aa_dels}", header=None, names=["nextclade_aa_dels"])
    nextclade_clade = pd.read_csv("~{nextclade_clade}", header=None, names=["nextclade_clade"])
    nextclade_lineage = pd.read_csv("~{nextclade_lineage}", header=None, names=["nextclade_lineage"])
    nextclade_qc = pd.read_csv("~{nextclade_qc}", header=None, names=["nextclade_qc"])
    nextclade_ds_tag_flu_ha = pd.read_csv("~{nextclade_ds_tag_flu_ha}", header=None, names=["nextclade_ds_tag_flu_ha"])
    nextclade_aa_subs_flu_ha = pd.read_csv("~{nextclade_aa_subs_flu_ha}", header=None, names=["nextclade_aa_subs_flu_ha"])
    nextclade_aa_dels_flu_ha = pd.read_csv("~{nextclade_aa_dels_flu_ha}", header=None, names=["nextclade_aa_dels_flu_ha"])
    nextclade_clade_flu_ha = pd.read_csv("~{nextclade_clade_flu_ha}", header=None, names=["nextclade_clade_flu_ha"])
    nextclade_qc_flu_ha = pd.read_csv("~{nextclade_qc_flu_ha}", header=None, names=["nextclade_qc_flu_ha"])
    nextclade_ds_tag_flu_na = pd.read_csv("~{nextclade_ds_tag_flu_na}", header=None, names=["nextclade_ds_tag_flu_na"])
    nextclade_aa_subs_flu_na = pd.read_csv("~{nextclade_aa_subs_flu_na}", header=None, names=["nextclade_aa_subs_flu_na"])
    nextclade_aa_dels_flu_na = pd.read_csv("~{nextclade_aa_dels_flu_na}", header=None, names=["nextclade_aa_dels_flu_na"])
    nextclade_clade_flu_na = pd.read_csv("~{nextclade_clade_flu_na}", header=None, names=["nextclade_clade_flu_na"])
    nextclade_qc_flu_na = pd.read_csv("~{nextclade_qc_flu_na}", header=None, names=["nextclade_qc_flu_na"])
    






    CODE



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