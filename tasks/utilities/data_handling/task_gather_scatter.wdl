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
    import os

    def load_json_data(file_path):
      if os.path.exists(file_path):
        with open(file_path, 'r') as file:
          json_data_from_file = json.load(file)
          df_from_file = pd.DataFrame(json_data_from_file, columns=[os.basename(file_path)])
          return df_from_file
      else:
        return None

    df = pd.DataFrame()
    
    taxon_ids_df = load_json_data("~{taxon_ids}")
    if taxon_ids_df is not None:
      df = pd.concat([df, taxon_ids_df], axis=1)
    
    print(df)


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