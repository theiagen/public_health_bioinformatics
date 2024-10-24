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
    File? metaspades_warning
    File? pilon_warning
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
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    Int disk_size = 50
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    python3<<CODE
    import pandas as pd
    import numpy as np
    import os
    import json

    def load_json_data(file_path, column_name, df):
      if os.path.exists(file_path):
        with open(file_path, 'r') as file:
          json_data_from_file = json.load(file)
          df_from_file = pd.DataFrame(json_data_from_file, columns=[column_name])
          df = pd.concat([df, df_from_file], axis=1)
          return df
      else:
        return None

    df = pd.DataFrame()
    
    df = load_json_data("~{taxon_ids}", "taxon_ids", df)
    df = load_json_data("~{organism}", "organism", df)
    df = load_json_data("~{extracted_read1}", "extracted_read1", df)
    df = load_json_data("~{extracted_read2}", "extracted_read2", df)
    df = load_json_data("~{krakentools_docker}", "krakentools_docker", df)
    df = load_json_data("~{fastq_scan_num_reads_binned1}", "fastq_scan_num_reads_binned1", df)
    df = load_json_data("~{fastq_scan_num_reads_binned2}", "fastq_scan_num_reads_binned2", df)
    df = load_json_data("~{fastq_scan_num_reads_binned_pairs}", "fastq_scan_num_reads_binned_pairs", df)
    df = load_json_data("~{fastq_scan_docker}", "fastq_scan_docker", df)
    df = load_json_data("~{fastq_scan_version}", "fastq_scan_version", df)
    df = load_json_data("~{metaspades_warning}", "metaspades_warning", df)
    df = load_json_data("~{pilon_warning}", "pilon_warning", df)
    df = load_json_data("~{pilon_assembly_fasta}", "pilon_assembly_fasta", df)
    df = load_json_data("~{quast_genome_length}", "quast_genome_length", df)
    df = load_json_data("~{quast_number_contigs}", "quast_number_contigs", df)
    df = load_json_data("~{quast_n50}", "quast_n50", df)
    df = load_json_data("~{quast_gc_percent}", "quast_gc_percent", df)
    df = load_json_data("~{number_N}", "number_N", df)
    df = load_json_data("~{number_ATCG}", "number_ATCG", df)
    df = load_json_data("~{number_Degenerate}", "number_Degenerate", df)
    df = load_json_data("~{number_Total}", "number_Total", df)
    df = load_json_data("~{percent_reference_coverage}", "percent_reference_coverage", df)
    df = load_json_data("~{pango_lineage}", "pango_lineage", df)
    df = load_json_data("~{pango_lineage_expanded}", "pango_lineage_expanded", df)
    df = load_json_data("~{pangolin_conflicts}", "pangolin_conflicts", df)
    df = load_json_data("~{pangolin_notes}", "pangolin_notes", df)
    df = load_json_data("~{pangolin_assignment_version}", "pangolin_assignment_version", df)
    df = load_json_data("~{pangolin_versions}", "pangolin_versions", df)
    df = load_json_data("~{pangolin_docker}", "pangolin_docker", df)
    df = load_json_data("~{nextclade_version}", "nextclade_version", df)
    df = load_json_data("~{nextclade_docker}", "nextclade_docker", df)
    df = load_json_data("~{nextclade_ds_tag}", "nextclade_ds_tag", df)
    df = load_json_data("~{nextclade_aa_subs}", "nextclade_aa_subs", df)
    df = load_json_data("~{nextclade_aa_dels}", "nextclade_aa_dels", df)
    df = load_json_data("~{nextclade_clade}", "nextclade_clade", df)
    df = load_json_data("~{nextclade_lineage}", "nextclade_lineage", df)
    df = load_json_data("~{nextclade_qc}", "nextclade_qc", df)
    df = load_json_data("~{nextclade_ds_tag_flu_ha}", "nextclade_ds_tag_flu_ha", df)
    df = load_json_data("~{nextclade_aa_subs_flu_ha}", "nextclade_aa_subs_flu_ha", df)
    df = load_json_data("~{nextclade_aa_dels_flu_ha}", "nextclade_aa_dels_flu_ha", df)
    df = load_json_data("~{nextclade_clade_flu_ha}", "nextclade_clade_flu_ha", df)
    df = load_json_data("~{nextclade_qc_flu_ha}", "nextclade_qc_flu_ha", df)
    df = load_json_data("~{nextclade_ds_tag_flu_na}", "nextclade_ds_tag_flu_na", df)
    df = load_json_data("~{nextclade_aa_subs_flu_na}", "nextclade_aa_subs_flu_na", df)
    df = load_json_data("~{nextclade_aa_dels_flu_na}", "nextclade_aa_dels_flu_na", df)
    df = load_json_data("~{nextclade_clade_flu_na}", "nextclade_clade_flu_na", df)
    df = load_json_data("~{nextclade_qc_flu_na}", "nextclade_qc_flu_na", df)
    
    print(df)
    df.to_csv("~{samplename}.results.tsv", sep='\t', index=False)


    organism_names = df["organism"].replace('', np.nan).dropna()
    organism_names.to_csv("~{samplename}.organism_names.tsv", index=False, header=False)
    CODE
  >>>
  output {
    File gathered_results = "~{samplename}.results.tsv"
    Array[String] organism_names = read_lines("~{samplename}.organism_names.tsv")

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