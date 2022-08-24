version 1.0

task export_taxon_tables {
  input {
    String terra_project
    String terra_workspace
    String sample_taxon
    File? taxon_tables
    String samplename 
    # TheiaProk Outputs
    File? reads
    File? read1
    File? read2
    File? read1_clean
    File? read2_clean
    String? run_id
    String? collection_date
    String? originating_lab
    String? city
    String? county
    String? zip
    String? theiaprok_illumina_pe_version
    String? theiaprok_illumina_pe_analysis_date
    String? theiaprok_illumina_se_version
    String? theiaprok_illumina_se_analysis_date
    String seq_platform
    Int num_reads_raw1
    Int? num_reads_raw2
    String? num_reads_raw_pairs
    String fastq_scan_version
    Int num_reads_clean1
    Int? num_reads_clean2
    String? num_reads_clean_pairs
    String trimmomatic_version
    String bbduk_docker
    Float r1_mean_q
    Float? r2_mean_q
    File assembly_fasta
    File? contigs_gfa
    String? shovill_pe_version
    String? shovill_se_version
    File quast_report
    String quast_version
    Int genome_length
    Int number_contigs
    Int n50_value
    File cg_pipeline_report
    String cg_pipeline_docker
    Float est_coverage
    File gambit_report
    String gambit_predicted_taxon
    String gambit_predicted_taxon_rank
    File gambit_closest_genomes
    String gambit_version
    String gambit_db_version
    String gambit_docker
    File amrfinderplus_all_report
    File amrfinderplus_amr_report
    File amrfinderplus_stress_report
    File amrfinderplus_virulence_report
    String amrfinderplus_amr_genes
    String amrfinderplus_stress_genes
    String amrfinderplus_virulence_genes
    String amrfinderplus_version
    String amrfinderplus_db_version
    String ts_mlst_results
    String ts_mlst_predicted_st
    String ts_mlst_pubmlst_scheme
    String ts_mlst_version
    File? serotypefinder_report
    String? serotypefinder_docker
    String? serotypefinder_serotype
    File? ectyper_results
    String? ectyper_version
    String? ectyper_predicted_serotype
    File? lissero_results
    String? lissero_version
    File? sistr_results
    File? sistr_allele_json
    File? sister_allele_fasta
    File? sistr_cgmlst
    String? sistr_version
    String? sistr_predicted_serotype
    File? seqsero2_report
    String? seqsero2_version
    String? seqsero2_predicted_antigenic_profile
    String? seqsero2_predicted_serotype
    String? seqsero2_predicted_contamination
    File? genotyphi_report_tsv
    File? genotyphi_mykrobe_json
    String? genotyphi_version
    String? genotyphi_species
    Float? genotyphi_st_probes_percent_coverage
    String? genotyphi_final_genotype
    String? genotyphi_genotype_confidence
    File? kleborate_output_file
    String? kleborate_version
    String? kleborate_key_resistance_genes
    String? kleborate_genomic_resistance_mutations
    String? kleborate_mlst_sequence_type
    File? tbprofiler_output_file
    File? tbprofiler_output_bam
    File? tbprofiler_output_bai
    String? tbprofiler_version
    String? tbprofiler_main_lineage
    String? tbprofiler_sub_lineage
    String? tbprofiler_dr_type
    String? tbprofiler_resistance_genes
    File? legsta_results
    String? legsta_predicted_sbt
    String? legsta_version
    
  }
  command <<<
  
    # capture taxon and corresponding table names from input taxon_tables
    taxon_array=($(cut -f1 ~{taxon_tables} | tail +2))
    echo "Taxon array: ${taxon_array[*]}"
    table_array=($(cut -f2 ~{taxon_tables} | tail +2))
    echo "Table array: ${table_array[*]}"
    # remove whitespace from sample_taxon
    sample_taxon=$(echo ~{sample_taxon} | tr ' ' '_')
    echo "Sample taxon: ${sample_taxon}"
    # set taxon and table vars
    echo "Checking if sample taxon should be exported to user-specified taxon table..."
    for index in ${!taxon_array[@]}; do
      taxon=${taxon_array[$index]}
      echo "Taxon: ${taxon}"
      table=${table_array[$index]}
      echo "Table: ${table}"
      if [[ "${sample_taxon}" == *"${taxon}"* ]]; then
        sample_table=${table}
        echo "Sample ~{samplename} identified as ~{sample_taxon}. As per user-defined taxon tables, ${taxon} samples will be exported to the ${table} terra data table"
        break
      else 
        echo "${sample_taxon} does not match ${taxon}."
      fi
    done
    if [ ! -z ${sample_table} ]; then
       # create single-entity sample data table
       ## header
      echo -e "entity:${sample_table}_id\treads\tread1\tread2\tread1_clean\tread2_clean\trun_id\tcollection_date\toriginating_lab\tcity\tcounty\tzip\ttheiaprok_illumina_pe_version\ttheiaprok_illumina_pe_analysis_date\ttheiaprok_illumina_se_version\ttheiaprok_illumina_se_analysis_date\tseq_platform\tnum_reads_raw1\tnum_reads_raw2\tnum_reads_raw_pairs\tfastq_scan_version\tnum_reads_clean1\tnum_reads_clean2\tnum_reads_clean_pairs\ttrimmomatic_version\tbbduk_docker\tr1_mean_q\tr2_mean_q\tassembly_fasta\tcontigs_gfa\tshovill_pe_version\tshovill_se_version\tquast_report\tquast_version\tgenome_length\tnumber_contigs\tn50_value\tcg_pipeline_report\tcg_pipeline_docker\test_coverage\tgambit_report\tgambit_predicted_taxon\tgambit_predicted_taxon_rank\tgambit_closest_genomes\tgambit_version\tgambit_db_version\tgambit_docker\tts_mlst_results\tts_mlst_predicted_st\tts_mlst_pubmlst_scheme\tts_mlst_version\tserotypefinder_report\tserotypefinder_docker\tserotypefinder_serotype\tectyper_results\tectyper_version\tectyper_predicted_serotype\tlissero_results\tlissero_version\tsistr_results\tsistr_allele_json\tsister_allele_fasta\tsistr_cgmlst\tsistr_version\tsistr_predicted_serotype\tseqsero2_report\tseqsero2_version\tseqsero2_predicted_antigenic_profile\tseqsero2_predicted_serotype\tseqsero2_predicted_contamination\tkleborate_output_file\tkleborate_version\tkleborate_key_resistance_genes\tkleborate_genomic_resistance_mutations\tkleborate_mlst_sequence_type\tlegsta_results\tlegsta_predicted_sbt\tlegsta_version\ttbprofiler_output_file\ttbprofiler_output_bam\ttbprofiler_output_bai\ttbprofiler_version\ttbprofiler_main_lineage\ttbprofiler_sub_lineage\ttbprofiler_dr_type\ttbprofiler_resistance_genes\tamrfinderplus_all_report\tamrfinderplus_amr_report\tamrfinderplus_stress_report\tamrfinderplus_virulence_report\tamrfinderplus_version\tamrfinderplus_db_version\tamrfinderplus_amr_genes\tamrfinderplus_stress_genes\tamrfinderplus_virulence_genes\tgenotyphi_report_tsv\tgenotyphi_mykrobe_json\tgenotyphi_version\tgenotyphi_species\tgenotyphi_st_probes_percent_coverage\tgenotyphi_final_genotype\tgenotyphi_genotype_confidence" > ~{samplename}_terra_table.tsv
      ## TheiaProk Outs
      echo -e "~{samplename}\t~{reads}\t~{read1}\t~{read2}\t~{read1_clean}\t~{read2_clean}\t~{run_id}\t~{collection_date}\t~{originating_lab}\t~{city}\t~{county}\t~{zip}\t~{theiaprok_illumina_pe_version}\t~{theiaprok_illumina_pe_analysis_date}\t~{theiaprok_illumina_se_version}\t~{theiaprok_illumina_se_analysis_date}\t~{seq_platform}\t~{num_reads_raw1}\t~{num_reads_raw2}\t~{num_reads_raw_pairs}\t~{fastq_scan_version}\t~{num_reads_clean1}\t~{num_reads_clean2}\t~{num_reads_clean_pairs}\t~{trimmomatic_version}\t~{bbduk_docker}\t~{r1_mean_q}\t~{r2_mean_q}\t~{assembly_fasta}\t~{contigs_gfa}\t~{shovill_pe_version}\t~{shovill_se_version}\t~{quast_report}\t~{quast_version}\t~{genome_length}\t~{number_contigs}\t~{n50_value}\t~{cg_pipeline_report}\t~{cg_pipeline_docker}\t~{est_coverage}\t~{gambit_report}\t~{gambit_predicted_taxon}\t~{gambit_predicted_taxon_rank}\t~{gambit_closest_genomes}\t~{gambit_version}\t~{gambit_db_version}\t~{gambit_docker}\t~{ts_mlst_results}\t~{ts_mlst_predicted_st}\t~{ts_mlst_pubmlst_scheme}\t~{ts_mlst_version}\t~{serotypefinder_report}\t~{serotypefinder_docker}\t~{serotypefinder_serotype}\t~{ectyper_results}\t~{ectyper_version}\t~{ectyper_predicted_serotype}\t~{lissero_results}\t~{lissero_version}\t~{sistr_results}\t~{sistr_allele_json}\t~{sister_allele_fasta}\t~{sistr_cgmlst}\t~{sistr_version}\t~{sistr_predicted_serotype}\t~{seqsero2_report}\t~{seqsero2_version}\t~{seqsero2_predicted_antigenic_profile}\t~{seqsero2_predicted_serotype}\t~{seqsero2_predicted_contamination}\t~{kleborate_output_file}\t~{kleborate_version}\t~{kleborate_key_resistance_genes}\t~{kleborate_genomic_resistance_mutations}\t~{kleborate_mlst_sequence_type}\t~{legsta_results}\t~{legsta_predicted_sbt}\t~{legsta_version}\t~{tbprofiler_output_file}\t~{tbprofiler_output_bam}\t~{tbprofiler_output_bai}\t~{tbprofiler_version}\t~{tbprofiler_main_lineage}\t~{tbprofiler_sub_lineage}\t~{tbprofiler_dr_type}\t~{tbprofiler_resistance_genes}\t~{amrfinderplus_all_report}\t~{amrfinderplus_amr_report}\t~{amrfinderplus_stress_report}\t~{amrfinderplus_virulence_report}\t~{amrfinderplus_version}\t~{amrfinderplus_db_version}\t~{amrfinderplus_amr_genes}\t~{amrfinderplus_stress_genes}\t~{amrfinderplus_virulence_genes}\t~{genotyphi_report_tsv}\t~{genotyphi_mykrobe_json}\t~{genotyphi_version}\t~{genotyphi_species}\t~{genotyphi_st_probes_percent_coverage}\t~{genotyphi_final_genotype}\t~{genotyphi_genotype_confidence}"  >> ~{samplename}_terra_table.tsv
      # modify file paths to GCP URIs
      sed -i 's/\/cromwell_root\//gs:\/\//g' ~{samplename}_terra_table.tsv
      # export table 
      python3 /scripts/import_large_tsv/import_large_tsv.py --project "~{terra_project}" --workspace "~{terra_workspace}" --tsv ~{samplename}_terra_table.tsv
    else
      echo "Table not defined for ~{sample_taxon}"
    fi
  >>>
  runtime {
    docker: "broadinstitute/terra-tools:tqdm"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk 100 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
  output {
    File? datatable1_tsv = "~{samplename}_terra_table.tsv"
  }
}
