version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/quality_control/basic_statistics/task_fastq_scan.wdl" as fastq_scan_task
import "../../tasks/taxon_id/task_ete4_taxon_id.wdl" as identify_taxon_id_task
import "../../tasks/taxon_id/task_krakentools.wdl" as krakentools_task
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken2_task
import "../../tasks/quality_control/read_filtering/task_bbduk.wdl" as bbduk_task
import "../../tasks/quality_control/read_filtering/task_fastp.wdl" as fastp_task
import "../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub_task
import "../../tasks/utilities/file_handling/task_cat_lanes.wdl" as cat_lanes_task
import "../../tasks/utilities/data_export/task_export_taxon_table.wdl" as export_taxon_table_task
import "../../tasks/utilities/file_handling/task_kraken_parser.wdl" as kraken_parser_task
import "../utilities/wf_host_decontaminate.wdl" as host_decontaminate_wf
import "wf_theiaviral_illumina_pe.wdl" as theiaviral_illumina_pe


workflow theiaviral_panel {
  input {
    File read1
    File read2
    String samplename
    # Default Taxon IDs 
    Array[String] taxon_ids = ["1618189", "37124", "46839", "12637", "1216928", "59301", "2169701", "118655", "11587", "11029", "11033", "11034", "1608084", "64286", 
      "11082", "11089", "64320", "10804", "12092", "3052230", "12475", "11676", "11709", "68887", "1980456", "3052470", "3052490", "169173", 
      "3052489", "238817", "1980442", "3052496", "3052499", "90961", "80935", "35305", "1221391", "38767", "11021", "38768", "2847089", 
      "3052223", "260964", "35511", "11072", "11577", "38766", "1474807", "12538", "11079", "3052225", "11083", "11292", "11580", "11080", "45270", "11084", 
      "11036", "11039", "1313215", "138948", "138949", "138950", "138951", "1239565", "1239570", "1239573", "142786", 
      "28875", "28876", "36427", "1348384", "1330524", "95341", "2849717", "1424613", "2010960", "565995", "3052302", "3052518", "3052477", "3052307", "3052480", 
      "2169991", "33743", "3052310", "3052148", "3052314", "3052303", "3052317", "33727", "12542", "3052493", "186539", "11588", "2907957", 
      "3052498", "1003835", "1452514", "186540", "186541", "3052503", "1891762", "10376", "10359", "333760", "333761", "337044", "337050", "1671798", "333754", "333767", "746830", "746831", "943908", "10632", "1891764", "1965344", "493803", "1203539", "1497391", "1891767", "1277649", "862909", "862909", "440266", 
      "11234", "152219", "10244", "2560602", "11041", "10335", "10255", "129875", "108098", "129951", "130310", "130308", "130309", "536079", "329641", "11137", "290028", "277944", 
      "31631", "162145", "12730", "2560525", "11216", "2560526", "1803956", "10798", "11250", "11320", "11520", "11552", "1335626", "147711", "147712", "463676", "2901879", 
      "2697049", "10404", "3052505", "337041", "337042", "333757", "337048", "333754", "333766", "337049"
    ]
    File output_taxon_table = "gs://theiagen-public-resources-rp/reference_data/family_agnostic/theiaviral_panel_taxon_table_20251111.tsv"
    String source_table_name
    String terra_project
    String terra_workspace
    File kraken_db = "gs://theiagen-public-resources-rp/reference_data/databases/kraken2/kraken2_humanGRCh38_viralRefSeq_20240828.tar.gz"
    Boolean extract_unclassified = false
    Int min_read_count = 1000
    Boolean call_metaviralspades = true
    String? host
  }  
  call versioning.version_capture {
    input:
  }

  # read QC, classification, extraction, and trimming
  call fastq_scan_task.fastq_scan_pe as fastq_scan_raw {
    input:
      read1 = read1,
      read2 = read2
  }
  # human read scrubbing
  call ncbi_scrub_task.ncbi_scrub_pe {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2
  }
  # adapter + read trimming
  call fastp_task.fastp {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_pe.read1_dehosted,
      read2 = ncbi_scrub_pe.read2_dehosted
  }
  call bbduk_task.bbduk {
    input:
      samplename = samplename,
      read1 = fastp.read1_trimmed,
      read2 = select_first([fastp.read2_trimmed])
  }
  # host read decontamination
  if (defined(host)) {
    call host_decontaminate_wf.host_decontaminate {
      input:
        samplename = samplename,
        read1 = bbduk.read1_clean,
        read2 = bbduk.read2_clean,
        host = select_first([host])
    }
  }
  # taxon-based read extraction
  call kraken2_task.kraken2_standalone as kraken2_standalone_raw {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      kraken2_db = select_first([kraken_db])
  }
  call kraken2_task.kraken2_standalone as kraken2_standalone_clean {
    input:
      samplename = samplename,
      read1 = select_first([host_decontaminate.dehost_read1, bbduk.read1_clean]),
      read2 = select_first([host_decontaminate.dehost_read2, bbduk.read2_clean]),
      kraken2_db = select_first([kraken_db])
  }
  # clean read stat gathering
  call fastq_scan_task.fastq_scan_pe as fastq_scan_clean {
    input:
      read1 = select_first([host_decontaminate.dehost_read1, bbduk.read1_clean]),
      read2 = select_first([host_decontaminate.dehost_read2, bbduk.read2_clean])
  }
  # parse the kraken report to get the taxon ids present in the sample, lowering scatter shards
  call kraken_parser_task.kraken_output_parser as kraken_parser {
    input:
      kraken2_report = select_first([kraken2_standalone_clean.kraken2_report]),
      taxon_ids = taxon_ids,
      read_count_threshold = min_read_count
  }
  # get kraken outputs ready will have to run it outside of read qc trim to fix
  scatter (taxon_id in kraken_parser.parsed_taxon_ids) {
    call krakentools_task.extract_kraken_reads as krakentools {
      input:
        kraken2_output = select_first([kraken2_standalone_clean.kraken2_classified_report]),
        kraken2_report = select_first([kraken2_standalone_clean.kraken2_report]),
        read1 = select_first([kraken2_standalone_clean.kraken2_classified_read1]),
        read2 = select_first([kraken2_standalone_clean.kraken2_classified_read2]),
        taxon_id = taxon_id
    }
    if (krakentools.status == "PASS") {
      if (extract_unclassified) {
        call cat_lanes_task.cat_lanes {
          input:
            samplename = samplename + "_" + taxon_id,
            read1_lane1 = select_first([kraken2_standalone_clean.kraken2_unclassified_read1]),
            read1_lane2 = select_first([krakentools.extracted_read1]),
            read2_lane1 = select_first([kraken2_standalone_clean.kraken2_unclassified_read2]),
            read2_lane2 = select_first([krakentools.extracted_read2])
        }
      }
      call kraken2_task.kraken2_standalone as kraken2 {
        input:
          read1 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
          read2 = select_first([cat_lanes.read2_concatenated, krakentools.extracted_read2]),
          kraken2_db = kraken_db,
          samplename = samplename + "_" + taxon_id
      }
      # get the taxon information from ncbi
      call identify_taxon_id_task.ete4_taxon_id as ete4_identify {
        input:
          taxon = taxon_id
      }
      call theiaviral_illumina_pe.theiaviral_illumina_pe {
        input:
          read1 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
          read2 = select_first([cat_lanes.read2_concatenated, krakentools.extracted_read2]),
          samplename = samplename + "_" + taxon_id,
          taxon = taxon_id,
          call_metaviralspades = call_metaviralspades,
          kraken_db = kraken_db,
          skip_qc = true,
          skip_screen = true,
      }
      # call export_taxon_table
      call export_taxon_table_task.export_taxon_table {
        input:
          samplename = samplename + "_" + sub(ete4_identify.taxon_name, " ", "_"),
          taxon_table = output_taxon_table,
          theiaviral_panel = true,
          gambit_predicted_taxon = ete4_identify.taxon_name,
          terra_project = terra_project,
          terra_workspace = terra_workspace,
          columns_to_export = {
            "samplename": samplename + "_" + taxon_id,
            "source_table": source_table_name,
            "kraken_extracted_taxon_name": ete4_identify.taxon_name,
            "assembly_fasta": theiaviral_illumina_pe.assembly_consensus_fasta,
            "theiaviral_panel_analysis_date": version_capture.date,
            "theiaviral_panel_version": version_capture.phb_version,
            "read1": select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
            "read2": select_first([cat_lanes.read2_concatenated, krakentools.extracted_read2]),
            "kraken2_report": kraken2.kraken2_report,
            "kraken2_classified_report": kraken2.kraken2_classified_report,
            "kraken2_percent_human": kraken2.kraken2_percent_human,
            "kraken2_version": kraken2.kraken2_version,
            "kraken2_docker": kraken2.kraken2_docker,
            "kraken2_database": kraken2.kraken2_database,
            "taxon_avg_genome_length": theiaviral_illumina_pe.taxon_avg_genome_length,
            "datasets_genome_length_docker": theiaviral_illumina_pe.datasets_genome_length_docker,
            "datasets_genome_length_version": theiaviral_illumina_pe.datasets_genome_length_version,
            "assembly_denovo_fasta": theiaviral_illumina_pe.assembly_denovo_fasta,
            "metaviralspades_status": theiaviral_illumina_pe.metaviralspades_status,
            "metaviralspades_version": theiaviral_illumina_pe.metaviralspades_version,
            "metaviralspades_docker": theiaviral_illumina_pe.metaviralspades_docker,
            "megahit_status": theiaviral_illumina_pe.megahit_status,
            "megahit_version": theiaviral_illumina_pe.megahit_version,
            "megahit_docker": theiaviral_illumina_pe.megahit_docker,
            "checkv_denovo_summary": theiaviral_illumina_pe.checkv_denovo_summary,
            "checkv_denovo_contamination": theiaviral_illumina_pe.checkv_denovo_contamination,
            "checkv_denovo_weighted_completeness": theiaviral_illumina_pe.checkv_denovo_weighted_completeness,
            "checkv_denovo_weighted_contamination": theiaviral_illumina_pe.checkv_denovo_weighted_contamination,
            "checkv_denovo_total_genes": theiaviral_illumina_pe.checkv_denovo_total_genes,
            "checkv_denovo_version": theiaviral_illumina_pe.checkv_denovo_version,
            "checkv_denovo_status": theiaviral_illumina_pe.checkv_denovo_status,
            "quast_denovo_report": theiaviral_illumina_pe.quast_denovo_report,
            "quast_denovo_genome_length": theiaviral_illumina_pe.quast_denovo_genome_length,
            "quast_denovo_number_contigs": theiaviral_illumina_pe.quast_denovo_number_contigs,
            "quast_denovo_n50_value": theiaviral_illumina_pe.quast_denovo_n50_value,
            "quast_denovo_largest_contig": theiaviral_illumina_pe.quast_denovo_largest_contig,
            "quast_denovo_gc_percent": theiaviral_illumina_pe.quast_denovo_gc_percent,
            "quast_denovo_uncalled_bases": theiaviral_illumina_pe.quast_denovo_uncalled_bases,
            "quast_denovo_version": theiaviral_illumina_pe.quast_denovo_version,
            "quast_denovo_docker": theiaviral_illumina_pe.quast_denovo_docker,
            "skani_reference_taxon_name": theiaviral_illumina_pe.skani_reference_taxon,
            "skani_report": theiaviral_illumina_pe.skani_report,
            "skani_warning": theiaviral_illumina_pe.skani_warning,
            "skani_status": theiaviral_illumina_pe.skani_status,
            "skani_top_accession": theiaviral_illumina_pe.skani_top_accession,
            "skani_top_query_coverage": theiaviral_illumina_pe.skani_top_query_coverage,
            "skani_top_score": theiaviral_illumina_pe.skani_top_score,
            "skani_top_ani": theiaviral_illumina_pe.skani_top_ani,
            "skani_database": theiaviral_illumina_pe.skani_database,
            "skani_version": theiaviral_illumina_pe.skani_version,
            "skani_docker": theiaviral_illumina_pe.skani_docker,
            "bwa_version": theiaviral_illumina_pe.bwa_version,
            "bwa_samtools_version": theiaviral_illumina_pe.bwa_samtools_version,
            "bwa_read1_aligned": theiaviral_illumina_pe.bwa_read1_aligned,
            "bwa_read2_aligned": theiaviral_illumina_pe.bwa_read2_aligned,
            "bwa_sorted_bam": theiaviral_illumina_pe.bwa_sorted_bam,
            "bwa_sorted_bai": theiaviral_illumina_pe.bwa_sorted_bai,
            "bwa_read1_unaligned": theiaviral_illumina_pe.bwa_read1_unaligned,
            "bwa_read2_unaligned": theiaviral_illumina_pe.bwa_read2_unaligned,
            "bwa_sorted_bam_unaligned": theiaviral_illumina_pe.bwa_sorted_bam_unaligned,
            "bwa_sorted_bai_unaligned": theiaviral_illumina_pe.bwa_sorted_bam_unaligned_bai,
            "ivar_tsv": theiaviral_illumina_pe.ivar_tsv,
            "ivar_vcf": theiaviral_illumina_pe.ivar_vcf,
            "ivar_variant_version": theiaviral_illumina_pe.ivar_variant_version,
            "read_mapping_report": theiaviral_illumina_pe.read_mapping_report,
            "read_mapping_statistics": theiaviral_illumina_pe.read_mapping_statistics,
            "read_mapping_cov_hist": theiaviral_illumina_pe.read_mapping_cov_hist,
            "read_mapping_cov_stats": theiaviral_illumina_pe.read_mapping_cov_stats,
            "read_mapping_flagstat": theiaviral_illumina_pe.read_mapping_flagstat,
            "read_mapping_coverage": theiaviral_illumina_pe.read_mapping_coverage,
            "read_mapping_depth": theiaviral_illumina_pe.read_mapping_depth,
            "read_mapping_meanbaseq": theiaviral_illumina_pe.read_mapping_meanbaseq,
            "read_mapping_meanmapq": theiaviral_illumina_pe.read_mapping_meanmapq,
            "read_mapping_percentage_mapped_reads": theiaviral_illumina_pe.read_mapping_percentage_mapped_reads,
            "read_mapping_samtools_version": theiaviral_illumina_pe.read_mapping_samtools_version,
            "consensus_qc_number_N": theiaviral_illumina_pe.consensus_qc_number_N,
            "consensus_qc_assembly_length_unambiguous": theiaviral_illumina_pe.consensus_qc_assembly_length_unambiguous,
            "consensus_qc_number_degenerate": theiaviral_illumina_pe.consensus_qc_number_Degenerate,
            "consensus_qc_number_total": theiaviral_illumina_pe.consensus_qc_number_Total,
            "consensus_qc_percent_reference_coverage": theiaviral_illumina_pe.consensus_qc_percent_reference_coverage,
            "checkv_consensus_summary": theiaviral_illumina_pe.checkv_consensus_summary,
            "checkv_consensus_contamination": theiaviral_illumina_pe.checkv_consensus_contamination,
            "checkv_consensus_weighted_completeness": theiaviral_illumina_pe.checkv_consensus_weighted_completeness,
            "checkv_consensus_weighted_contamination": theiaviral_illumina_pe.checkv_consensus_weighted_contamination,
            "checkv_consensus_total_genes": theiaviral_illumina_pe.checkv_consensus_total_genes,
            "checkv_consensus_version": theiaviral_illumina_pe.checkv_consensus_version,
            "checkv_consensus_status": theiaviral_illumina_pe.checkv_consensus_status,
            "pango_lineage": theiaviral_illumina_pe.pango_lineage,
            "pango_lineage_expanded": theiaviral_illumina_pe.pango_lineage_expanded,
            "pangolin_conflicts": theiaviral_illumina_pe.pangolin_conflicts,
            "pangolin_notes": theiaviral_illumina_pe.pangolin_notes,
            "pangolin_assignment_version": theiaviral_illumina_pe.pangolin_assignment_version,
            "pango_lineage_report": theiaviral_illumina_pe.pango_lineage_report,
            "pangolin_docker": theiaviral_illumina_pe.pangolin_docker,
            "pangolin_versions": theiaviral_illumina_pe.pangolin_versions,
            "nextclade_version": theiaviral_illumina_pe.nextclade_version,
            "nextclade_docker": theiaviral_illumina_pe.nextclade_docker,
            "nextclade_json": theiaviral_illumina_pe.nextclade_json,
            "auspice_json": theiaviral_illumina_pe.auspice_json,
            "nextclade_tsv": theiaviral_illumina_pe.nextclade_tsv,
            "nextclade_ds_tag": theiaviral_illumina_pe.nextclade_ds_tag,
            "nextclade_aa_subs": theiaviral_illumina_pe.nextclade_aa_subs,
            "nextclade_aa_dels": theiaviral_illumina_pe.nextclade_aa_dels,
            "nextclade_clade": theiaviral_illumina_pe.nextclade_clade,
            "nextclade_lineage": theiaviral_illumina_pe.nextclade_lineage,
            "nextclade_qc": theiaviral_illumina_pe.nextclade_qc,
            "nextclade_json_flu_ha": theiaviral_illumina_pe.nextclade_json_flu_ha,
            "auspice_json_flu_ha": theiaviral_illumina_pe.auspice_json_flu_ha,
            "nextclade_tsv_flu_ha": theiaviral_illumina_pe.nextclade_tsv_flu_ha,
            "nextclade_ds_tag_flu_ha": theiaviral_illumina_pe.nextclade_ds_tag_flu_ha,
            "nextclade_aa_subs_flu_ha": theiaviral_illumina_pe.nextclade_aa_subs_flu_ha,
            "nextclade_aa_dels_flu_ha": theiaviral_illumina_pe.nextclade_aa_dels_flu_ha,
            "nextclade_clade_flu_ha": theiaviral_illumina_pe.nextclade_clade_flu_ha,
            "nextclade_qc_flu_ha": theiaviral_illumina_pe.nextclade_qc_flu_ha,
            "nextclade_json_flu_na": theiaviral_illumina_pe.nextclade_json_flu_na,
            "auspice_json_flu_na": theiaviral_illumina_pe.auspice_json_flu_na,
            "nextclade_tsv_flu_na": theiaviral_illumina_pe.nextclade_tsv_flu_na,
            "nextclade_ds_tag_flu_na": theiaviral_illumina_pe.nextclade_ds_tag_flu_na,
            "nextclade_aa_subs_flu_na": theiaviral_illumina_pe.nextclade_aa_subs_flu_na,
            "nextclade_aa_dels_flu_na": theiaviral_illumina_pe.nextclade_aa_dels_flu_na,
            "nextclade_clade_flu_na": theiaviral_illumina_pe.nextclade_clade_flu_na,
            "nextclade_qc_flu_na": theiaviral_illumina_pe.nextclade_qc_flu_na,
            "irma_version": theiaviral_illumina_pe.irma_version,
            "irma_docker": theiaviral_illumina_pe.irma_docker,
            "irma_type": theiaviral_illumina_pe.irma_type,
            "irma_subtype": theiaviral_illumina_pe.irma_subtype,
            "irma_subtype_notes": theiaviral_illumina_pe.irma_subtype_notes,
            "irma_read1_aligned": theiaviral_illumina_pe.irma_read1_aligned,
            "irma_read2_aligned": theiaviral_illumina_pe.irma_read2_aligned,
            "genoflu_version": theiaviral_illumina_pe.genoflu_version,
            "genoflu_genotype": theiaviral_illumina_pe.genoflu_genotype,
            "genoflu_all_segments": theiaviral_illumina_pe.genoflu_all_segments,
            "genoflu_output_tsv": theiaviral_illumina_pe.genoflu_output_tsv,
            "abricate_flu_type": theiaviral_illumina_pe.abricate_flu_type,
            "abricate_flu_subtype":  theiaviral_illumina_pe.abricate_flu_subtype,
            "abricate_flu_results": theiaviral_illumina_pe.abricate_flu_results,
            "abricate_flu_database":  theiaviral_illumina_pe.abricate_flu_database,
            "abricate_flu_version": theiaviral_illumina_pe.abricate_flu_version
          }
      }
    }  
  }
  output {
    # Number of assembled viruses
    Int assembled_viruses = length(select_all(theiaviral_illumina_pe.assembly_consensus_fasta))
    # Taxon of organisms identified
    Array[String] identified_organisms = select_all(ete4_identify.taxon_name)
    # All assembled FASTA files
    Array[File] assemblies = select_all(theiaviral_illumina_pe.assembly_consensus_fasta)
    # Workflow Versioning
    String theiaviral_panel_version = version_capture.phb_version
    String theiaviral_panel_analysis_date = version_capture.date
    # Standalone Kraken2 outputs
    String kraken2_version = kraken2_standalone_raw.kraken2_version
    String kraken2_database = kraken2_standalone_raw.kraken2_database
    String kraken2_docker = kraken2_standalone_raw.kraken2_docker
    File kraken2_report_raw = kraken2_standalone_raw.kraken2_report
    Float? kraken2_percent_human_raw = kraken2_standalone_clean.kraken2_percent_human
    File kraken2_report_clean = select_first([kraken2_standalone_clean.kraken2_report])
    File kraken2_classified_report = select_first([kraken2_standalone_clean.kraken2_classified_report])
    Float? kraken2_percent_human_clean = kraken2_standalone_clean.kraken2_percent_human
    # raw read quality control
    Int? fastq_scan_num_reads_raw1 = fastq_scan_raw.read1_seq
    Int? fastq_scan_num_reads_raw2 = fastq_scan_raw.read2_seq
    String? fastq_scan_raw_pairs = fastq_scan_raw.read_pairs
    String? fastq_scan_version = fastq_scan_raw.version
    String? fastq_scan_docker = fastq_scan_raw.fastq_scan_docker
    File? fastq_scan_raw1_json = fastq_scan_raw.read1_fastq_scan_json
    File? fastq_scan_raw2_json = fastq_scan_raw.read2_fastq_scan_json
    # trimming data
    String? fastp_version = fastp.fastp_version
    String? fastp_docker = fastp.fastp_docker
    File? fastp_html_report = fastp.fastp_stats_html
    File? fastp_json_report = fastp.fastp_stats_json
    # host decontamination outputs
    File? dehost_wf_dehost_read1 = host_decontaminate.dehost_read1
    File? dehost_wf_dehost_read2 = host_decontaminate.dehost_read2
    String? dehost_wf_host_accession = host_decontaminate.host_genome_accession
    File? dehost_wf_host_mapped_bam = host_decontaminate.host_mapped_sorted_bam
    File? dehost_wf_host_mapped_bai = host_decontaminate.host_mapped_sorted_bai
    File? dehost_wf_host_fasta = host_decontaminate.host_genome_fasta
    File? dehost_wf_host_mapping_stats = host_decontaminate.host_mapping_stats
    File? dehost_wf_host_mapping_cov_hist = host_decontaminate.host_mapping_cov_hist
    File? dehost_wf_host_flagstat = host_decontaminate.host_flagstat
    Float? dehost_wf_host_mapping_coverage = host_decontaminate.host_mapping_coverage
    Float? dehost_wf_host_mapping_mean_depth = host_decontaminate.host_mapping_mean_depth
    Float? dehost_wf_host_percent_mapped_reads = host_decontaminate.host_percent_mapped_reads
    File? dehost_wf_host_mapping_metrics = host_decontaminate.host_mapping_metrics
    # NCBI scrubber
    File? ncbi_scrub_read1_dehosted = ncbi_scrub_pe.read1_dehosted
    File? ncbi_scrub_read2_dehosted = ncbi_scrub_pe.read2_dehosted
    Int? ncbi_scrub_human_spots_removed = ncbi_scrub_pe.human_spots_removed
    String? ncbi_scrub_docker = ncbi_scrub_pe.ncbi_scrub_docker
    # bbduk
    File read1_clean = bbduk.read1_clean
    File read2_clean = bbduk.read2_clean
    String bbduk_docker = bbduk.bbduk_docker
    # clean read quality control
    Int? fastq_scan_num_reads_clean1 = fastq_scan_clean.read1_seq
    Int? fastq_scan_num_reads_clean2 = fastq_scan_clean.read2_seq
    String? fastq_scan_clean_pairs = fastq_scan_clean.read_pairs
    File? fastq_scan_clean1_json = fastq_scan_clean.read1_fastq_scan_json
    File? fastq_scan_clean2_json = fastq_scan_clean.read2_fastq_scan_json
  }
}