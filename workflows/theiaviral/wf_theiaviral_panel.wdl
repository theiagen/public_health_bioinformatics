version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc
import "../../tasks/assembly/task_spades.wdl" as spades_task
import "../../tasks/assembly/task_megahit.wdl" as megahit_task
import "../../tasks/quality_control/advanced_metrics/task_checkv.wdl" as checkv_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/taxon_id/task_skani.wdl" as skani_task
import "../../tasks/utilities/data_import/task_ncbi_datasets.wdl" as ncbi_datasets_task
import "../../tasks/taxon_id/task_identify_taxon_id.wdl" as identify_taxon_id_task
import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/assembly/task_ivar_consensus.wdl" as ivar_consensus
import "../../tasks/gene_typing/variant_detection/task_ivar_variant_call.wdl" as variant_call_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../utilities/wf_morgana_magic.wdl" as morgana_magic_wf
import "../../tasks/taxon_id/task_krakentools.wdl" as krakentools_task
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken2_task
import "../../tasks/quality_control/basic_statistics/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/utilities/file_handling/task_cat_lanes.wdl" as cat_lanes
import "../../tasks/utilities/data_export/task_export_taxon_table.wdl" as export_taxon_table_task
import "wf_theiaviral_illumina_pe.wdl" as theiaviral_illumina_pe

workflow theiaviral_panel {
  input {
    File read1
    File read2
    String samplename
    # Default Taxon IDs for compatibility with VSP2 
    Array[String] taxon_ids = ["1618189", "37124", "46839", "12637", "1216928", "59301", "2169701", "118655", "11587", "11029", "28292", "11033", "11034", "1608084", "64286", 
      "11082", "11089", "64320", "10804", "12092", "10407", "3052230", "12475", "291484", "11676", "11709", "68887", "1980456", "308159", "3052470", "3052490", "169173", "42097", 
      "261204", "430511", "3052489", "238817", "2259728", "47301", "1980442", "3052496", "3052499", "90961", "80935", "35305", "1221391", "38767", "11021", "38768", "2847089", 
      "3052223", "260964", "35511", "11072", "11577", "38766", "2971765", "1474807", "12538", "11079", "3052225", "11083", "11292", "11580", "11080", "45270", "11084", "11590", 
      "11036", "11039", "1313215", "33758", "138948", "138949", "138950", "138951", "1239565", "1239570", "1239572", "1239573", "142786", 
      "28875", "28876", "36427", "1348384", "1330524", "95341", "2849717", "1424613", "2010960", "565995", "3052302", "3052518", "3052477", "1570291", "3052307", "3052480", 
      "2169991", "33743", "3052310", "3052148", "3052314", "3052303", "3052317", "1708253", "33727", "1708252", "12542", "3052493", "378809", "186539", "11588", "2907957", 
      "3052498", "1003835", "1452514", "186540", "186541", "3052503", "1891762", "10376", "10359", "10580", "333760", "333761", "337044", "10585", "10586", "10587", "10588", 
      "10615", "10590", "10591", "10592", "10593", "10595", "10618", "337050", "1671798", "10596", "10598", "37115", "333754", "333767", "37119", "45240", "37121", "39457", "51033", 
      "129724", "746830", "746831", "943908", "10632", "1891764", "1965344", "493803", "1203539", "1497391", "1891767", "1277649", "862909", "862909", "440266", "10298", "10310", 
      "11234", "152219", "10244", "2560602", "11041", "10335", "10255", "129875", "108098", "129951", "130310", "130308", "130309", "536079", "329641", "11137", "290028", "277944", 
      "31631", "162145", "12730", "2560525", "11216", "2560526", "1803956", "10798", "208893", "208895", "11320", "11520", "11552", "1335626", "147711", "147712", "463676", "2901879", 
      "2697049", "10404"
    ]
    File output_taxon_table
    String source_table_name

    String terra_project
    String terra_workspace
    String kraken_db = "gs://theiagen-public-resources-rp/reference_data/databases/kraken2/kraken2_humanGRCh38_viralRefSeq_20240828.tar.gz"
    File? skani_db
    File? checkv_db
    Int min_map_quality = 20
    Int min_depth = 10
    Float min_allele_freq = 0.6
    String? read_extraction_rank
    Boolean concatenate_unclassified = false
    Boolean skip_theiaviral_screen = true
    Int min_read_count = 1000
    Boolean call_metaviralspades = true
  }  
  call versioning.version_capture {
    input:
  }
  # read QC, classification, extraction, and trimming
  call read_qc.read_QC_trim_pe as read_QC_trim {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      kraken_db = kraken_db,
      workflow_series = "theiaviral_panel"
  }
  # get kraken outputs ready will have to run it outside of read qc trim to fix
  scatter (taxon_id in taxon_ids) {
    call krakentools_task.extract_kraken_reads as krakentools {
      input:
        kraken2_output = select_first([read_QC_trim.kraken_classified_report]),
        kraken2_report = select_first([read_QC_trim.kraken_report_clean]),
        read1 = select_first([read_QC_trim.kraken_classified_read1]),
        read2 = select_first([read_QC_trim.kraken_classified_read2]),
        taxon_id = taxon_id
    }
    if (krakentools.success) {
      if (concatenate_unclassified) {
        call cat_lanes.cat_lanes {
          input:
            samplename = samplename + "_" + taxon_id,
            read1_lane1 = select_first([read_QC_trim.kraken_unclassified_read1]),
            read1_lane2 = select_first([krakentools.extracted_read1]),
            read2_lane1 = select_first([read_QC_trim.kraken_unclassified_read2]),
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
      call fastq_scan.fastq_scan_pe as fastq_scan_binned {
        input:
          read1 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
          read2 = select_first([cat_lanes.read2_concatenated, krakentools.extracted_read2])
      }
      if (fastq_scan_binned.read1_seq > min_read_count) {
        # get the taxon information from ncbi
        call identify_taxon_id_task.identify_taxon_id as ncbi_identify {
          input:
            taxon = taxon_id,
            rank = read_extraction_rank,
            use_ncbi_virus = true
        }
        call theiaviral_illumina_pe.theiaviral_illumina_pe as theiaviral {
          input:
            read1 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
            read2 = select_first([cat_lanes.read2_concatenated, krakentools.extracted_read2]),
            samplename = samplename + "_" + taxon_id,
            taxon = taxon_id,
            call_metaviralspades = call_metaviralspades,
            kraken_db = kraken_db,
            skip_qc = true,
            skip_screen = skip_theiaviral_screen,
            skani_db = skani_db,
            checkv_db = checkv_db,
            genome_length = ncbi_identify.avg_genome_length,
            min_map_quality = min_map_quality,
            min_depth = min_depth,
            min_allele_freq = min_allele_freq,
            read_extraction_rank = read_extraction_rank
        }
        if (defined(output_taxon_table)){
          # call export_taxon_table
          call export_taxon_table_task.export_taxon_table_vsp {
            input:
              samplename = samplename + "_" + sub(ncbi_identify.taxon_name, " ", "-"),
              taxon_table = output_taxon_table,
              gambit_predicted_taxon = ncbi_identify.taxon_name,
              terra_project = terra_project,
              terra_workspace = terra_workspace,
              columns_to_export = {
                "samplename": samplename + "_" + taxon_id,
                "source_table": source_table_name,
                "taxon_name": ncbi_identify.taxon_name,
                "assembly_fasta": theiaviral.assembly_consensus_fasta,
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
                "fastq_scan_num_reads_binned1": fastq_scan_binned.read1_seq,
                "fastq_scan_num_reads_binned2": fastq_scan_binned.read2_seq,
                "fastq_scan_num_reads_binned_pairs": fastq_scan_binned.read_pairs,
                "fastq_scan_docker": fastq_scan_binned.fastq_scan_docker,
                "fastq_scan_version": fastq_scan_binned.version,
                "fastq_scan_binned1_json": fastq_scan_binned.read1_fastq_scan_json,
                "fastq_scan_binned2_json": fastq_scan_binned.read2_fastq_scan_json,
                "ncbi_taxon_summary_tsv": ncbi_identify.taxon_summary_tsv,
                "ncbi_taxon_summary_avg_genome_length": ncbi_identify.avg_genome_length,
                "assembly_denovo_fasta": theiaviral.assembly_denovo_fasta,
                "metaviralspades_status": theiaviral.metaviralspades_status,
                "metaviralspades_version": theiaviral.metaviralspades_version,
                "metaviralspades_docker": theiaviral.metaviralspades_docker,
                "megahit_status": theiaviral.megahit_status,
                "megahit_version": theiaviral.megahit_version,
                "megahit_docker": theiaviral.megahit_docker,
                "checkv_denovo_summary": theiaviral.checkv_denovo_summary,
                "checkv_denovo_contamination": theiaviral.checkv_denovo_contamination,
                "checkv_denovo_weighted_completeness": theiaviral.checkv_denovo_weighted_completeness,
                "checkv_denovo_weighted_contamination": theiaviral.checkv_denovo_weighted_contamination,
                "checkv_denovo_total_genes": theiaviral.checkv_denovo_total_genes,
                "checkv_denovo_version": theiaviral.checkv_denovo_version,
                "quast_denovo_report": theiaviral.quast_denovo_report,
                "quast_denovo_genome_length": theiaviral.quast_denovo_genome_length,
                "quast_denovo_number_contigs": theiaviral.quast_denovo_number_contigs,
                "quast_denovo_n50_value": theiaviral.quast_denovo_n50_value,
                "quast_denovo_largest_contig": theiaviral.quast_denovo_largest_contig,
                "quast_denovo_gc_percent": theiaviral.quast_denovo_gc_percent,
                "quast_denovo_uncalled_bases": theiaviral.quast_denovo_uncalled_bases,
                "quast_denovo_version": theiaviral.quast_denovo_version,
                "quast_denovo_docker": theiaviral.quast_denovo_docker,
                "skani_report": theiaviral.skani_report,
                "skani_warning": theiaviral.skani_warning,
                "skani_status": theiaviral.skani_status,
                "skani_top_accession": theiaviral.skani_top_accession,
                "skani_top_query_coverage": theiaviral.skani_top_query_coverage,
                "skani_top_score": theiaviral.skani_top_score,
                "skani_top_ani": theiaviral.skani_top_ani,
                "skani_database": theiaviral.skani_database,
                "skani_version": theiaviral.skani_version,
                "skani_docker": theiaviral.skani_docker,
                "skani_top_ani_fasta": theiaviral.skani_top_ani_fasta,
                "bwa_version": theiaviral.bwa_version,
                "bwa_samtools_version": theiaviral.bwa_samtools_version,
                "bwa_read1_aligned": theiaviral.bwa_read1_aligned,
                "bwa_read2_aligned": theiaviral.bwa_read2_aligned,
                "bwa_sorted_bam": theiaviral.bwa_sorted_bam,
                "bwa_sorted_bai": theiaviral.bwa_sorted_bai,
                "bwa_read1_unaligned": theiaviral.bwa_read1_unaligned,
                "bwa_read2_unaligned": theiaviral.bwa_read2_unaligned,
                "bwa_sorted_bam_unaligned": theiaviral.bwa_sorted_bam_unaligned,
                "bwa_sorted_bai_unaligned": theiaviral.bwa_sorted_bam_unaligned_bai,
                "ivar_tsv": theiaviral.ivar_tsv,
                "ivar_vcf": theiaviral.ivar_vcf,
                "ivar_variant_version": theiaviral.ivar_variant_version,
                "read_mapping_report": theiaviral.read_mapping_report,
                "read_mapping_statistics": theiaviral.read_mapping_statistics,
                "read_mapping_cov_hist": theiaviral.read_mapping_cov_hist,
                "read_mapping_cov_stats": theiaviral.read_mapping_cov_stats,
                "read_mapping_flagstat": theiaviral.read_mapping_flagstat,
                "read_mapping_coverage": theiaviral.read_mapping_coverage,
                "read_mapping_depth": theiaviral.read_mapping_depth,
                "read_mapping_meanbaseq": theiaviral.read_mapping_meanbaseq,
                "read_mapping_meanmapq": theiaviral.read_mapping_meanmapq,
                "read_mapping_percentage_mapped_reads": theiaviral.read_mapping_percentage_mapped_reads,
                "read_mapping_samtools_version": theiaviral.read_mapping_samtools_version,
                "consensus_qc_number_N": theiaviral.consensus_qc_number_N,
                "consensus_qc_assembly_length_unambiguous": theiaviral.consensus_qc_assembly_length_unambiguous,
                "consensus_qc_number_degenerate": theiaviral.consensus_qc_number_Degenerate,
                "consensus_qc_number_total": theiaviral.consensus_qc_number_Total,
                "consensus_qc_percent_reference_coverage": theiaviral.consensus_qc_percent_reference_coverage,
                "checkv_consensus_summary": theiaviral.checkv_consensus_summary,
                "checkv_consensus_contamination": theiaviral.checkv_consensus_contamination,
                "checkv_consensus_weighted_completeness": theiaviral.checkv_consensus_weighted_completeness,
                "checkv_consensus_weighted_contamination": theiaviral.checkv_consensus_weighted_contamination,
                "checkv_consensus_total_genes": theiaviral.checkv_consensus_total_genes,
                "checkv_consensus_version": theiaviral.checkv_consensus_version,
                "pango_lineage": theiaviral.pango_lineage,
                "pango_lineage_expanded": theiaviral.pango_lineage_expanded,
                "pangolin_conflicts": theiaviral.pangolin_conflicts,
                "pangolin_notes": theiaviral.pangolin_notes,
                "pangolin_assignment_version": theiaviral.pangolin_assignment_version,
                "pango_lineage_report": theiaviral.pango_lineage_report,
                "pangolin_docker": theiaviral.pangolin_docker,
                "pangolin_versions": theiaviral.pangolin_versions,
                "nextclade_version": theiaviral.nextclade_version,
                "nextclade_docker": theiaviral.nextclade_docker,
                "nextclade_json_mpxv": theiaviral.nextclade_json_mpxv,
                "auspice_json_mpxv": theiaviral.auspice_json_mpxv,
                "nextclade_tsv_mpxv": theiaviral.nextclade_tsv_mpxv,
                "nextclade_ds_tag": theiaviral.nextclade_ds_tag,
                "nextclade_aa_subs_mpxv": theiaviral.nextclade_aa_subs_mpxv,
                "nextclade_aa_dels_mpxv": theiaviral.nextclade_aa_dels_mpxv,
                "nextclade_clade_mpxv": theiaviral.nextclade_clade_mpxv,
                "nextclade_lineage_mpxv": theiaviral.nextclade_lineage_mpxv,
                "nextclade_qc_mpxv": theiaviral.nextclade_qc_mpxv,
                "nextclade_json_flu_ha": theiaviral.nextclade_json_flu_ha,
                "auspice_json_flu_ha": theiaviral.auspice_json_flu_ha,
                "nextclade_tsv_flu_ha": theiaviral.nextclade_tsv_flu_ha,
                "nextclade_ds_tag_flu_ha": theiaviral.nextclade_ds_tag_flu_ha,
                "nextclade_aa_subs_flu_ha": theiaviral.nextclade_aa_subs_flu_ha,
                "nextclade_aa_dels_flu_ha": theiaviral.nextclade_aa_dels_flu_ha,
                "nextclade_clade_flu_ha": theiaviral.nextclade_clade_flu_ha,
                "nextclade_qc_flu_ha": theiaviral.nextclade_qc_flu_ha,
                "nextclade_json_flu_na": theiaviral.nextclade_json_flu_na,
                "auspice_json_flu_na": theiaviral.auspice_json_flu_na,
                "nextclade_tsv_flu_na": theiaviral.nextclade_tsv_flu_na,
                "nextclade_ds_tag_flu_na": theiaviral.nextclade_ds_tag_flu_na,
                "nextclade_aa_subs_flu_na": theiaviral.nextclade_aa_subs_flu_na,
                "nextclade_aa_dels_flu_na": theiaviral.nextclade_aa_dels_flu_na,
                "nextclade_clade_flu_na": theiaviral.nextclade_clade_flu_na,
                "nextclade_qc_flu_na": theiaviral.nextclade_qc_flu_na,
                "irma_version": theiaviral.irma_version,
                "irma_docker": theiaviral.irma_docker,
                "irma_type": theiaviral.irma_type,
                "irma_subtype": theiaviral.irma_subtype,
                "irma_subtype_notes": theiaviral.irma_subtype_notes,
                "genoflu_version": theiaviral.genoflu_version,
                "genoflu_genotype": theiaviral.genoflu_genotype,
                "genoflu_all_segments": theiaviral.genoflu_all_segments,
                "genoflu_output_tsv": theiaviral.genoflu_output_tsv,
                "abricate_flu_type": theiaviral.abricate_flu_type,
                "abricate_flu_subtype":  theiaviral.abricate_flu_subtype,
                "abricate_flu_results": theiaviral.abricate_flu_results,
                "abricate_flu_database":  theiaviral.abricate_flu_database,
                "abricate_flu_version": theiaviral.abricate_flu_version
              }
          }
        }
      }
    }  
  }
 
  output {
    # Number of assembled viruses
    Int assembled_viruses = length(select_all(theiaviral.assembly_consensus_fasta))
    # Taxon of organisms identified
    Array[String] identified_organisms = select_all(ncbi_identify.taxon_name)
    # All assembled FASTA files
    Array[File] assemblies = select_all(theiaviral.assembly_consensus_fasta)
    # Workflow Versioning
    String theiaviral_panel_version = version_capture.phb_version
    String theiaviral_panel_analysis_date = version_capture.date
    # Standalone Kraken2 outputs
    String kraken2_version = read_QC_trim.kraken_version
    String kraken2_database = read_QC_trim.kraken_database
    String kraken2_docker = read_QC_trim.kraken_docker
    File kraken2_report_raw = read_QC_trim.kraken_report
    File kraken2_report_clean = select_first([read_QC_trim.kraken_report_clean])
    File kraken2_classified_report = select_first([read_QC_trim.kraken_classified_report])
    Float? kraken_percent_human_raw = read_QC_trim.kraken_human
    Float? kraken_percent_human_clean = read_QC_trim.kraken_human_clean
    # fastq_scan
    Int? fastq_scan_num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? fastq_scan_num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? fastq_scan_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    Int? fastq_scan_num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? fastq_scan_num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? fastq_scan_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    String? fastq_scan_docker = read_QC_trim.fastq_scan_docker
    File? fastq_scan_raw1_json = read_QC_trim.fastq_scan_raw1_json
    File? fastq_scan_raw2_json = read_QC_trim.fastq_scan_raw2_json
    File? fastq_scan_clean1_json = read_QC_trim.fastq_scan_clean1_json
    File? fastq_scan_clean2_json = read_QC_trim.fastq_scan_clean2_json
    # fastqc
    Int? fastqc_num_reads_raw1 = read_QC_trim.fastqc_raw1
    Int? fastqc_num_reads_raw2 = read_QC_trim.fastqc_raw2
    String? fastqc_num_reads_raw_pairs = read_QC_trim.fastqc_raw_pairs
    Int? fastqc_num_reads_clean1 = read_QC_trim.fastqc_clean1
    Int? fastqc_num_reads_clean2 = read_QC_trim.fastqc_clean2
    String? fastqc_num_reads_clean_pairs = read_QC_trim.fastqc_clean_pairs
    String? fastqc_version = read_QC_trim.fastqc_version
    String? fastqc_docker = read_QC_trim.fastqc_docker
    File? fastqc_raw1_html = read_QC_trim.fastqc_raw1_html
    File? fastqc_raw2_html = read_QC_trim.fastqc_raw2_html
    File? fastqc_clean1_html = read_QC_trim.fastqc_clean1_html
    File? fastqc_clean2_html = read_QC_trim.fastqc_clean2_html
    # trimming versioning
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
    String? fastp_version = read_QC_trim.fastp_version
    File? fastp_html_report = read_QC_trim.fastp_html_report
    # host decontamination outputs
    File? dehost_wf_dehost_read1 = read_QC_trim.dehost_wf_dehost_read1
    File? dehost_wf_dehost_read2 = read_QC_trim.dehost_wf_dehost_read2
    String? dehost_wf_host_accession = read_QC_trim.dehost_wf_host_accession
    File? dehost_wf_host_mapped_bam = read_QC_trim.dehost_wf_host_mapped_bam
    File? dehost_wf_host_mapped_bai = read_QC_trim.dehost_wf_host_mapped_bai
    File? dehost_wf_host_fasta = read_QC_trim.dehost_wf_host_fasta
    File? dehost_wf_host_mapping_stats = read_QC_trim.dehost_wf_host_mapping_stats
    File? dehost_wf_host_mapping_cov_hist = read_QC_trim.dehost_wf_host_mapping_cov_hist
    File? dehost_wf_host_flagstat = read_QC_trim.dehost_wf_host_flagstat
    Float? dehost_wf_host_mapping_coverage = read_QC_trim.dehost_wf_host_mapping_coverage
    Float? dehost_wf_host_mapping_mean_depth = read_QC_trim.dehost_wf_host_mapping_mean_depth
    Float? dehost_wf_host_percent_mapped_reads = read_QC_trim.dehost_wf_host_percent_mapped_reads
    File? dehost_wf_host_mapping_metrics = read_QC_trim.dehost_wf_host_mapping_metrics
    # NCBI scrubber
    File? read1_dehosted = read_QC_trim.read1_dehosted
    File? read2_dehosted = read_QC_trim.read2_dehosted
    Int? ncbi_scrub_human_spots_removed = read_QC_trim.ncbi_scrub_human_spots_removed
    String? ncbi_scrub_docker = read_QC_trim.ncbi_scrub_docker
    # bbduk
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    String bbduk_docker = read_QC_trim.bbduk_docker
  }
}