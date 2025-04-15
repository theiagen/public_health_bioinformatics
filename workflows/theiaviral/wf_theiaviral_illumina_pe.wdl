version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc
import "../../tasks/quality_control/comparisons/task_screen.wdl" as read_screen_task
# assembly
import "../../tasks/assembly/task_spades.wdl" as spades_task
import "../../tasks/assembly/task_megahit.wdl" as megahit_task
import "../../tasks/quality_control/advanced_metrics/task_checkv.wdl" as checkv_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/taxon_id/task_skani.wdl" as skani_task
import "../../tasks/utilities/data_import/task_ncbi_datasets.wdl" as ncbi_datasets_task
# bwa
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics_task
import "../../tasks/utilities/data_handling/task_fasta_utilities.wdl" as fasta_utilities_task
# consensus module
import "../../tasks/assembly/task_bcftools_consensus.wdl" as bcftools_consensus_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/task_versioning.wdl" as versioning

workflow theiaviral_illumina_pe {
  meta {
    description: "De novo assembly, dynamic reference selection, and reference-based consensus calling for viral metagenomic/genomic data generated on Illumina paired-end NGS platforms."
  }
  input {
    File read1
    File read2
    String taxon_id
    String samplename
    Float min_mask_depth = 20 # minimum depth for masking low coverage regions
    Boolean skip_screen = true # if false, run clean read screening
    File? reference_fasta # optional, if provided, will be used instead of dynamic reference selection

    Boolean call_metaviralspades = false
    Boolean call_metaspades = false
    Boolean call_spades = false
    Boolean call_megahit = false
  }
  # get the PHB version
  call versioning.version_capture {
    input:
  }
  # read QC, classification, extraction, and trimming
  # NEED to expose theiaviral specific parameters, e.g. exclusion_extraction, extract_unclassified
  call read_qc.read_QC_trim_pe as read_QC_trim {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      taxon_id = taxon_id,
      workflow_series = "theiaviral"
  }
  # clean read screening
  if (! skip_screen) {
    call read_screen_task.check_reads as clean_check_reads {
      input:
        read1 = select_first([read_QC_trim.kraken2_extracted_read1]),
        read2 = select_first([read_QC_trim.kraken2_extracted_read2]),
        min_reads = 0,
        min_basepairs = 0,
        min_genome_length = 0,
        max_genome_length = 10000000000,
        min_coverage = 0,
        min_proportion = 0
    }
  }

  if (select_first([clean_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
    # run de novo if no reference genome is provided so we can select a reference
    if (! defined(reference_fasta)) {
      # de novo assembly benchmarking
      if (call_metaviralspades) {
        call spades_task.metaviralspades_pe {
          input:
            read1_cleaned = select_first([read_QC_trim.kraken2_extracted_read1]),
            read2_cleaned = select_first([read_QC_trim.kraken2_extracted_read2]),
            samplename = samplename
        }
      }
      if (call_metaspades) {
        call spades_task.metaspades_pe {
          input:
            read1_cleaned = select_first([read_QC_trim.kraken2_extracted_read1]),
            read2_cleaned = select_first([read_QC_trim.kraken2_extracted_read2]),
            samplename = samplename
        }
      }
      if (call_spades) {
        call spades_task.spades_pe {
          input:
            read1_cleaned = select_first([read_QC_trim.kraken2_extracted_read1]),
            read2_cleaned = select_first([read_QC_trim.kraken2_extracted_read2]),
            samplename = samplename
        }
      }
      if (call_megahit) {
        call megahit_task.megahit_pe {
          input:
            read1_cleaned = select_first([read_QC_trim.kraken2_extracted_read1]),
            read2_cleaned = select_first([read_QC_trim.kraken2_extracted_read2]),
            samplename = samplename
        }
      }
  
      # quality control metrics for de novo assembly (ie. completeness, viral gene count, contamination)
      call checkv_task.checkv as checkv_denovo {
        input:
          assembly = select_first([metaviralspades_pe.assembly_fasta, metaspades_pe.assembly_fasta, spades_pe.assembly_fasta, megahit_pe.assembly_fasta]),
          samplename = samplename
      }
      # quality control metrics for de novo assembly (ie. contigs, n50, GC content, genome length)
      call quast_task.quast as quast_denovo {
        input:
          assembly = select_first([metaviralspades_pe.assembly_fasta, metaspades_pe.assembly_fasta, spades_pe.assembly_fasta, megahit_pe.assembly_fasta]),
          samplename = samplename
      }
      # ANI-based reference genome selection
      call skani_task.skani as skani {
        input:
          assembly_fasta = select_first([metaviralspades_pe.assembly_fasta, metaspades_pe.assembly_fasta, spades_pe.assembly_fasta, megahit_pe.assembly_fasta]),
          samplename = samplename
      }
      # download the best reference determined from skani
      call ncbi_datasets_task.ncbi_datasets_download_genome_accession as ncbi_datasets {
        input:
          ncbi_accession = skani.skani_top_accession,
          use_ncbi_virus = true
      }
    }
    # align assembly to reference genome
  
    # generate bam file from sam output
  #  call parse_mapping_task.sam_to_sorted_bam as parse_mapping {
   #   input:
    #    sam = minimap2.minimap2_out,
     #   samplename = samplename
    #}
    # quality control metrics for reads mapping to reference (ie. coverage, depth, base/map quality)
  #  call assembly_metrics_task.stats_n_coverage as assembly_metrics {
   #   input:
    #    bamfile = parse_mapping.bam,
     #   samplename = samplename
  #  }
    # Index the reference genome for Clair3
   # call fasta_utilities_task.samtools_faidx as fasta_utilities{
    #  input:
     #   fasta = select_first([ncbi_datasets.ncbi_datasets_assembly_fasta, reference_fasta])
  #  }
    # variant calling
  
    # mask low coverage regions with Ns
  #  call parse_mapping_task.mask_low_coverage {
   #   input:
    #    bam = parse_mapping.bam,
     #   bai = parse_mapping.bai,
      #  reference_fasta = select_first([ncbi_datasets.ncbi_datasets_assembly_fasta, reference_fasta]),
       # min_depth = min_mask_depth
  #  }
    # create consensus genome based on variant calls
   # call bcftools_consensus_task.bcftools_consensus as bcftools_consensus {
    #  input:
     #   reference_fasta = mask_low_coverage.mask_reference_fasta,
      #  input_vcf = clair3.clair3_variants_vcf,
       # samplename = samplename
  #  }
    # quality control metrics for consensus (ie. number of bases, degenerate bases, genome length)
   # call consensus_qc_task.consensus_qc as consensus_qc {
    #  input:
     #   assembly_fasta = bcftools_consensus.bcftools_consensus_fasta,
      #  reference_genome = ncbi_datasets.ncbi_datasets_assembly_fasta
  #  }
    # quality control metrics for consensus (ie. completeness, viral gene count, contamination)
  #  call checkv_task.checkv as checkv_consensus {
   #   input:
    #    assembly = bcftools_consensus.bcftools_consensus_fasta,
     #   samplename = samplename
  #  }
    # quality control metrics for consensus (ie. contigs, n50, GC content, genome length)
   # call quast_task.quast as quast_consensus {
    #  input:
     #   assembly = bcftools_consensus.bcftools_consensus_fasta,
      #  samplename = samplename
#  }
  }
  output {
    # versioning outputs
    String? theiaviral_illumina_pe_version = version_capture.phb_version
    String? theiaviral_illumina_pe_date = version_capture.date
    # raw read quality control
    Int? fastq_scan_num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? fastq_scan_num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? fastq_scan_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    String? fastq_scan_docker = read_QC_trim.fastq_scan_docker
    File? fastq_scan_raw1_json = read_QC_trim.fastq_scan_raw1_json
    File? fastq_scan_raw2_json = read_QC_trim.fastq_scan_raw2_json
    # NCBI scrubbing
    File? read1_dehosted = read_QC_trim.read1_dehosted
    File? read2_dehosted = read_QC_trim.read2_dehosted
    Int? ncbi_scrub_human_spots_removed = read_QC_trim.ncbi_scrub_human_spots_removed
    String? ncbi_scrub_docker = read_QC_trim.ncbi_scrub_docker
    # trimming outputs - adapter trimming
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
    String? fastp_version = read_QC_trim.fastp_version
    File? fastp_html_report = read_QC_trim.fastp_html_report
    # bbduk outputs
    String? bbduk_docker = read_QC_trim.bbduk_docker
    File? bbduk_read1_clean = read_QC_trim.read1_clean
    File? bbduk_read2_clean = read_QC_trim.read2_clean
    # kraken2 outputs - taxonomic classification and read extraction
    String? kraken_version = read_QC_trim.kraken_version
    String? kraken_docker = read_QC_trim.kraken_docker
    String? kraken_database = read_QC_trim.kraken_database
    String? kraken_report = read_QC_trim.kraken_report
    File? kraken2_extracted_read1 = read_QC_trim.kraken2_extracted_read1
    File? kraken2_extracted_read2 = read_QC_trim.kraken2_extracted_read2
    # clean read quality control
    Int? fastq_scan_num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? fastq_scan_num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? fastq_scan_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    File? fastq_scan_clean1_json = read_QC_trim.fastq_scan_clean1_json
    File? fastq_scan_clean2_json = read_QC_trim.fastq_scan_clean2_json
    # denovo genome assembly
    File? assembly_fasta = select_first([spades_pe.assembly_fasta, metaspades_pe.assembly_fasta, metaviralspades_pe.assembly_fasta, megahit_pe.assembly_fasta])
    # checkv_denovo outputs - denovo assembly quality control
    File? checkv_denovo_summary = checkv_denovo.checkv_summary
    File? checkv_denovo_contamination = checkv_denovo.checkv_contamination
    Float? checkv_denovo_total_contamination = checkv_denovo.total_contamination
    Float? checkv_denovo_total_completeness = checkv_denovo.total_completeness
    Int? checkv_denovo_total_genes = checkv_denovo.total_genes
    String? checkv_denovo_version = checkv_denovo.checkv_version
    # quast_denovo outputs - denovo assembly quality control
    File? quast_denovo_report = quast_denovo.quast_report
    Int? quast_denovo_genome_length = quast_denovo.genome_length
    Int? quast_denovo_number_contigs = quast_denovo.number_contigs
    Int? quast_denovo_n50_value = quast_denovo.n50_value
    Int? quast_denovo_largest_contig = quast_denovo.largest_contig
    Float? quast_denovo_gc_percent = quast_denovo.gc_percent
    Float? quast_denovo_uncalled_bases = quast_denovo.uncalled_bases
    String? quast_denovo_pipeline_date = quast_denovo.pipeline_date
    String? quast_denovo_version = quast_denovo.version
    String? quast_denovo_docker = quast_denovo.quast_docker
    # skani outputs - ANI-based reference genome selection
    File? skani_report = skani.skani_report
    String? skani_top_accession = skani.skani_top_accession
    Float? skani_top_ani = skani.skani_top_ani
    Float? skani_top_ref_coverage = skani.skani_top_ref_coverage
    String? skani_database = skani.skani_database
    String? skani_version = skani.skani_version
    String? skani_docker = skani.skani_docker
    # ncbi_datasets outputs - download reference genome
    File? skani_top_ani_fasta = ncbi_datasets.ncbi_datasets_assembly_fasta
    String? ncbi_datasets_version = ncbi_datasets.ncbi_datasets_version
    String? ncbi_datasets_docker = ncbi_datasets.ncbi_datasets_docker
    # bwa outputs - reads aligned to best reference

    # parse_mapping outputs - sam to sorted bam conversion
  #  File assembly_to_ref_bam = parse_mapping.bam
   # File assembly_to_ref_bai = parse_mapping.bai
    #String parse_mapping_samtools_version = parse_mapping.samtools_version
#    String parse_mapping_samtools_docker = parse_mapping.samtools_docker
    # assembly_metrics outputs - read mapping quality control
 #   File assembly_metrics_report = assembly_metrics.metrics_txt
  #  File assembly_metrics_stats = assembly_metrics.stats
   # File assembly_metrics_cov_hist = assembly_metrics.cov_hist
    #File assembly_metrics_cov_stats = assembly_metrics.cov_stats
#    File assembly_metrics_flagstat = assembly_metrics.flagstat
 #   Float assembly_metrics_coverage = assembly_metrics.coverage
  #  Float assembly_metrics_depth = assembly_metrics.depth
   # Float assembly_metrics_meanbaseq = assembly_metrics.meanbaseq
    #Float assembly_metrics_meanmapq = assembly_metrics.meanmapq
#    Float assembly_metrics_percentage_mapped_reads = assembly_metrics.percentage_mapped_reads
 #   String assembly_metrics_date = assembly_metrics.date
  #  String assembly_metrics_samtools_version = assembly_metrics.samtools_version
    # fasta_utilities outputs - samtools faidx reference genome
   # File fasta_utilities_fai = fasta_utilities.fai
    #String fasta_utilities_samtools_version = fasta_utilities.samtools_version
#    String fasta_utilities_samtools_docker = fasta_utilities.samtools_docker
    # variant calling

    # coverage_mask outputs - low coverage regions
#    File mask_low_coverage_bed = mask_low_coverage.low_coverage_regions_bed
 #   File mask_low_coverage_all_coverage_bed = mask_low_coverage.all_coverage_regions_bed
  #  File mask_low_coverage_reference_fasta = mask_low_coverage.mask_reference_fasta
   # String mask_low_coverage_bedtools_version = mask_low_coverage.bedtools_version
    #String mask_low_coverage_bedtools_docker = mask_low_coverage.bedtools_docker
    # bcftools_consensus outputs - consensus genome
#    File bcftools_consensus_fasta = bcftools_consensus.bcftools_consensus_fasta
 #   File bcftools_norm_vcf = bcftools_consensus.bcftools_norm_vcf
  #  String bcftools_version = bcftools_consensus.bcftools_version
   # String bcftools_docker = bcftools_consensus.bcftools_docker
    # consensus assembly statistics
    #Int consensus_qc_number_N = consensus_qc.number_N
#    Int consensus_qc_assembly_length_unambiguous = consensus_qc.number_ATCG
 #   Int consensus_qc_number_Degenerate = consensus_qc.number_Degenerate
  #  Int consensus_qc_number_Total = consensus_qc.number_Total
   # Float consensus_qc_percent_reference_coverage = consensus_qc.percent_reference_coverage
    # checkv_consensus outputs - consensus assembly quality control
    #File checkv_consensus_summary = checkv_consensus.checkv_summary
#    File checkv_consensus_contamination = checkv_consensus.checkv_contamination
 #   Float checkv_consensus_total_contamination = checkv_consensus.total_contamination
  #  Float checkv_consensus_total_completeness = checkv_consensus.total_completeness
   # Int checkv_consensus_total_genes = checkv_consensus.total_genes
    #String checkv_consensus_version = checkv_consensus.checkv_version
    # quast_consensus outputs - consensus assembly quality control
#    File quast_consensus_report = quast_consensus.quast_report
 #   Int quast_consensus_genome_length = quast_consensus.genome_length
  #  Int quast_consensus_number_contigs = quast_consensus.number_contigs
   # Int quast_consensus_n50_value = quast_consensus.n50_value
    #Int quast_consensus_largest_contig = quast_consensus.largest_contig
#    Float quast_consensus_gc_percent = quast_consensus.gc_percent
 #   Float quast_consensus_uncalled_bases = quast_consensus.uncalled_bases
  #  String quast_consensus_pipeline_date = quast_consensus.pipeline_date
   # String quast_consensus_version = quast_consensus.version
    #String quast_consensus_docker = quast_consensus.quast_docker
  }
}