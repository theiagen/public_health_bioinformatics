version 1.0
import "../../tasks/quality_control/comparisons/task_screen.wdl" as screen
import "../../tasks/quality_control/read_filtering/task_porechop.wdl" as porechop_task
import "../../tasks/quality_control/read_filtering/task_nanoq.wdl" as nanoq_task
import "../../tasks/taxon_id/contamination/task_metabuli.wdl" as metabuli_task
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task
import "../../tasks/quality_control/basic_statistics/task_nanoplot.wdl" as nanoplot_task
import "../../tasks/assembly/task_raven.wdl" as raven_task
import "../../tasks/quality_control/advanced_metrics/task_checkv.wdl" as checkv_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/taxon_id/task_skani.wdl" as skani_task
import "../../tasks/utilities/data_import/task_ncbi_datasets.wdl" as ncbi_datasets_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics_task
import "../../tasks/utilities/data_handling/task_fasta_utilities.wdl" as fasta_utilities_task
import "../../tasks/gene_typing/variant_detection/task_clair3_variants.wdl" as clair3_task
import "../../tasks/assembly/task_bcftools_consensus.wdl" as bcftools_consensus_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/task_versioning.wdl" as versioning

workflow theiaviral_ont {
  meta {
    description: "De novo assembly, dynamic reference selection, and reference-based consensus calling for viral metagenomic/genomic data generated on ONT NGS platforms."
  }
  input {
    File read1
    String taxon_id
    String samplename
    Boolean trim_adapters = false
    Boolean call_rasusa = false
    String genome_length = "0" # required for RASUSA, could be optional otherwise # delete later
    File? reference_fasta # optional, if provided, will be used instead of dynamic reference selection
  }
  # raw read quality check, genome_length is not required for the rest of the workflow, so estimated coverage is not exposed to the outputs
  # incorporating nanoplot's estimated coverage output without requiring genome_length would require a conditional within the task to skip coverage if no genome length is provided
  call nanoplot_task.nanoplot as nanoplot_raw {
    input:
      read1 = read1,
      samplename = samplename,
      est_genome_length = select_first([genome_length, 1])
  }
  if (trim_adapters) {
    call porechop_task.porechop as porechop {
      input:
        read1 = read1,
        samplename = samplename
    }
  }
  # read filtering
  call nanoq_task.nanoq as nanoq {
    input:
      read1 = select_first([porechop.trimmed_reads, read1]),
      samplename = samplename
  }
  # taxonomic classification and read extraction
  call metabuli_task.metabuli as metabuli {
    input:
      read1 = select_first([porechop.trimmed_reads, nanoq.filtered_read1, read1]),
      samplename = samplename,
      taxon_id = taxon_id
  }
  if (call_rasusa) {
    # rasusa downsampling reads to specified coverage level
    call rasusa_task.rasusa as rasusa {
      input:
        read1 = metabuli.metabuli_read1_extract,
        samplename = samplename,
        genome_length = genome_length
    }
  }
  # clean read quality control
  call nanoplot_task.nanoplot as nanoplot_clean {
    input:
      read1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
      samplename = samplename,
      est_genome_length = select_first([genome_length, 1])
  }
  if (! defined(reference_fasta)) {
    # de novo assembly with raven
    call raven_task.raven as raven {
      input:
        read1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
        samplename = samplename
    }
    # quality control metrics for de novo assembly (ie. completeness, viral gene count, contamination)
    call checkv_task.checkv as checkv_denovo {
      input:
        assembly = raven.assembly_fasta,
        samplename = samplename
    }
    # quality control metrics for de novo assembly (ie. contigs, n50, GC content, genome length)
    call quast_task.quast as quast_denovo {
      input:
        assembly = raven.assembly_fasta,
        samplename = samplename
    }
    # ANI-based reference genome selection
    call skani_task.skani as skani {
      input:
        assembly_fasta = raven.assembly_fasta,
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
  call minimap2_task.minimap2 as minimap2 {
    input:
      query1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
      reference = select_first([ncbi_datasets.ncbi_datasets_assembly_fasta, reference_fasta]),
      samplename = samplename,
      mode = "map-ont",
      output_sam = true,
      long_read_flags = true
  }
  # generate bam file from sam output
  call parse_mapping_task.sam_to_sorted_bam as parse_mapping {
    input:
      sam = minimap2.minimap2_out,
      samplename = samplename
  }
  # quality control metrics for reads mapping to reference (ie. coverage, depth, base/map quality)
  call assembly_metrics_task.stats_n_coverage as assembly_metrics {
    input:
      bamfile = parse_mapping.bam,
      samplename = samplename
  }
  # Index the reference genome for Clair3
  call fasta_utilities_task.samtools_faidx as fasta_utilities{
    input:
      fasta = select_first([ncbi_datasets.ncbi_datasets_assembly_fasta, reference_fasta])
  }
  # variant calling with Clair3
  call clair3_task.clair3_variants as clair3 {
    input:
      alignment_bam_file = parse_mapping.bam,
      alignment_bam_file_index = parse_mapping.bai,
      reference_genome_file = select_first([ncbi_datasets.ncbi_datasets_assembly_fasta, reference_fasta]),
      reference_genome_file_index = fasta_utilities.fai,
      sequencing_platform = "ont",
      samplename = samplename
  }
  call parse_mapping_task.mask_low_coverage {
    input:
      bam = parse_mapping.bam,
      bai = parse_mapping.bai,
      reference_fasta = select_first([ncbi_datasets.ncbi_datasets_assembly_fasta, reference_fasta]),
  }
  # create consensus genome based on variant calls
  call bcftools_consensus_task.bcftools_consensus as bcftools_consensus {
    input:
      reference_fasta = mask_low_coverage.mask_reference_fasta,
      input_vcf = clair3.clair3_variants_vcf,
      samplename = samplename
  }
  # quality control metrics for consensus (ie. number of bases, degenerate bases, genome length)
  call consensus_qc_task.consensus_qc as consensus_qc {
    input:
      assembly_fasta = bcftools_consensus.bcftools_consensus_fasta,
      reference_genome = ncbi_datasets.ncbi_datasets_assembly_fasta
  }
  # quality control metrics for consensus (ie. completeness, viral gene count, contamination)
  call checkv_task.checkv as checkv_consensus {
    input:
      assembly = bcftools_consensus.bcftools_consensus_fasta,
      samplename = samplename
  }
  # quality control metrics for consensus (ie. contigs, n50, GC content, genome length)
  call quast_task.quast as quast_consensus {
    input:
      assembly = bcftools_consensus.bcftools_consensus_fasta,
      samplename = samplename
  }
  call versioning.version_capture {
    input:
  }
  output {
    # raw read quality control
    File nanoplot_html_raw = nanoplot_raw.nanoplot_html
    File nanoplot_tsv_raw = nanoplot_raw.nanoplot_tsv
    Int nanoplot_num_reads_raw1 = nanoplot_raw.num_reads
    Float nanoplot_r1_median_readlength_raw = nanoplot_raw.median_readlength
    Float nanoplot_r1_mean_readlength_raw = nanoplot_raw.mean_readlength
    Float nanoplot_r1_stdev_readlength_raw = nanoplot_raw.stdev_readlength
    Float nanoplot_r1_n50_raw = nanoplot_raw.n50
    Float nanoplot_r1_mean_q_raw = nanoplot_raw.mean_q
    Float nanoplot_r1_median_q_raw = nanoplot_raw.median_q
    # porechop outputs - adapter trimming
    File? porechop_trimmed_read1 = porechop.trimmed_reads
    String? porechop_version = porechop.porechop_version
    # nanoq outputs - read filtering
    File nanoq_filtered_read1 = nanoq.filtered_read1
    String nanoq_version = nanoq.version
    # metabuli outputs - taxonomic classification and read extraction
    File metabuli_report = metabuli.metabuli_report
    File metabuli_classified = metabuli.metabuli_classified
    File metabuli_read1_extract = metabuli.metabuli_read1_extract
    File metabuli_krona_report = metabuli.metabuli_krona_report
    String metabuli_database = metabuli.metabuli_database
    String metabuli_version = metabuli.metabuli_version
    String metabuli_docker = metabuli.metabuli_docker
    # rasusa outputs - downsampled reads
    File? rasusa_read1_subsampled = rasusa.read1_subsampled
    File? rasusa_read2_subsampled = rasusa.read2_subsampled
    String? rasusa_version = rasusa.rasusa_version
    # clean read quality control
    File nanoplot_html_clean = nanoplot_clean.nanoplot_html
    File nanoplot_tsv_clean = nanoplot_clean.nanoplot_tsv
    Int nanoplot_num_reads_clean1 = nanoplot_clean.num_reads
    Float nanoplot_r1_median_readlength_clean = nanoplot_clean.median_readlength
    Float nanoplot_r1_mean_readlength_clean = nanoplot_clean.mean_readlength
    Float nanoplot_r1_stdev_readlength_clean = nanoplot_clean.stdev_readlength
    Float nanoplot_r1_n50_clean = nanoplot_clean.n50
    Float nanoplot_r1_mean_q_clean = nanoplot_clean.mean_q
    Float nanoplot_r1_median_q_clean = nanoplot_clean.median_q
    # raven outputs - denovo genome assembly
    File? raven_denovo_assembly = raven.assembly_fasta
    String? raven_denovo_version = raven.raven_version
    String? raven_denovo_docker = raven.raven_docker
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
    # minimap2 outputs - reads aligned to best reference
    File minimap2_out = minimap2.minimap2_out
    String minimap2_version = minimap2.minimap2_version
    String minimap2_docker = minimap2.minimap2_docker
    # parse_mapping outputs - sam to sorted bam conversion
    File assembly_to_ref_bam = parse_mapping.bam
    File assembly_to_ref_bai = parse_mapping.bai
    String parse_mapping_samtools_version = parse_mapping.samtools_version
    String parse_mapping_samtools_docker = parse_mapping.samtools_docker
    # assembly_metrics outputs - read mapping quality control
    File assembly_metrics_report = assembly_metrics.metrics_txt
    File assembly_metrics_stats = assembly_metrics.stats
    File assembly_metrics_cov_hist = assembly_metrics.cov_hist
    File assembly_metrics_cov_stats = assembly_metrics.cov_stats
    File assembly_metrics_flagstat = assembly_metrics.flagstat
    Float assembly_metrics_coverage = assembly_metrics.coverage
    Float assembly_metrics_depth = assembly_metrics.depth
    Float assembly_metrics_meanbaseq = assembly_metrics.meanbaseq
    Float assembly_metrics_meanmapq = assembly_metrics.meanmapq
    Float assembly_metrics_percentage_mapped_reads = assembly_metrics.percentage_mapped_reads
    String assembly_metrics_date = assembly_metrics.date
    String assembly_metrics_samtools_version = assembly_metrics.samtools_version
    # fasta_utilities outputs - samtools faidx reference genome
    File fasta_utilities_fai = fasta_utilities.fai
    String fasta_utilities_samtools_version = fasta_utilities.samtools_version
    String fasta_utilities_samtools_docker = fasta_utilities.samtools_docker
    # clair3 outputs - variant calling
    File clair3_vcf = clair3.clair3_variants_vcf
    File? clair3_gvcf = clair3.clair3_variants_gvcf
    String clair3_model = clair3.clair3_model_used
    String clair3_version = clair3.clair3_version
    String clair3_docker = clair3.clair3_variants_docker_image
    # coverage_mask outputs - low coverage regions
    File mask_low_coverage_bed = mask_low_coverage.low_coverage_regions_bed
    File mask_low_coverage_all_coverage_bed = mask_low_coverage.all_coverage_regions_bed
    File mask_low_coverage_reference_fasta = mask_low_coverage.mask_reference_fasta
    String mask_low_coverage_bedtools_version = mask_low_coverage.bedtools_version
    String mask_low_coverage_bedtools_docker = mask_low_coverage.bedtools_docker
    # bcftools_consensus outputs - consensus genome
    File bcftools_consensus_fasta = bcftools_consensus.bcftools_consensus_fasta
    File bcftools_norm_vcf = bcftools_consensus.bcftools_norm_vcf
    String bcftools_version = bcftools_consensus.bcftools_version
    String bcftools_docker = bcftools_consensus.bcftools_docker
    # consensus assembly statistics
    Int consensus_qc_number_N = consensus_qc.number_N
    Int consensus_qc_assembly_length_unambiguous = consensus_qc.number_ATCG
    Int consensus_qc_number_Degenerate = consensus_qc.number_Degenerate
    Int consensus_qc_number_Total = consensus_qc.number_Total
    Float consensus_qc_percent_reference_coverage = consensus_qc.percent_reference_coverage
    # checkv_consensus outputs - consensus assembly quality control
    File checkv_consensus_summary = checkv_consensus.checkv_summary
    File checkv_consensus_contamination = checkv_consensus.checkv_contamination
    Float checkv_consensus_total_contamination = checkv_consensus.total_contamination
    Float checkv_consensus_total_completeness = checkv_consensus.total_completeness
    Int checkv_consensus_total_genes = checkv_consensus.total_genes
    String checkv_consensus_version = checkv_consensus.checkv_version
    # quast_consensus outputs - consensus assembly quality control
    File quast_consensus_report = quast_consensus.quast_report
    Int quast_consensus_genome_length = quast_consensus.genome_length
    Int quast_consensus_number_contigs = quast_consensus.number_contigs
    Int quast_consensus_n50_value = quast_consensus.n50_value
    Int quast_consensus_largest_contig = quast_consensus.largest_contig
    Float quast_consensus_gc_percent = quast_consensus.gc_percent
    Float quast_consensus_uncalled_bases = quast_consensus.uncalled_bases
    String quast_consensus_pipeline_date = quast_consensus.pipeline_date
    String quast_consensus_version = quast_consensus.version
    String quast_consensus_docker = quast_consensus.quast_docker
    # versioning outputs
    String theiaviral_ont_version = version_capture.phb_version
    String theiaviral_ont_date = version_capture.date
  }
}