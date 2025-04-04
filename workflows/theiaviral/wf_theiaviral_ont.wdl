version 1.0
import "../../tasks/quality_control/read_filtering/task_porechop.wdl" as porechop_task
import "../../tasks/quality_control/read_filtering/task_nanoq.wdl" as nanoq_task
import "../../tasks/taxon_id/contamination/task_metabuli.wdl" as metabuli_task
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task
import "../../tasks/assembly/task_raven.wdl" as raven_task
import "../../tasks/quality_control/advanced_metrics/task_checkv.wdl" as checkv_task
import "../../tasks/phylogenetic_inference/utilities/task_skani.wdl" as skani_task
import "../../tasks/utilities/data_import/task_ncbi_datasets.wdl" as ncbi_datasets_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/gene_typing/variant_detection/task_clair3_variants.wdl" as clair3_task
import "../../tasks/utilities/data_handling/task_fasta_utilities.wdl" as fasta_utilities_task
import "../../tasks/assembly/task_bcftools_consensus.wdl" as bcftools_consensus_task
import "../../tasks/task_versioning.wdl" as versioning

workflow theiaviral_ont {
  meta {
    description: "..."
  }
  input {
    File read1
    String taxon_of_interest
    String samplename
    Boolean trim_adapters = false
    Boolean call_rasusa = false
    String? genome_length # required for RASUSA, could be optional otherwise # delete later
    Float downsampling_coverage = 150
  }
  if (trim_adapters) {
    call porechop_task.porechop as porechop {
      input:
        read1 = nanoq.filtered_read1,
        samplename = samplename
    }
  }
  call nanoq_task.nanoq as nanoq {
    input:
      read1 = read1,
      samplename = samplename
  }
  call metabuli_task.metabuli as metabuli {
    input:
      read1 = select_first([porechop.trimmed_reads, nanoq.filtered_read1]),
      samplename = samplename,
      taxon_of_interest = taxon_of_interest
  }
  if (call_rasusa) {
    if (! defined(genome_length)) {
      call ncbi_datasets_task.ncbi_datasets_viral_taxon_summary as ncbi_taxon_summary {
        input:
          taxon_id = taxon_of_interest
    }
  }
    call rasusa_task.rasusa as rasusa {
      input:
        read1 = metabuli.metabuli_read1_extract,
        samplename = samplename,
        coverage = downsampling_coverage,
        genome_length = select_first([genome_length, ncbi_taxon_summary.ncbi_datasets_avg_genome_length])
    }
  }
  call raven_task.raven as raven {
    input:
      read1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
      samplename = samplename
  }
  call checkv_task.checkv as checkv_denovo {
    input:
      assembly = raven.assembly_fasta,
      samplename = samplename
  }
  call skani_task.skani as skani {
    input:
      assembly_fasta = raven.assembly_fasta,
      samplename = samplename
  }
  # download the best reference determined from skani
  call ncbi_datasets_task.ncbi_datasets_download_genome_accession as ncbi_datasets {
    input:
      ncbi_accession = skani.skani_top_ani_accession,
      use_ncbi_virus = true
  }
  call consensus_qc_task.consensus_qc as denovo_qc {
    input:
      assembly_fasta = raven.assembly_fasta,
      reference_genome = ncbi_datasets.ncbi_datasets_assembly_fasta
  }
  call minimap2_task.minimap2 as minimap2 {
    input:
      query1 = read1,
      reference = ncbi_datasets.ncbi_datasets_assembly_fasta,
      samplename = samplename,
      mode = "map-ont",
      output_sam = true,
      long_read_flags = true
  }
  call parse_mapping_task.sam_to_sorted_bam as parse_mapping {
    input:
      sam = minimap2.minimap2_out,
      samplename = samplename
  }
  # Index the reference genome for Clair3
  call fasta_utilities_task.samtools_faidx as fasta_utilities{
    input:
      fasta = ncbi_datasets.ncbi_datasets_assembly_fasta,
  }
  call clair3_task.clair3_variants as clair3 {
    input:
      alignment_bam_file = parse_mapping.bam,
      alignment_bam_file_index = parse_mapping.bai,
      reference_genome_file = ncbi_datasets.ncbi_datasets_assembly_fasta,
      reference_genome_file_index = fasta_utilities.fai,
      sequencing_platform = "ont",
      samplename = samplename
  }
  call bcftools_consensus_task.bcftools_consensus as bcftools_consensus {
    input:
      reference_fasta = ncbi_datasets.ncbi_datasets_assembly_fasta,
      input_vcf = clair3.clair3_variants_vcf,
      samplename = samplename
  }
  call checkv_task.checkv as checkv_consensus {
    input:
      assembly = bcftools_consensus.bcftools_consensus_fasta,
      samplename = samplename
  }
  call consensus_qc_task.consensus_qc as consensus_qc {
    input:
      assembly_fasta = bcftools_consensus.bcftools_consensus_fasta,
      reference_genome = ncbi_datasets.ncbi_datasets_assembly_fasta
  }
  call versioning.version_capture {
    input:
  }
  output {
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
    String metabuli_database = metabuli.metabuli_database
    String metabuli_version = metabuli.metabuli_version
    String metabuli_docker = metabuli.metabuli_docker
    # ncbi datasets - taxon summary
    File? ncbi_taxon_summary_tsv = ncbi_taxon_summary.ncbi_datasets_taxon_summary_tsv
    String? ncbi_taxon_summary_avg_genome_length = ncbi_taxon_summary.ncbi_datasets_avg_genome_length
    String? ncbi_taxon_summary_version = ncbi_taxon_summary.ncbi_datasets_version
    String? ncbi_taxon_summary_docker = ncbi_taxon_summary.ncbi_datasets_docker
    # rasusa outputs - downsampled reads
    File? rasusa_read1_subsampled = rasusa.read1_subsampled
    File? rasusa_read2_subsampled = rasusa.read2_subsampled
    String? rasusa_version = rasusa.rasusa_version
    # raven outputs - denovo genome assembly
    File raven_denovo_assembly = raven.assembly_fasta
    String raven_denovo_version = raven.raven_version
    String raven_denovo_docker = raven.raven_docker
    # checkv denovo outputs - quality control
    File checkv_denovo_summary = checkv_denovo.checkv_summary
    String checkv_denovo_version = checkv_denovo.checkv_version
    # skani outputs - ANI-based reference genome selection
    File skani_report = skani.skani_report
    String skani_top_ani_accession = skani.skani_top_ani_accession
    String skani_database = skani.skani_database
    String skani_version = skani.skani_version
    String skani_docker = skani.skani_docker
    # ncbi_datasets outputs - download reference genome
    File skani_top_ani_fasta = ncbi_datasets.ncbi_datasets_assembly_fasta
    String ncbi_datasets_version = ncbi_datasets.ncbi_datasets_version
    String ncbi_datasets_docker = ncbi_datasets.ncbi_datasets_docker
    # denovo assembly statistics
    Int denovo_qc_number_N = denovo_qc.number_N
    Int denovo_qc_assembly_length_unambiguous = denovo_qc.number_ATCG
    Int denovo_qc_number_Degenerate = denovo_qc.number_Degenerate
    Int denovo_qc_number_Total = denovo_qc.number_Total
    Float denovo_qc_percent_reference_coverage = denovo_qc.percent_reference_coverage
    # minimap2 outputs - reads aligned to best reference
    File minimap2_out = minimap2.minimap2_out
    String minimap2_version = minimap2.minimap2_version
    String minimap2_docker = minimap2.minimap2_docker
    # parse_mapping outputs - sam to sorted bam conversion
    File assembly_to_ref_bam = parse_mapping.bam
    File assembly_to_ref_bai = parse_mapping.bai
    String parse_mapping_samtools_version = parse_mapping.samtools_version
    String parse_mapping_samtools_docker = parse_mapping.samtools_docker
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
    # bcftools_consensus outputs - consensus genome
    File bcftools_consensus_fasta = bcftools_consensus.bcftools_consensus_fasta
    File bcftools_norm_vcf = bcftools_consensus.bcftools_norm_vcf
    String bcftools_version = bcftools_consensus.bcftools_version
    String bcftools_docker = bcftools_consensus.bcftools_docker
    # checkv consensus outputs - quality control
    File checkv_consensus_summary = checkv_consensus.checkv_summary
    String checkv_consensus_version = checkv_consensus.checkv_version
    # consensus assembly statistics
    Int consensus_qc_number_N = consensus_qc.number_N
    Int consensus_qc_assembly_length_unambiguous = consensus_qc.number_ATCG
    Int consensus_qc_number_Degenerate = consensus_qc.number_Degenerate
    Int consensus_qc_number_Total = consensus_qc.number_Total
    Float consensus_qc_percent_reference_coverage = consensus_qc.percent_reference_coverage
    # versioning outputs
    String theiaviral_ont_version = version_capture.phb_version
    String theiaviral_ont_date = version_capture.date
  }
}