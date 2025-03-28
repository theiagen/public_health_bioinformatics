version 1.0
import "../../tasks/quality_control/read_filtering/task_nanoq.wdl" as nanoq_task
import "../../tasks/taxon_id/contamination/task_metabuli.wdl" as metabuli_task
import "../../tasks/assembly/task_flye.wdl" as flye_task
import "../../tasks/assembly/task_raven.wdl" as raven_task
import "../../tasks/phylogenetic_inference/utilities/task_skani.wdl" as skani_task
import "../../tasks/utilities/data_import/task_ncbi_datasets.wdl" as ncbi_datasets_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/assembly/task_ivar_consensus.wdl" as ivar_consensus_task
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/quality_control/advanced_metrics/task_checkv.wdl" as checkv_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task

workflow theiaviral_ont{
  meta {
    description: "..."
  }
  input {
    File read1
    String taxon_of_interest
    String samplename
    File metabuli_db # delete this later only used for miniwdl testing
    File taxonomy_path # delete this later only used for miniwdl testing
    File skani_db # delete this later only used for miniwdl testing
    File checkv_db # delete this later only used for miniwdl testing

    # rasusa downsampling inputs
    Float downsampling_coverage = 0
    String? rasusa_bases # delete after miniwdl testing
    Int? rasusa_seed # delete after
    Float? rasusa_fraction_of_reads
    Int? rasusa_number_of_reads
    String? genome_length # required for RASUSA, could be optional otherwise

    # assembler inputs
    Boolean run_raven = false
  }
  call nanoq_task.nanoq as nanoq {
    input:
      read1 = read1,
      samplename = samplename
  }
  call metabuli_task.metabuli as metabuli {
    input:
      read1 = nanoq.filtered_read1,
      samplename = samplename,
      taxon_of_interest = taxon_of_interest,
      metabuli_db = metabuli_db,
      taxonomy_path = taxonomy_path
  }

  if (downsampling_coverage > 0.0) {
    call rasusa_task.rasusa {
      input:
        read1 = metabuli.metabuli_read1_extract,
        samplename = samplename,
        coverage = downsampling_coverage,
        genome_length = select_first([genome_length, ""]),
        num = rasusa_number_of_reads,
        frac = rasusa_fraction_of_reads,
        seed = rasusa_seed,
        bases = rasusa_bases
    }  
  }
  if (run_raven) {
    call raven_task.raven as raven {
      input:
        read1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
        samplename = samplename
    }
  } 
  if (!run_raven) {
    call flye_task.flye as flye {
      input:
        read1 = select_first([rasusa.read1_subsampled, metabuli.metabuli_read1_extract]),
        samplename = samplename,
        asm_coverage = 50,
        genome_length = select_first([genome_length, ""])
    }
  }

  call skani_task.skani as skani {
    input:
      assembly_fasta = select_first([raven.assembly_fasta, flye.assembly_fasta]),
      samplename = samplename,
      skani_db = skani_db
  }
  # download the best reference determined from skani
  call ncbi_datasets_task.ncbi_datasets_download_genome_accession as ncbi_datasets {
    input:
      ncbi_accession = skani.skani_top_ani_accession,
      use_ncbi_virus = true
  }

  # NEED to make QC optional for de novo assembly
  call checkv_task.checkv as checkv_denovo {
    input:
      assembly = select_first([raven.assembly_fasta, flye.assembly_fasta]),
      samplename = samplename,
      checkv_db = checkv_db
  }
  call consensus_qc_task.consensus_qc as consensus_qc_denovo {
    input:
      assembly_fasta = select_first([raven.assembly_fasta, flye.assembly_fasta]),
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
  call ivar_consensus_task.consensus as ivar {
    input:
      bamfile = parse_mapping.bam,
      samplename = samplename,
      reference_genome = ncbi_datasets.ncbi_datasets_assembly_fasta
  }
  call checkv_task.checkv as checkv_consensus {
    input:
      assembly = ivar.consensus_seq,
      samplename = samplename,
      checkv_db = checkv_db
  }
  call consensus_qc_task.consensus_qc as consensus_qc_consensus {
    input:
      assembly_fasta = ivar.consensus_seq,
      reference_genome = ncbi_datasets.ncbi_datasets_assembly_fasta
  }

  call versioning.version_capture {
    input:
  }

  output {
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
    # raven outputs - denovo genome assembly
    File? raven_assembly = raven.assembly_fasta
    String? raven_version = raven.raven_version
    String? raven_docker = raven.raven_docker
    # flye outputs - denovo genome assembly
    File? flye_assembly = flye.assembly_fasta
    File? flye_assembly_graph = flye.assembly_graph_gfa
    File? flye_assembly_info = flye.assembly_info
    String? flye_version = flye.flye_version
    String? flye_docker = flye.flye_docker
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
    # minimap2 outputs - reads aligned to best reference
    File minimap2_out = minimap2.minimap2_out
    String minimap2_version = minimap2.minimap2_version
    String minimap2_docker = minimap2.minimap2_docker
    # parse_mapping outputs - sam to sorted bam conversion
    File assembly_to_ref_bam = parse_mapping.bam
    File assembly_to_ref_bai = parse_mapping.bai
    String parse_mapping_samtools_version = parse_mapping.samtools_version
    String parse_mapping_samtools_docker = parse_mapping.samtools_docker
    # ivar outputs - consensus genome
    File ivar_consensus_seq = ivar.consensus_seq
    String ivar_version = ivar.ivar_version
    String ivar_samtools_version = ivar.samtools_version
    # versioning outputs
    String theiaviral_ont_version = version_capture.phb_version
    String theiaviral_ont_date = version_capture.date
    # checkv outputs - quality control
    File? checkv_denovo_summary = checkv_denovo.checkv_summary
    String? checkv_denovo_version = checkv_denovo.checkv_version
    File checkv_consensus_summary = checkv_consensus.checkv_summary
    String checkv_consensus_version = checkv_consensus.checkv_version
    # basic assembly statistics
    Int? denovo_number_N = consensus_qc_denovo.number_N
    Int? denovo_assembly_length_unambiguous = consensus_qc_denovo.number_ATCG
    Int? denovo_number_Degenerate = consensus_qc_denovo.number_Degenerate
    Int? denovo_number_Total = consensus_qc_denovo.number_Total
    Int? consensus_number_N = consensus_qc_consensus.number_N
    Int? consensus_assembly_length_unambiguous = consensus_qc_consensus.number_ATCG
    Int? consensus_number_Degenerate = consensus_qc_consensus.number_Degenerate
    Int? consensus_number_Total = consensus_qc_consensus.number_Total
    # rasusa outputs
    String? rasusa_version = rasusa.rasusa_version
  }
}