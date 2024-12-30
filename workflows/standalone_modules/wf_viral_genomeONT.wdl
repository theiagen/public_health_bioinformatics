version 1.0

import "../../tasks/quality_control/basic_statistics/task_fastqc.wdl" as fastqc_task
import "../../tasks/quality_control/read_filtering/task_seqtk_trim.wdl" as seqtk_task
import "../../tasks/utilities/file_handling/task_plot_read_lengths.wdl" as plot_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/file_handling/task_samtools_process.wdl" as samtools_task
import "../../tasks/utilities/file_handling/task_filter_coverage.wdl" as filter_task
import "../../tasks/polishing/task_medaka.wdl" as medaka_task

workflow refnaap_workflow {
  meta {
    description: "A generalizable workflow for reference-guided genome assembly and analysis of an ONT sequencing file."
    version: "1.0"
  }
  input {
    File fastq_file               # input FASTQ file
    File reference_file           # Reference genome file
    String sample_name            # Sample name for output file naming
    Int? left_trim           
    Int? right_trim
    Int? min_length
    Int min_coverage = 10         # Minimum coverage for filtering
  }
  # Step 1: Run FastQC
  call fastqc_task.fastqc_se {
    input:
      read1 = fastq_file
  }
  # Step 2: Trimming and Filtering
  call seqtk_task.seqtk_trim {
    input:
      fastq_file = fastq_file,
      left_trim = left_trim,
      right_trim = right_trim,
      min_length = min_length,
      sample_name = sample_name
  }
  # Step 3: Plot Read Length Distribution
  call plot_task.plot_read_length_distribution {
    input:
      filtered_fastq = seqtk_trim.trimmed_fastq,
      sample_name = sample_name
  }
  # Step 4: Align to Reference minimap2
  call minimap2_task.minimap2_align_ont {
    input:
      fastq = seqtk_trim.trimmed_fastq,
      reference = reference_file,
      sample_name = sample_name
  }
  # Step 5: Process SAM/BAM and Calculate Coverage and statistics
  call samtools_task.samtools_process {
    input:
      sam_file = minimap2_align_ont.aligned_sam,
      sample_name = sample_name
  }
  # Step 6: Filter Coverage and Extract Scaffolds
  call filter_task.filter_coverage {
    input:
      bam_file = samtools_process.sorted_bam_file,
      reference = reference_file,
      min_coverage = min_coverage,
      sample_name = sample_name
  }
  # Step 7: Generate Consensus with Medaka
  call medaka_task.medaka_consensus {
    input:
      unpolished_fasta = filter_coverage.scaffolds_fasta,
      samplename = sample_name,
      read1 = fastq_file
  }
  output {
    File fastqc_report = fastqc_se.read1_fastqc_html
    File trimmed_fastq_file = seqtk_trim.trimmed_fastq
    File read_length_plot = plot_read_length_distribution.read_length_plot
    File aligned_sam_file = minimap2_align_ont.aligned_sam
    File sorted_bam_file = samtools_process.sorted_bam_file
    File coverage_file = samtools_process.coverage_file
    File scaffolds_fasta = filter_coverage.scaffolds_fasta
    File medaka_consensus_fasta = medaka_consensus.medaka_fasta
    String medaka_version = medaka_consensus.medaka_version
    String medaka_model = medaka_consensus.medaka_model
  }
}
