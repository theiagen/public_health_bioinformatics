version 1.0

import "../../tasks/assembly/task_metaspades.wdl" as metaspades_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/quality_control/task_pilon.wdl" as pilon_task

workflow metaspades_assembly_pe {
  meta {
    description: "De novo assembly using metaspades and correction with pilon"
  }
    input {
    String samplename
    File read1
    File read2
    String? kmers
    String? metaspades_opts
  }
  call metaspades_task.metaspades_pe as metaspades {
    input:
      read1_cleaned = read1,
      read2_cleaned = read2,
      samplename = samplename,
      kmers = kmers,
      metaspades_opts = metaspades_opts
  }

  call minimap2_task.minimap2 {
    input:
      query1 = read1,
      query2 = read2, 
      reference = metaspades.assembly_fasta,
      samplename = samplename,
      mode = "sr",
      output_sam = true
  }
  call parse_mapping_task.sam_to_sorted_bam {
    input:
      sam = minimap2.minimap2_out,
      samplename = samplename
  }
  call pilon_task.pilon {
    input:
      assembly = metaspades.assembly_fasta,
      bam = sam_to_sorted_bam.bam,
      bai = sam_to_sorted_bam.bai,
      samplename = samplename
  }
  output {
    # metaspades output
    String metaspades_version = metaspades.metaspades_version
    String metaspades_docker = metaspades.metaspades_docker
    # minimap2 output
    String minimap2_version = minimap2.minimap2_version
    String minimap2_docker = minimap2.minimap2_docker
    # samtools output
    String samtools_version = sam_to_sorted_bam.samtools_version
    String samtools_docker = sam_to_sorted_bam.samtools_docker
    # pilon output
    File assembly_fasta = pilon.assembly_fasta
    String pilon_version = pilon.pilon_version
    String pilon_docker = pilon.pilon_docker
  }

}