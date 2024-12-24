version 1.0

import "../../tasks/gene_typing/variant_detection/task_clair3_variants.wdl" as clair3
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/utilities/data_handling/task_fasta_utilities.wdl" as fasta_utilities_task 
import "../../tasks/task_versioning.wdl" as versioning

# MVP Clair3 variant calling workflow for ONT data
workflow clair3_variants_ont {
  meta {
    description: "Call variants using Clair3 for ONT data"
  }
  input {
    File read1
    File reference_genome_file
    String samplename
    String model = "r941_prom_hac_g360+g422"
    String? docker
    Int? memory
    Int? cpu
    Int? disk_size
    Int variant_quality = 2
    Boolean enable_gvcf = false
  }
  # Call the minimap2 task with recommended options for ONT data, when we expand to PacBio need to refactor logic
  call minimap2_task.minimap2 {
    input:
      query1 = read1,
      reference = reference_genome_file,
      samplename = samplename,
      mode = "map-ont",
      output_sam = true,
      additional_options = "-L --cs --MD"
  }
  # Parse the minimap2 output to a sorted BAM file, and index it, expected by clair3
  call parse_mapping_task.sam_to_sorted_bam {
    input:
      sam = minimap2.minimap2_out,
      samplename = samplename
  }
  # Index the reference genome for Clair3
  call fasta_utilities_task.samtools_faidx {
    input:
      fasta = reference_genome_file
  }

  call clair3.clair3_variants {
    input:
      alignment_bam_file = sam_to_sorted_bam.bam,
      alignment_bam_file_index = sam_to_sorted_bam.bai,
      reference_genome_file = reference_genome_file,
      reference_genome_file_index = samtools_faidx.fai,
      sequencing_platform = "ont",
      samplename = samplename,
      model = model,
      docker = docker,
      memory = memory,
      cpu = cpu,
      disk_size = disk_size,
      qual = variant_quality,
      enable_gvcf = enable_gvcf
  }

  call versioning.version_capture {
    input:
  }
  
  output {
    # Data handling - Read Alignment- samtools version
    String samtools_version = sam_to_sorted_bam.samtools_version
    File aligned_bam = sam_to_sorted_bam.bam
    File aligned_bai = sam_to_sorted_bam.bai
    # Data handling - Reference genome - samtools version
    File aligned_fai = samtools_faidx.fai
    #Clair3 variant calling
    String clair3_variants_wf_version = version_capture.phb_version
    File clair3_variants_final_vcf = clair3_variants.clair3_variants_final_vcf
    File clair3_variants_pileup_vcf = clair3_variants.clair3_variants_pileup_vcf
    File clair3_variants_full_alignment_vcf = clair3_variants.clair3_variants_full_alignment_vcf
    File? clair3_variants_gvcf = clair3_variants.clair3_variants_gvcf
    String clair3_docker_image = clair3_variants.clair3_variants_docker_image
    String clair3_model = clair3_variants.clair3_model
  }

}