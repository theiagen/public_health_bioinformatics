version 1.0

import "../tasks/species_typing/task_cauris_cladetyper.wdl" as cauris_cladetyper
# import "../tasks/phylogenetic_inference/task_snippy.wdl" as snippy
# import "../tasks/assembly/task_mycosnp_consensus_assembly.wdl" as mycosnp
  # ADD SNIPPY TASK
  # Maybe add Mycosnp task?

workflow merlin_magic {
  meta {
    description: "Workflow for fungal species typing; based on the Bactopia subworkflow Merlin (https://bactopia.github.io/bactopia-tools/merlin/)"
  }
  input {
    String samplename
    String merlin_tag
    File assembly
    # File read1
    # File? read2
    # Boolean paired_end = true
  }
  if (merlin_tag == "Candida") {
    call cauris_cladetyper.cauris_cladetyper as cladetyper {
      input: 
        assembly_fasta = assembly,
        samplename = samplename
    }
    }
  output {
  # Candida Typing
  String? clade_type = cladetyper.gambit_cladetype
  String? clade_type_analysis_date = cladetyper.date
  String? clade_typer_version = cladetyper.version
  String? cladetyper_docker_image = cladetyper.gambit_cladetyper_docker_image
  String? cladetype_annotated_ref = cladetyper.clade_spec_ref
#   String snippy_version = read_string("VERSION")
#   File? snippy_aligned_fasta = snippy_task.snippy_aligned_fasta
#   File? snippy_bam = snippy_task.snippy_bam
#   File? snippy_bai = snippy_task.snippy_bai
#   File? snippy_bed = snippy_task.snippy_bed
#   File? snippy_consensus_fasta = snippy_task.snippy_consensus_fasta
#   File? snippy_subs_fasta = snippy_task.snippy_subs_fasta
#   File? snippy_csv = snippy_task.snippy_csv
#   File? snippy_filtered_vcf = snippy_task.snippy_filtered_vcf
#   File? snippy_gff = snippy_task.snippy_gff
#   File? snippy_html_report = snippy_task.snippy_html_report
#   File? snippy_log = snippy_task.snippy_log
#   File? snippy_raw_vcf = snippy_task.snippy_raw_vcf
#   File? snippy_subs_vcf = snippy_task.snippy_subs_vcf
#   File? snippy_tsv = snippy_task.snippy_tsv
#   File? snippy_txt = snippy_task.snippy_txt
#   File? snippy_vcf = snippy_task.snippy_vcf
#   File? snippy_vcf_gz = snippy_task.snippy_vcf_gz
#   File? snippy_vcf_gz_csi = snippy_task.snippy_vcf_gz_csi
 }
}