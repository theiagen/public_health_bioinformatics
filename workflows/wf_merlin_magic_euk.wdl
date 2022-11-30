version 1.0

import "../tasks/species_typing/task_cauris_cladetyper.wdl" as cauris_cladetyper
import "../tasks/gene_typing/task_snippy_variants.wdl" as snippy
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
    File read1
    File? read2
    File? ann_ref
    # Boolean paired_end = true
  }
  if (merlin_tag == "Candida auris") {
    call cauris_cladetyper.cauris_cladetyper as cladetyper {
      input: 
        assembly_fasta = assembly,
        samplename = samplename
        }
    call snippy.snippy_variants {
      input:
        reference = cladetyper.clade_spec_ref,
        read1 = read1,
        read2 = read2,
        query_gene = "FKS1",
        samplename = samplename
        }
    }
  output {
  # Candida Typing
  String? clade_type = cladetyper.gambit_cladetype
  String? cladetyper_analysis_date = cladetyper.date
  String? cladetyper_version = cladetyper.version
  String? cladetyper_docker_image = cladetyper.gambit_cladetyper_docker_image
  String? cladetype_annotated_ref = cladetyper.clade_spec_ref
  String? snippy_variants_version = snippy_variants.snippy_variants_version
  String? snippy_variants_query = snippy_variants.snippy_variants_query
  String? snippy_variants_hits = snippy_variants.snippy_variants_hits
  File? snippy_variants_gene_query_results = snippy_variants.snippy_variants_gene_query_results
  Array[File]? snippy_outputs = snippy_variants.snippy_outputs
  File? snippy_variants_results = snippy_variants.snippy_variants_results
  File? snippy_variants_bam = snippy_variants.snippy_variants_bam
  File? snippy_variants_bai = snippy_variants.snippy_variants_bai
  File? snippy_variants_summary = snippy_variants.snippy_variants_summary
 }
}