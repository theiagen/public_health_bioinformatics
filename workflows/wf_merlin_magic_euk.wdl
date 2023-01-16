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
    # Boolean paired_end = true
  }
  if (merlin_tag == "Candida auris") {
    call cauris_cladetyper.cauris_cladetyper as cladetyper {
      input: 
        assembly_fasta = assembly,
        samplename = samplename
        }
    call snippy.snippy_variants as snippy_cauris {
      input:
        reference = cladetyper.clade_spec_ref,
        read1 = read1,
        read2 = read2,
        query_gene = "FKS1,ERG11,FUR1",
        samplename = samplename
        }
    }
  if (merlin_tag == "Candida albicans") {
    call snippy.snippy_variants as snippy_calbicans {
      input:
        reference = "gs://theiagen-public-files/terra/theiaeuk_files/Candida_albicans_GCF_000182965.3_ASM18296v3_genomic.gbff",
        read1 = read1,
        read2 = read2,
        query_gene = "ERG11,FKS1,FUR1,RTA2",
        samplename = samplename
        }
    }
  if (merlin_tag == "Aspergillus fumigatus") {
    call snippy.snippy_variants as snippy_afumigatus {
      input:
        reference = "gs://theiagen-public-files/terra/theiaeuk_files/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.gbff",
        read1 = read1,
        read2 = read2,
        query_gene = "CYP51a,HAPE,COX10",
        samplename = samplename
        }
    }
  if (merlin_tag == "Cryptococcus neoformans") {
    call snippy.snippy_variants as snippy_crypto {
      input:
        reference = "gs://theiagen-public-files/terra/theiaeuk_files/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.gbff",
        read1 = read1,
        read2 = read2,
        query_gene = "ERG11",
        samplename = samplename
    }
  }
  output {
  # Typing
  String? clade_type = cladetyper.gambit_cladetype
  String? cladetyper_analysis_date = cladetyper.date
  String? cladetyper_version = cladetyper.version
  String? cladetyper_docker_image = cladetyper.gambit_cladetyper_docker_image
  String? cladetype_annotated_ref = cladetyper.clade_spec_ref
  String? snippy_variants_version = select_first([snippy_cauris.snippy_variants_version, snippy_calbicans.snippy_variants_version, snippy_afumigatus.snippy_variants_version, snippy_crypto.snippy_variants_version])
  String? snippy_variants_query = select_first([snippy_cauris.snippy_variants_query, snippy_calbicans.snippy_variants_query, snippy_afumigatus.snippy_variants_query, snippy_crypto.snippy_variants_query])
  String? snippy_variants_hits = select_first([snippy_cauris.snippy_variants_hits, snippy_calbicans.snippy_variants_hits, snippy_afumigatus.snippy_variants_hits, snippy_crypto.snippy_variants_hits])
  File? snippy_variants_gene_query_results = select_first([snippy_cauris.snippy_variants_gene_query_results, snippy_calbicans.snippy_variants_gene_query_results, snippy_afumigatus.snippy_variants_gene_query_results, snippy_crypto.snippy_variants_gene_query_results])
  Array[File]? snippy_outputs = select_first([snippy_cauris.snippy_outputs, snippy_calbicans.snippy_outputs, snippy_afumigatus.snippy_outputs, snippy_crypto.snippy_outputs])
  File? snippy_variants_results = select_first([snippy_cauris.snippy_variants_results, snippy_calbicans.snippy_variants_results, snippy_afumigatus.snippy_variants_results, snippy_crypto.snippy_variants_results])
  File? snippy_variants_bam = select_first([snippy_cauris.snippy_variants_bam, snippy_calbicans.snippy_variants_bam, snippy_afumigatus.snippy_variants_bam, snippy_crypto.snippy_variants_bam])
  File? snippy_variants_bai = select_first([snippy_cauris.snippy_variants_bai, snippy_calbicans.snippy_variants_bai, snippy_afumigatus.snippy_variants_bai, snippy_crypto.snippy_variants_bai])
  File? snippy_variants_summary = select_first([snippy_cauris.snippy_variants_summary, snippy_calbicans.snippy_variants_summary, snippy_afumigatus.snippy_variants_summary, snippy_crypto.snippy_variants_summary])
 }
}