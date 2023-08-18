version 1.0

import "../../tasks/gene_typing/task_snippy_gene_query.wdl" as snippy_gene_query_task
import "../../tasks/gene_typing/task_snippy_variants.wdl" as snippy
import "../../tasks/task_versioning.wdl" as versioning

workflow snippy_variants_wf {
  meta {
    description: "Perform SNP analysis using snippy"
  }
  input {
    File reference_genome_file
    File read1
    File? read2
    String samplename
    String? query_gene
    # optional inputs are exposed here so that they can be set by user in higher-level workflows. example: snippy_streamline
    # a bit of duplicate code, but done for a reason!
    String? docker
    Int? cpus
    Int? memory
    Int? map_qual
    Int? base_quality
    Int? min_coverage
    Float? min_frac
    Int? min_quality
    Int? maxsoft
  }
  call snippy.snippy_variants {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      reference_genome_file = reference_genome_file,
      docker = docker,
      cpus = cpus,
      memory = memory,
      map_qual = map_qual,
      base_quality = base_quality,
      min_coverage = min_coverage,
      min_frac = min_frac,
      min_quality = min_quality,
      maxsoft = maxsoft
  }
  if ("~{query_gene}" != "") {
    call snippy_gene_query_task.snippy_gene_query {
      input:
        samplename = samplename,
        snippy_variants_results = snippy_variants.snippy_variants_results,
        reference = reference_genome_file,
        query_gene = query_gene
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    String snippy_variants_wf_version = version_capture.phb_version
    String snippy_variants_version = snippy_variants.snippy_variants_version
    String snippy_variants_docker = snippy_variants.snippy_variants_docker
    File snippy_variants_results = snippy_variants.snippy_variants_results
    File snippy_variants_bam = snippy_variants.snippy_variants_bam
    File snippy_variants_bai = snippy_variants.snippy_variants_bai
    File snippy_variants_summary = snippy_variants.snippy_variants_summary
    File snippy_variants_outdir_tarball = snippy_variants.snippy_variants_outdir_tarball
    # snippy gene query outputs
    String? snippy_variants_query = snippy_gene_query.snippy_variants_query
    String? snippy_variants_query_check = snippy_gene_query.snippy_variants_query_check
    String? snippy_variants_hits = snippy_gene_query.snippy_variants_hits
    File? snippy_variants_gene_query_results = snippy_gene_query.snippy_variants_gene_query_results
  }
}
