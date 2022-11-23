version 1.0

import "../tasks/gene_typing/task_snippy.wdl" as snippy
import "../tasks/task_versioning.wdl" as versioning

workflow snippy_variants_wf {
  meta {
    description: "Perform SNP analysis using snippy"
  }
  input {
    File reference
    File read1
    File? read2
    String samplename
  }
  call snippy.snippy_variants {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      reference = reference
  }
  call versioning.version_capture{
    input:
  }
  output {
    String snippy_variants_wf_version = version_capture.phbg_version
    String snippy_variants_version = snippy_variants.snippy_version
    String snippy_variant_query = snippy_variants.snippy_variant_query
    String snippy_variant_hits = snippy_variants.snippy_variant_hits
    File snippy_variant_gene_query_results = snippy_variants.snippy_variant_gene_query_results
    File snippy_variants_full_results = "~{samplename}/~{samplename}.csv"
  }
}
