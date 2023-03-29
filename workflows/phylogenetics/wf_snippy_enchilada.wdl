version 1.0

import "../../workflows/standalone_modules/wf_snippy_variants.wdl" as snippy_variants_workflow
import "../../workflows/phylogenetics/wf_snippy_tree.wdl" as snippy_tree_workflow
import "../../tasks/task_versioning.wdl" as versioning

# input is arrays
# centroid takes in an array of assemblies
# reference seeker takes in one sample -- result from centroid
# get reference seeker genome
# scatter all input read arrays to snippy_variants
# gather outputs of snippy_variants as input to snippy_tree to run with reference seeker genome

workflow snippy_enchilada {
  input {
    Array[File] read1
    Array[File] read2
    Array[File] assembly_fasta
    Array[String] samplenames
    String tree_name
  }

  # call centroid {
  #   # add that here once merged
  # }
  # call reference_seeker {
  #   # add once merged
  # }
  # call ncbi_download {
  #   # add once merged 
  # }

  # see https://github.com/openwdl/wdl/issues/279 for syntax
  scatter (triplet in zip(zip(read1, read2), samplenames)) {
    call snippy_variants_workflow.snippy_variants_wf {
      input:
        read1 = triplet.left.left, # access the left-most object (read 1)
        read2 = triplet.left.right, # access the right-side object on the left (read 2)
        reference = ncbi_download.reference, 
        samplename = triplet.right # access the right-most object (samplename)
    }
  }
  call snippy_tree_workflow.snippy_tree_wf {
    input:
      tree_name = tree_name,
      snippy_variants_outdir_tarball = snippy_variants_wf.snippy_variants_outdir_tarball,
      samplenames = samplenames,
      reference = ncbi_download.reference
  }
  call versioning.version_capture {
    input:
  }
  output {
    # version capture
    String snippy_enchilada_version = version_capture.phb_version
    String snippy_enchilada_analysis_date = version_capture.date

    String snippy_enchilada_snippy_version = snippy_tree_wf.snippy_tree_snippy_version
    File snippy_enchilada_core_alignment = snippy_tree_wf.snippy_tree_core_alignment
    File snippy_enchilada_full_alignment = snippy_tree_wf.snippy_tree_full_alignment
    File snippy_enchilada_clean_full_alignment = snippy_tree_wf.snippy_tree_clean_full_alignment
    File snippy_enchilada_ref = snippy_tree_wf.snippy_tree_ref
    File snippy_enchilada_all_snps = snippy_tree_wf.snippy_tree_all_snps
    File snippy_enchilada_snps_summary = snippy_tree_wf.snippy_tree_snps_summary
    File snippy_enchilada_vcf = snippy_tree_wf.snippy_tree_vcf
    # iqtree outputs
    String snippy_enchilada_iqtree_version = snippy_tree_wf.snippy_tree_iqtree_version
    # snp_sites outputs
    String snp_sites_version = snippy_tree_wf.snp_sites_version
    # gubbins outputs
    String? snippy_enchilada_gubbins_version = snippy_tree_wf.snippy_tree_gubbins_version
    File? snippy_enchilada_gubbins_labelled_tree = snippy_tree_wf.snippy_tree_gubbins_labelled_tree
    File? snippy_enchilada_gubbins_polymorphic_fasta = snippy_tree_wf.snippy_tree_gubbins_polymorphic_fasta
    File? snippy_enchilada_gubbins_recombination_gff = snippy_tree_wf.snippy_tree_gubbins_recombination_gff
    File? snippy_enchilada_gubbins_branch_stats = snippy_tree_wf.snippy_tree_gubbins_branch_stats
    File? snippy_enchilada_gubbins_timetree = snippy_tree_wf.snippy_tree_gubbins_timetree
    File? snippy_enchilada_gubbins_timetree_stats = snippy_tree_wf.snippy_tree_gubbins_timetree_stats
    # snpdists outputs
    String snippy_enchilada_snpdists_version = snippy_tree_wf.snippy_tree_snpdists_version
    # reorder matrix outputs
    File snippy_enchilada_matrix = snippy_tree_wf.snippy_tree_matrix
    File snippy_enchilada_tree = snippy_tree_wf.snippy_tree_tree
    # data summary outputs
    File? snippy_enchilada_summarized_data = snippy_tree_wf.snippy_tree_summarized_data
  }

}