version 1.0

import "../../workflows/standalone_modules/wf_snippy_variants.wdl" as snippy_variants_workflow
import "../../workflows/phylogenetics/wf_snippy_tree.wdl" as snippy_tree_workflow
import "../../tasks/phylogenetic_inference/task_centroid.wdl" as centroid_task
import "../../tasks/phylogenetic_inference/task_referenceseeker.wdl" as referenceseeker_task
import "../../tasks/utilities/task_ncbi_datasets.wdl" as ncbi_datasets_task
import "../../tasks/task_versioning.wdl" as versioning

# input is arrays
# centroid takes in an array of assemblies
# reference seeker takes in one sample -- result from centroid
# get reference seeker genome
# scatter all input read arrays to snippy_variants
# gather outputs of snippy_variants as input to snippy_tree to run with reference seeker genome

workflow snippy_streamline {
  input {
    Array[File] read1
    Array[File] read2
    Array[File] assembly_fasta
    Array[String] samplenames
    String tree_name
    # this input file can be a FASTA or GBK
    File? reference_genome_file
  }
  # if user does not provide reference genome fasta, determine one for the user by running, centroid, referenceseeker and ncbi datasets to acquire one
  if(! defined(reference_genome_file)){
  call centroid_task.centroid {
    input:
      assembly_fasta = assembly_fasta
  }
  call referenceseeker_task.referenceseeker {
    input:
      assembly_fasta = centroid.centroid_genome_fasta,
      samplename = centroid.centroid_genome_samplename
  }
  call ncbi_datasets_task.ncbi_datasets_download_genome_accession {
    input:
      ncbi_accession = referenceseeker.referenceseeker_top_hit_ncbi_accession
  }
  }
  # see https://github.com/openwdl/wdl/issues/279 for syntax explanation
  # see also https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#arraypairxy-ziparrayx-arrayy for zip explanation
  scatter (triplet in zip(zip(read1, read2), samplenames)) {
    call snippy_variants_workflow.snippy_variants_wf {
      input:
        read1 = triplet.left.left, # access the left-most object (read 1)
        read2 = triplet.left.right, # access the right-side object on the left (read 2)
        reference_genome_file = select_first([reference_genome_file,ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_fasta]),
        samplename = triplet.right # access the right-most object (samplename)
    }
  }
  call snippy_tree_workflow.snippy_tree_wf {
    input:
      tree_name = tree_name,
      snippy_variants_outdir_tarball = snippy_variants_wf.snippy_variants_outdir_tarball,
      samplenames = samplenames,
      reference_genome_file = select_first([reference_genome_file,ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_fasta])
  }
  call versioning.version_capture {
    input:
  }
  output {
    # version capture
    String snippy_streamline_version = version_capture.phb_version
    String snippy_streamline_analysis_date = version_capture.date
    # centroid outputs
    String? snippy_streamline_centroid_genome_samplename = centroid.centroid_genome_samplename
    File? snippy_streamline_centroid_genome_fasta = centroid.centroid_genome_fasta
    File? snippy_streamline_centroid_mash_tsv = centroid.centroid_mash_tsv
    # referenceseeker outputs
    String? snippy_streamline_referenceseeker_top_hit_ncbi_accession = referenceseeker.referenceseeker_top_hit_ncbi_accession
    String? snippy_streamline_referenceseeker_version = referenceseeker.referenceseeker_version
    File? snippy_streamline_referenceseeker_tsv = referenceseeker.referenceseeker_tsv
    String? snippy_streamline_referenceseeker_docker = referenceseeker.referenceseeker_docker
    String? snippy_streamline_referenceseeker_database = referenceseeker.referenceseeker_database
    # ncbi datasets outputs
    File? snippy_streamline_ncbi_datasets_assembly_data_report_json = ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_data_report_json
    File? snippy_streamline_ncbi_datasets_assembly_fasta = ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_fasta 
    String? snippy_streamline_ncbi_datasets_version = ncbi_datasets_download_genome_accession.ncbi_datasets_version
    String? snippy_streamline_ncbi_datasets_docker = ncbi_datasets_download_genome_accession.ncbi_datasets_docker
    # snippy_variants output
    Array[File] snippy_streamline_snippy_variants_outdir_tarball = snippy_variants_wf.snippy_variants_outdir_tarball
    # snippy_tree version
    String snippy_streamline_snippy_version = snippy_tree_wf.snippy_tree_snippy_version
    # snippy_tree alignments
    File snippy_streamline_core_alignment = snippy_tree_wf.snippy_tree_core_alignment
    File snippy_streamline_full_alignment = snippy_tree_wf.snippy_tree_full_alignment
    File snippy_streamline_clean_full_alignment = snippy_tree_wf.snippy_tree_clean_full_alignment
    # snippy_tree reference file
    File snippy_streamline_ref = snippy_tree_wf.snippy_tree_ref
    # snippy_tree variant outputs
    File snippy_streamline_all_snps = snippy_tree_wf.snippy_tree_all_snps
    File snippy_streamline_snps_summary = snippy_tree_wf.snippy_tree_snps_summary
    File snippy_streamline_vcf = snippy_tree_wf.snippy_tree_vcf
    # iqtree outputs
    String snippy_streamline_iqtree_version = snippy_tree_wf.snippy_tree_iqtree_version
    # snp_sites outputs
    String snp_sites_version = snippy_tree_wf.snp_sites_version
    # gubbins outputs
    String? snippy_streamline_gubbins_version = snippy_tree_wf.snippy_tree_gubbins_version
    File? snippy_streamline_gubbins_labelled_tree = snippy_tree_wf.snippy_tree_gubbins_labelled_tree
    File? snippy_streamline_gubbins_polymorphic_fasta = snippy_tree_wf.snippy_tree_gubbins_polymorphic_fasta
    File? snippy_streamline_gubbins_recombination_gff = snippy_tree_wf.snippy_tree_gubbins_recombination_gff
    File? snippy_streamline_gubbins_branch_stats = snippy_tree_wf.snippy_tree_gubbins_branch_stats
    File? snippy_streamline_gubbins_timetree = snippy_tree_wf.snippy_tree_gubbins_timetree
    File? snippy_streamline_gubbins_timetree_stats = snippy_tree_wf.snippy_tree_gubbins_timetree_stats
    # snpdists outputs
    String snippy_streamline_snpdists_version = snippy_tree_wf.snippy_tree_snpdists_version
    # reorder matrix outputs
    File snippy_streamline_matrix = snippy_tree_wf.snippy_tree_matrix
    File snippy_streamline_tree = snippy_tree_wf.snippy_tree_tree
    # data summary outputs
    File? snippy_streamline_summarized_data = snippy_tree_wf.snippy_tree_summarized_data
  }

}