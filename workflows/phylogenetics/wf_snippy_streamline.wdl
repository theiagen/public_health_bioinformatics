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
    ### version capture ###
    String snippy_streamline_version = version_capture.phb_version
    String snippy_streamline_analysis_date = version_capture.date

    ### centroid outputs ###
    String? ss_centroid_genome_samplename = centroid.centroid_genome_samplename
    File? snippy_streamline_centroid_genome_fasta = centroid.centroid_genome_fasta
    File? snippy_streamline_centroid_mash_tsv = centroid.centroid_mash_tsv
    
    ### referenceseeker outputs ###
    String? snippy_streamline_referenceseeker_top_hit_ncbi_accession = referenceseeker.referenceseeker_top_hit_ncbi_accession
    String? snippy_streamline_referenceseeker_version = referenceseeker.referenceseeker_version
    File? snippy_streamline_referenceseeker_tsv = referenceseeker.referenceseeker_tsv
    String? snippy_streamline_referenceseeker_docker = referenceseeker.referenceseeker_docker
    String? snippy_streamline_referenceseeker_database = referenceseeker.referenceseeker_database

    ### ncbi datasets outputs ###
    # output reference FASTA file now consolidated to File snippy_ref
    File? snippy_streamline_ref_metadata_json = ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_data_report_json
    String? snippy_streamline_ncbi_datasets_version = ncbi_datasets_download_genome_accession.ncbi_datasets_version
    String? snippy_streamline_ncbi_datasets_docker = ncbi_datasets_download_genome_accession.ncbi_datasets_docker

    ### snippy_variants wf outputs ###
    Array[File] snippy_variants_outdir_tarball = snippy_variants_wf.snippy_variants_outdir_tarball
    Array[String] snippy_variants_snippy_version = snippy_variants_wf.snippy_version
    Array[String] snippy_variants_snippy_docker = snippy_variants_wf.snippy_docker

    ### snippy_tree wf outputs ###
    String snippy_tree_snippy_version = snippy_tree_wf.snippy_tree_snippy_version
    String snippy_tree_snippy_docker = snippy_tree_wf.snippy_tree_snippy_docker
    String snippy_iqtree2_docker = snippy_tree_wf.snippy_iqtree2_docker
    String snippy_iqtree2_version = snippy_tree_wf.snippy_iqtree2_version
    String snippy_iqtree2_model_used = snippy_tree_wf.snippy_iqtree2_model_used
    File snippy_final_alignment = snippy_tree_wf.snippy_final_alignment
    File snippy_final_tree = snippy_tree_wf.snippy_final_tree
    File snippy_ref = select_first([snippy_tree_wf.snippy_ref, ncbi_datasets_download_genome_accession.ncbi_datasets_assembly_fasta])
    File snippy_msa_snps_summary = snippy_tree_wf.snippy_msa_snps_summary
    String? snippy_snp_sites_version = snippy_tree_wf.snippy_snp_sites_version
    String? snippy_snp_sites_docker = snippy_tree_wf.snippy_snp_sites_docker
    String? snippy_gubbins_version = snippy_tree_wf.snippy_gubbins_version
    String? snippy_gubbins_docker = snippy_tree_wf.snippy_gubbins_docker
    File? snippy_gubbins_recombination_gff = snippy_tree_wf.snippy_gubbins_recombination_gff
    File? snippy_gubbins_branch_stats = snippy_tree_wf.snippy_gubbins_branch_stats
    String snippy_snp_dists_version = snippy_tree_wf.snippy_snp_dists_version
    String snippy_snp_dists_docker = snippy_tree_wf.snippy_snp_dists_docker
    File snippy_snp_matrix = snippy_tree_wf.snippy_snp_matrix
    File? snippy_summarized_data = snippy_tree_wf.snippy_summarized_data
  }
}