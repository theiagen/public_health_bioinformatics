version 1.0


import "../../tasks/assembly/task_irma.wdl" as irma_task
import "../../tasks/gene_typing/drug_resistance/task_abricate.wdl" as abricate
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics
import "../../tasks/species_typing/orthomyxoviridae/task_genoflu.wdl" as genoflu_task
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../utilities/wf_influenza_antiviral_substitutions.wdl" as flu_antiviral
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults

workflow flu_track {
  meta {
    description: "This subworkflow contains all of the flu-specific modules to help organize the TheiaCoV workflows"
  }
  input {
    File read1
    File? read2
    String samplename

    String? flu_subtype

    String seq_method
    String standardized_organism

    # optional inputs to the tasks are available here due to Terra hiding them since they're within a subworkflow
    # IRMA inputs
    Boolean? irma_keep_ref_deletions
    # irma_min_consensus_support is defaulted to 30 for ILMN PE, 50 for ONT and both are set at workflow level WDLs and passed into flu_track subwf
    Int? irma_min_consensus_support
    Int? irma_min_read_length
    Int? irma_min_avg_consensus_allele_quality
    Float? irma_min_ambiguous_threshold
    String? irma_docker_image
    Int? irma_memory
    Int? irma_cpu
    Int? irma_disk_size

    # Assembly metrics inputs
    Int? assembly_metrics_memory
    Int? assembly_metrics_cpu
    Int? assembly_metrics_disk_size
    String? assembly_metrics_docker

    # GenoFLU inputs
    Float? genoflu_min_percent_identity
    File? genoflu_cross_reference
    Int? genoflu_cpu
    Int? genoflu_disk_size
    String? genoflu_docker
    Int? genoflu_memory

    # Abricate inputs
    Int? abricate_flu_min_percent_identity
    Int? abricate_flu_min_percent_coverage
    String? abricate_flu_docker
    Int? abricate_flu_memory
    Int? abricate_flu_cpu
    Int? abricate_flu_disk_size

    # flu antiviral substitutions subworkflow inputs
    File? flu_h1_ha_ref
    File? flu_h3_ha_ref
    File? flu_n1_na_ref
    File? flu_n2_na_ref
    File? flu_pa_ref
    File? flu_pb1_ref
    File? flu_pb2_ref
    File? flu_h1n1_m2_ref
    File? flu_h3n2_m2_ref
    String? antiviral_aa_subs

    # nextclade inputs
    String? nextclade_docker_image
    Int? nextclade_cpu
    Int? nextclade_memory
    Int? nextclade_disk_size
    File? nextclade_custom_input_dataset

    # nextclade output parser inputs
    String? nextclade_output_parser_docker
    Int? nextclade_output_parser_cpu
    Int? nextclade_output_parser_memory
    Int? nextclade_output_parser_disk_size
  }
  # IRMA will run if no assembly is provided (as in the case of TheiaCoV_FASTA)
  call irma_task.irma {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      seq_method = seq_method,
      minimum_consensus_support = irma_min_consensus_support,
      minimum_read_length = irma_min_read_length,
      minimum_average_consensus_allele_quality = irma_min_avg_consensus_allele_quality,
      minimum_ambiguous_threshold = irma_min_ambiguous_threshold,
      keep_ref_deletions = irma_keep_ref_deletions,
      docker = irma_docker_image,
      memory = irma_memory,
      cpu = irma_cpu,
      disk_size = irma_disk_size
  }
  # can be redone later to accomodate processing of HA and NA bams together in the task, perhaps with an organism flag
  if (defined(irma.seg_ha_bam)) {
    call assembly_metrics.stats_n_coverage as ha_assembly_coverage {
      input:
        bamfile = select_first([irma.seg_ha_bam]),
        samplename = samplename,
        memory = assembly_metrics_memory,
        cpu = assembly_metrics_cpu,
        disk_size = assembly_metrics_disk_size,
        docker = assembly_metrics_docker
    }
  }
  if (defined(irma.seg_na_bam)) {
    call assembly_metrics.stats_n_coverage as na_assembly_coverage {
      input:
        bamfile = select_first([irma.seg_na_bam]),
        samplename = samplename,
        memory = assembly_metrics_memory,
        cpu = assembly_metrics_cpu,
        disk_size = assembly_metrics_disk_size,
        docker = assembly_metrics_docker
    }
  }
  # combine HA & NA assembly coverages
  String ha_na_assembly_coverage_string = "HA: " + select_first([ha_assembly_coverage.depth, ""]) + ", NA: " + select_first([na_assembly_coverage.depth, ""])
  
  # combine HA & NA mapped reads percentages
  String ha_na_percentage_mapped_reads = "HA: " + select_first([ha_assembly_coverage.percentage_mapped_reads, ""]) + ", NA: " + select_first([na_assembly_coverage.percentage_mapped_reads, ""])

  # ABRICATE will run if assembly is provided, or was generated with IRMA
  if (defined(irma.irma_plurality_consensus_assemblies) && defined(irma.irma_assembly_fasta)){
    call abricate.abricate_flu {
      input:
        assembly = select_first([irma.irma_assembly_fasta]),
        samplename = samplename,
        min_percent_identity = abricate_flu_min_percent_identity,
        min_percent_coverage = abricate_flu_min_percent_coverage,
        cpu = abricate_flu_cpu,
        memory = abricate_flu_memory,
        docker = abricate_flu_docker,
        disk_size = abricate_flu_disk_size
    }
    # check IRMA subtype content if IRMA was run
    if (defined(irma.irma_subtype)) {
      # if IRMA cannot predict a subtype (like with Flu B samples), then set the flu_subtype to the abricate_flu_subtype String output (e.g. "Victoria" for Flu B)
      String algorithmic_flu_subtype = if irma.irma_subtype == "No subtype predicted by IRMA" then abricate_flu.abricate_flu_subtype else irma.irma_subtype
    }
    call set_organism_defaults.organism_parameters as set_flu_na_nextclade_values {
      input:
        organism = standardized_organism,
        flu_segment = "NA",
        flu_subtype = select_first([flu_subtype, algorithmic_flu_subtype, abricate_flu.abricate_flu_subtype, "N/A"]),
    }
    call set_organism_defaults.organism_parameters as set_flu_ha_nextclade_values {
      input:
        organism = standardized_organism,
        flu_segment = "HA",
        flu_subtype = select_first([flu_subtype, algorithmic_flu_subtype, abricate_flu.abricate_flu_subtype, "N/A"])
    }
    # these are necessary because these are optional values and cannot be directly compared in before the nextclade task. 
    # checking for variable definition can be done though, which is why we create variables here
    if (set_flu_na_nextclade_values.nextclade_dataset_tag == "NA") { # this "NA" is Not Applicable, not the NA segment
      Boolean do_not_run_flu_na_nextclade = true
    }
    if (set_flu_ha_nextclade_values.nextclade_dataset_tag == "NA") {
      Boolean do_not_run_flu_ha_nextclade = true
    }
  }       
  # if IRMA was run successfully, run the flu_antiviral substitutions task 
  # this block must be placed beneath the previous block because it is used in this subworkflow
  if (defined(irma.irma_plurality_consensus_assemblies)) {
    call flu_antiviral.flu_antiviral_substitutions {
      input:
        na_segment_assembly = irma.seg_na_assembly_padded,
        ha_segment_assembly = irma.seg_ha_assembly_padded,
        pa_segment_assembly = irma.seg_pa_assembly_padded,
        pb1_segment_assembly = irma.seg_pb1_assembly_padded,
        pb2_segment_assembly = irma.seg_pb2_assembly_padded,
        mp_segment_assembly = irma.seg_mp_assembly_padded,
        abricate_flu_subtype = select_first([abricate_flu.abricate_flu_subtype, ""]),
        irma_flu_subtype = irma.irma_subtype,
        antiviral_aa_subs = antiviral_aa_subs,
        flu_h1_ha_ref = flu_h1_ha_ref,
        flu_h3_ha_ref = flu_h3_ha_ref,
        flu_n1_na_ref = flu_n1_na_ref,
        flu_n2_na_ref = flu_n2_na_ref,
        flu_pa_ref = flu_pa_ref,
        flu_pb1_ref = flu_pb1_ref,
        flu_pb2_ref = flu_pb2_ref,
        flu_h1n1_m2_ref = flu_h1n1_m2_ref,
        flu_h3n2_m2_ref = flu_h3n2_m2_ref
    }
  }
  if (defined(irma.seg_ha_assembly) && ! defined(do_not_run_flu_ha_nextclade)) {
    call nextclade_task.nextclade_v3 as nextclade_flu_ha {
      input:
        genome_fasta = select_first([irma.seg_ha_assembly]),
        dataset_name = select_first([set_flu_ha_nextclade_values.nextclade_dataset_name]),
        dataset_tag = select_first([set_flu_ha_nextclade_values.nextclade_dataset_tag]),
        docker = nextclade_docker_image,
        cpu = nextclade_cpu,
        memory = nextclade_memory,
        disk_size = nextclade_disk_size
    }
    call nextclade_task.nextclade_output_parser as nextclade_output_parser_flu_ha {
      input:
        nextclade_tsv = nextclade_flu_ha.nextclade_tsv,
        organism = standardized_organism,
        docker = nextclade_output_parser_docker,
        cpu = nextclade_output_parser_cpu,
        memory = nextclade_output_parser_memory,
        disk_size = nextclade_output_parser_disk_size
    }
  }
  if (defined(irma.seg_na_assembly) && ! defined(do_not_run_flu_na_nextclade)) {
    call nextclade_task.nextclade_v3 as nextclade_flu_na {
      input:
        genome_fasta = select_first([irma.seg_na_assembly]),
        dataset_name = select_first([set_flu_na_nextclade_values.nextclade_dataset_name]),
        dataset_tag = select_first([set_flu_na_nextclade_values.nextclade_dataset_tag]),
        docker = nextclade_docker_image,
        cpu = nextclade_cpu,
        memory = nextclade_memory,
        disk_size = nextclade_disk_size
    }
    call nextclade_task.nextclade_output_parser as nextclade_output_parser_flu_na {
      input:
        nextclade_tsv = nextclade_flu_na.nextclade_tsv,
        organism = standardized_organism,
        docker = nextclade_output_parser_docker,
        cpu = nextclade_output_parser_cpu,
        memory = nextclade_output_parser_memory,
        disk_size = nextclade_output_parser_disk_size
    }
  }
  # only run GenoFLU and custom nextclade dataset if the subtype is H5N1 and the clade is 2.3.4.4b as they are specific to this subtype and clade.
  if (select_first([flu_subtype, algorithmic_flu_subtype, abricate_flu.abricate_flu_subtype, "N/A"]) == "H5N1" && select_first([nextclade_output_parser_flu_ha.nextclade_clade, ""]) == "2.3.4.4b") {
    call genoflu_task.genoflu {
      input:
        assembly_fasta = select_first([irma.irma_assembly_fasta]),
        samplename = samplename,
        min_percent_identity = genoflu_min_percent_identity,
        cross_reference = genoflu_cross_reference,
        cpu = genoflu_cpu,
        disk_size = genoflu_disk_size,
        docker = genoflu_docker,
        memory = genoflu_memory
    }
    if (genoflu.genoflu_genotype == "B3.13" && defined(nextclade_custom_input_dataset)) {
      call nextclade_task.nextclade_v3 as nextclade_flu_h5n1 {
        input:
          genome_fasta = select_first([irma.irma_assembly_fasta_concatenated]),
          custom_input_dataset = nextclade_custom_input_dataset,
          docker = nextclade_docker_image,
          cpu = nextclade_cpu,
          memory = nextclade_memory,
          disk_size = nextclade_disk_size
      }
      call nextclade_task.nextclade_output_parser as nextclade_output_parser_flu_h5n1 {
        input:
          nextclade_tsv = nextclade_flu_h5n1.nextclade_tsv,
          organism = standardized_organism,
          docker = nextclade_output_parser_docker,
          cpu = nextclade_output_parser_cpu,
          memory = nextclade_output_parser_memory,
          disk_size = nextclade_output_parser_disk_size
      }
    }
  }
  output {
    # IRMA outputs 
    String irma_version = irma.irma_version
    String irma_docker = irma.irma_docker
    Int irma_minimum_consensus_support = irma.irma_minimum_consensus_support
    String irma_type = irma.irma_type
    String irma_subtype = irma.irma_subtype
    String irma_subtype_notes = irma.irma_subtype_notes
    File? irma_assembly_fasta = irma.irma_assembly_fasta
    File? irma_assembly_fasta_concatenated = irma.irma_assembly_fasta_concatenated
    File? irma_assembly_fasta_padded = irma.irma_assembly_fasta_padded
    File? irma_assembly_fasta_concatenated_padded = irma.irma_assembly_fasta_concatenated_padded
    File? irma_ha_segment_fasta = irma.seg_ha_assembly
    File? irma_na_segment_fasta = irma.seg_na_assembly
    File? irma_pa_segment_fasta = irma.seg_pa_assembly
    File? irma_pb1_segment_fasta = irma.seg_pb1_assembly
    File? irma_pb2_segment_fasta = irma.seg_pb2_assembly
    File? irma_mp_segment_fasta = irma.seg_mp_assembly
    File? irma_np_segment_fasta = irma.seg_np_assembly
    File? irma_ns_segment_fasta = irma.seg_ns_assembly
    Array[File] irma_assemblies = irma.irma_plurality_consensus_assemblies
    Array[File] irma_vcfs = irma.irma_vcfs
    Array[File] irma_bams = irma.irma_bams
    File? irma_ha_bam = irma.seg_ha_bam
    File? irma_na_bam = irma.seg_na_bam
    String ha_na_assembly_coverage = ha_na_assembly_coverage_string
     # calulate mapped reads percentage for flu samples
    String percentage_mapped_reads = ha_na_percentage_mapped_reads
    # GenoFLU outputs
    String? genoflu_version = genoflu.genoflu_version
    String? genoflu_genotype = genoflu.genoflu_genotype
    String? genoflu_all_segments = genoflu.genoflu_all_segments
    File? genoflu_output_tsv = genoflu.genoflu_output_tsv
    # Abricate outputs
    String? abricate_flu_type = abricate_flu.abricate_flu_type
    String? abricate_flu_subtype =  abricate_flu.abricate_flu_subtype
    File? abricate_flu_results = abricate_flu.abricate_flu_results
    String? abricate_flu_database =  abricate_flu.abricate_flu_database
    String? abricate_flu_version = abricate_flu.abricate_flu_version
    # Nextclade outputs
    String? nextclade_version = nextclade_flu_ha.nextclade_version
    String? nextclade_docker = nextclade_flu_ha.nextclade_docker
    # Nextclade H5N1 outputs
    File? nextclade_json_flu_h5n1 = nextclade_flu_h5n1.nextclade_json
    File? auspice_json_flu_h5n1 = nextclade_flu_h5n1.auspice_json
    File? nextclade_tsv_flu_h5n1 = nextclade_flu_h5n1.nextclade_tsv
    String? nextclade_aa_subs_flu_h5n1 = nextclade_output_parser_flu_h5n1.nextclade_aa_subs
    String? nextclade_aa_dels_flu_h5n1 = nextclade_output_parser_flu_h5n1.nextclade_aa_dels
    String? nextclade_clade_flu_h5n1 = nextclade_output_parser_flu_h5n1.nextclade_clade
    String? nextclade_qc_flu_h5n1 = nextclade_output_parser_flu_h5n1.nextclade_qc
    # Nextclade HA outputs
    File? nextclade_json_flu_ha = nextclade_flu_ha.nextclade_json
    File? auspice_json_flu_ha =  nextclade_flu_ha.auspice_json
    File? nextclade_tsv_flu_ha = nextclade_flu_ha.nextclade_tsv
    String? nextclade_ds_tag_flu_ha = set_flu_ha_nextclade_values.nextclade_dataset_tag
    String? nextclade_aa_subs_flu_ha = nextclade_output_parser_flu_ha.nextclade_aa_subs
    String? nextclade_aa_dels_flu_ha = nextclade_output_parser_flu_ha.nextclade_aa_dels
    String? nextclade_clade_flu_ha = nextclade_output_parser_flu_ha.nextclade_clade
    String? nextclade_qc_flu_ha = nextclade_output_parser_flu_ha.nextclade_qc
    # Nextclade NA outputs
    File? nextclade_json_flu_na = nextclade_flu_na.nextclade_json
    File? auspice_json_flu_na = nextclade_flu_na.auspice_json
    File? nextclade_tsv_flu_na = nextclade_flu_na.nextclade_tsv
    String? nextclade_ds_tag_flu_na = set_flu_na_nextclade_values.nextclade_dataset_tag
    String? nextclade_aa_subs_flu_na = nextclade_output_parser_flu_na.nextclade_aa_subs
    String? nextclade_aa_dels_flu_na = nextclade_output_parser_flu_na.nextclade_aa_dels
    String? nextclade_clade_flu_na = nextclade_output_parser_flu_na.nextclade_clade
    String? nextclade_qc_flu_na = nextclade_output_parser_flu_na.nextclade_qc
    # Flu Antiviral Substitution Outputs
    String? flu_A_315675_resistance = flu_antiviral_substitutions.flu_A_315675_resistance
    String? flu_amantadine_resistance = flu_antiviral_substitutions.flu_amantadine_resistance
    String? flu_compound_367_resistance = flu_antiviral_substitutions.flu_compound_367_resistance
    String? flu_favipiravir_resistance = flu_antiviral_substitutions.flu_favipiravir_resistance
    String? flu_fludase_resistance = flu_antiviral_substitutions.flu_fludase_resistance
    String? flu_L_742_001_resistance = flu_antiviral_substitutions.flu_L_742_001_resistance
    String? flu_laninamivir_resistance = flu_antiviral_substitutions.flu_laninamivir_resistance
    String? flu_peramivir_resistance = flu_antiviral_substitutions.flu_peramivir_resistance
    String? flu_pimodivir_resistance = flu_antiviral_substitutions.flu_pimodivir_resistance
    String? flu_rimantadine_resistance = flu_antiviral_substitutions.flu_rimantadine_resistance
    String? flu_oseltamivir_resistance = flu_antiviral_substitutions.flu_oseltamivir_resistance
    String? flu_xofluza_resistance = flu_antiviral_substitutions.flu_xofluza_resistance
    String? flu_zanamivir_resistance = flu_antiviral_substitutions.flu_zanamivir_resistance
  }
}
