version 1.0

# theiaprok
import "../../tasks/species_typing/task_serotypefinder.wdl" as serotypefinder_task
import "../../tasks/species_typing/task_ectyper.wdl" as ectyper_task
import "../../tasks/species_typing/task_shigatyper.wdl" as shigatyper_task
import "../../tasks/species_typing/task_shigeifinder.wdl" as shigeifinder_task
import "../../tasks/species_typing/task_sonneityping.wdl" as sonneityping_task
import "../../tasks/species_typing/task_lissero.wdl" as lissero_task
import "../../tasks/species_typing/task_sistr.wdl" as sistr_task
import "../../tasks/species_typing/task_seqsero2.wdl" as seqsero2_task
import "../../tasks/species_typing/task_kleborate.wdl" as kleborate_task
import "../../tasks/species_typing/task_tbprofiler.wdl" as tbprofiler_task
import "../../tasks/species_typing/task_legsta.wdl" as legsta_task
import "../../tasks/species_typing/task_genotyphi.wdl" as genotyphi
import "../../tasks/species_typing/task_kaptive.wdl" as kaptive_task
import "../../tasks/species_typing/task_ngmaster.wdl" as ngmaster_task
import "../../tasks/species_typing/task_meningotype.wdl" as meningotype_task
import "../../tasks/species_typing/task_spatyper.wdl" as spatyper_task
import "../../tasks/species_typing/task_staphopiasccmec.wdl" as staphopia_sccmec_task
import "../../tasks/species_typing/task_agrvate.wdl" as agrvate_task
import "../../tasks/species_typing/task_seroba.wdl" as seroba
import "../../tasks/species_typing/task_pbptyper.wdl" as pbptyper
import "../../tasks/species_typing/task_poppunk_streppneumo.wdl" as poppunk_spneumo
import "../../tasks/species_typing/task_pasty.wdl" as pasty_task
import "../../tasks/gene_typing/task_abricate.wdl" as abricate_task
import "../../tasks/species_typing/task_emmtypingtool.wdl" as emmtypingtool_task
import "../../tasks/species_typing/task_hicap.wdl" as hicap_task
import "../../tasks/species_typing/task_srst2_vibrio.wdl" as srst2_vibrio_task
import "../../tasks/species_typing/task_virulencefinder.wdl" as virulencefinder_task

# theiaeuk
import "../../tasks/species_typing/task_cauris_cladetyper.wdl" as cauris_cladetyper
import "../../tasks/gene_typing/task_snippy_variants.wdl" as snippy
import "../../tasks/gene_typing/task_snippy_gene_query.wdl" as snippy_gene_query

workflow merlin_magic {
  meta {
    description: "Workflow for bacterial and fungal species typing; based on the Bactopia subworkflow Merlin (https://bactopia.github.io/bactopia-tools/merlin/)"
  }
  input {
    String samplename
    String merlin_tag
    File assembly
    File? read1
    File? read2
    Int? pasty_min_pident
    Int? pasty_min_coverage
    String? hicap_docker_image
    String? pasty_docker_image
    String? emmtypingtool_docker_image
    String? shigeifinder_docker_image
    String? staphopia_sccmec_docker_image
    String? agrvate_docker_image
    String? virulencefinder_docker_image
    Boolean paired_end = true
    Boolean call_poppunk = true
    Boolean ont_data = false
    Boolean call_shigeifinder_reads_input = false
    Boolean assembly_only = false
    Boolean theiaeuk = false
    Boolean tbprofiler_additional_outputs = false
    String output_seq_method_type = "WGS"
    String? snippy_query_gene
    Int srst2_min_cov = 80
    Int srst2_max_divergence = 20
    Int srst2_min_depth = 5
    Int srst2_min_edge_depth = 2
    Int srst2_gene_max_mismatch = 2000
    Float? virulencefinder_coverage_threshold
    Float? virulencefinder_identity_threshold
    String? virulencefinder_database
  }
  # theiaprok
  if (merlin_tag == "Acinetobacter baumannii") {
    call kaptive_task.kaptive {
      input:
        assembly = assembly,
        samplename = samplename
    }
    call abricate_task.abricate {
      input:
        assembly = assembly,
        samplename = samplename,
        database = "AcinetobacterPlasmidTyping",
        minid = 95 # strict threshold of 95% identity for typing purposes
    }
  }
  if (merlin_tag == "Escherichia" || merlin_tag == "Shigella_sonnei" ) {
    # tools specific to all Escherichia and Shigella species
    call serotypefinder_task.serotypefinder {
      input:
        assembly = assembly,
        samplename = samplename
    }
    call ectyper_task.ectyper {
      input:
        assembly = assembly,
        samplename = samplename
    }
    if (!assembly_only) {
      call shigatyper_task.shigatyper { # test ONT compatibility
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          read1_is_ont = ont_data
      }
    }
    call shigeifinder_task.shigeifinder {
      input:
        assembly = assembly,
        samplename = samplename,
        docker = shigeifinder_docker_image
    }
    if (call_shigeifinder_reads_input && !assembly_only && !ont_data) { # illumina only
      call shigeifinder_task.shigeifinder_reads as shigeifinder_reads {
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          docker = shigeifinder_docker_image,
          paired_end = paired_end
      }
    }
  }
  if (merlin_tag == "Escherichia" ) {
    # E coli specific tasks
    call virulencefinder_task.virulencefinder {
      input:
      #  read1 = read1,
      #  read2 = read2,
        assembly = assembly,
        samplename = samplename,
      #  paired_end = paired_end,
      #  assembly_only = assembly_only,
      #  ont_data = ont_data,
        docker = virulencefinder_docker_image,
        coverage_threshold = virulencefinder_coverage_threshold,
        identity_threshold = virulencefinder_identity_threshold,
        database = virulencefinder_database
    }
  }
  if (merlin_tag == "Shigella_sonnei") {
    # Shigella sonnei specific tasks
    if (!assembly_only) {
      call sonneityping_task.sonneityping { # test ONT compatibility
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          ont_data = ont_data
      }
    }
  }
  if (merlin_tag == "Listeria") {
    call lissero_task.lissero {
      input:
        assembly = assembly,
        samplename = samplename
    }
  }
  if (merlin_tag == "Salmonella") {
    call sistr_task.sistr {
      input: 
        assembly = assembly,
        samplename = samplename
    }
    if (!ont_data && !assembly_only) {
      call seqsero2_task.seqsero2 { 
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          paired_end = paired_end
      }
    }
    if (ont_data || assembly_only) {
      call seqsero2_task.seqsero2_assembly {
        input:
          assembly_fasta = assembly,
          samplename = samplename
      }
    }
    if ((select_first([seqsero2.seqsero2_predicted_serotype, seqsero2_assembly.seqsero2_predicted_serotype]) == "Typhi" || sistr.sistr_predicted_serotype == "Typhi") && !assembly_only) {
      call genotyphi.genotyphi as genotyphi_task { # needs testing
        input: 
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          ont_data = ont_data
      }
    }
  }
  if (merlin_tag == "Klebsiella") {
    call kleborate_task.kleborate {
      input:
        assembly = assembly,
        samplename = samplename
    }
  }
  if (merlin_tag == "Neisseria gonorrhoeae") {
    call ngmaster_task.ngmaster {
      input:
        assembly = assembly,
        samplename = samplename
    }
  }
  if (merlin_tag == "Neisseria meningitidis") {
    call meningotype_task.meningotype {
      input:
        assembly = assembly,
        samplename = samplename
    }
  }
  if (merlin_tag == "Pseudomonas aeruginosa") {
    call pasty_task.pasty {
      input:
        assembly = assembly,
        samplename = samplename,
        min_pident = pasty_min_pident,
        min_coverage = pasty_min_coverage,
        docker = pasty_docker_image
    }
  }
  if (merlin_tag == "Mycobacterium tuberculosis") {
    if (!assembly_only) {
      call tbprofiler_task.tbprofiler { # needs testing
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          ont_data = ont_data,
          tbprofiler_additional_outputs = tbprofiler_additional_outputs,
          output_seq_method_type = output_seq_method_type 
      }
    }
  }
  if (merlin_tag == "Legionella pneumophila") {
    call legsta_task.legsta {
      input:
        assembly = assembly,
        samplename = samplename
    }
  }
  if (merlin_tag == "Staphylococcus aureus") {
      call spatyper_task.spatyper {
        input:
          assembly = assembly,
          samplename = samplename
      }

    call staphopia_sccmec_task.staphopiasccmec {
        input:
          assembly = assembly,
          samplename = samplename,
          docker = staphopia_sccmec_docker_image
      }
    call agrvate_task.agrvate {
        input:
          assembly = assembly,
          samplename = samplename,
          docker = agrvate_docker_image
      }
  }
  if (merlin_tag == "Streptococcus pneumoniae") {
    if (paired_end && !ont_data) {
      call seroba.seroba as seroba_task {
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename
      }
    }
    call pbptyper.pbptyper as pbptyper_task {
      input:
        assembly = assembly,
        samplename = samplename
    }      
    if (call_poppunk) {
      call poppunk_spneumo.poppunk as poppunk_task {
        input:
          assembly = assembly,
          samplename = samplename
      }  
    }
  }
  if (merlin_tag == "Streptococcus pyogenes") {
    if (paired_end && !ont_data) {
      call emmtypingtool_task.emmtypingtool {
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          docker = emmtypingtool_docker_image
      }
    }
  }
  if (merlin_tag == "Haemophilus influenzae") {
    call hicap_task.hicap {
      input:
        assembly = assembly,
        samplename = samplename,
        docker = hicap_docker_image
    }
  }
  if (merlin_tag == "Vibrio") {
    call srst2_vibrio_task.srst2_vibrio {
      input:
        reads1 = select_first([read1]),
        reads2 = read2,
        samplename = samplename,
        srst2_min_cov = srst2_min_cov,
        srst2_max_divergence = srst2_max_divergence,
        srst2_min_depth = srst2_min_depth,
        srst2_min_edge_depth = srst2_min_edge_depth,
        srst2_gene_max_mismatch = srst2_gene_max_mismatch
    }
  }
  
  # theiaeuk
  if (theiaeuk) {
    if (merlin_tag == "Candida auris") {
      call cauris_cladetyper.cauris_cladetyper as cladetyper {
        input: 
          assembly_fasta = assembly,
          samplename = samplename
      }
      if (!assembly_only && !ont_data) {
        call snippy.snippy_variants as snippy_cauris { # no ONT support right now
          input:
            reference_genome_file = cladetyper.clade_spec_ref,
            read1 = select_first([read1]),
            read2 = read2,
            samplename = samplename
        }
        call snippy_gene_query.snippy_gene_query as snippy_gene_query_cauris {
          input:
            samplename = samplename,
            snippy_variants_results = snippy_cauris.snippy_variants_results,
            reference = cladetyper.clade_spec_ref,
            query_gene = select_first([snippy_query_gene,"FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase"]),
        }
      }
    }
    if (merlin_tag == "Candida albicans") {
      if (!assembly_only && !ont_data) {
        call snippy.snippy_variants as snippy_calbicans {
          input:
            reference_genome_file = "gs://theiagen-public-files/terra/theiaeuk_files/Candida_albicans_GCF_000182965.3_ASM18296v3_genomic.gbff",
            read1 = select_first([read1]),
            read2 = read2,
            samplename = samplename
        }
        call snippy_gene_query.snippy_gene_query as snippy_gene_query_calbicans {
          input:
            samplename = samplename,
            snippy_variants_results = snippy_calbicans.snippy_variants_results,
            reference = "gs://theiagen-public-files/terra/theiaeuk_files/Candida_albicans_GCF_000182965.3_ASM18296v3_genomic.gbff",
            query_gene = select_first([snippy_query_gene,"GCS1,ERG11,FUR1,RTA2"]), # GCS1 is another name for FKS1
        }
      }
    }
    if (merlin_tag == "Aspergillus fumigatus") {
      if (!assembly_only && !ont_data) {
        call snippy.snippy_variants as snippy_afumigatus {
          input:
            reference_genome_file = "gs://theiagen-public-files/terra/theiaeuk_files/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.gbff",
            read1 = select_first([read1]),
            read2 = read2,
            samplename = samplename
        }
        call snippy_gene_query.snippy_gene_query as snippy_gene_query_afumigatus {
          input:
            samplename = samplename,
            snippy_variants_results = snippy_afumigatus.snippy_variants_results,
            reference = "gs://theiagen-public-files/terra/theiaeuk_files/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.gbff",
            query_gene = select_first([snippy_query_gene,"Cyp51A,HapE,AFUA_4G08340"]), # AFUA_4G08340 is COX10 according to MARDy
        }
      }
    }
    if (merlin_tag == "Cryptococcus neoformans") {
      if (!assembly_only && !ont_data) {
        call snippy.snippy_variants as snippy_crypto {
          input:
            reference_genome_file = "gs://theiagen-public-files/terra/theiaeuk_files/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.gbff",
            read1 = select_first([read1]),
            read2 = read2,
            samplename = samplename
        }
        call snippy_gene_query.snippy_gene_query as snippy_gene_query_crypto {
          input:
            samplename = samplename,
            snippy_variants_results = snippy_crypto.snippy_variants_results,
            reference = "gs://theiagen-public-files/terra/theiaeuk_files/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.gbff",
            query_gene = select_first([snippy_query_gene,"CNA00300"]), # CNA00300 is ERG11 for this reference genome
        }
      }
    }
  }
  output {
    # theiaprok
    # Ecoli Typing
    File? serotypefinder_report = serotypefinder.serotypefinder_report
    String? serotypefinder_docker = serotypefinder.serotypefinder_docker
    String? serotypefinder_serotype = serotypefinder.serotypefinder_serotype
    File? ectyper_results = ectyper.ectyper_results
    String? ectyper_version = ectyper.ectyper_version
    String? ectyper_predicted_serotype = ectyper.ectyper_predicted_serotype
    String? shigatyper_predicted_serotype = shigatyper.shigatyper_predicted_serotype
    String? shigatyper_ipaB_presence_absence = shigatyper.shigatyper_ipaB_presence_absence
    String? shigatyper_notes = shigatyper.shigatyper_notes
    File? shigatyper_hits_tsv = shigatyper.shigatyper_hits_tsv
    File? shigatyper_summary_tsv = shigatyper.shigatyper_summary_tsv
    String? shigatyper_version = shigatyper.shigatyper_version
    String? shigatyper_docker = shigatyper.shigatyper_docker
    File? shigeifinder_report = shigeifinder.shigeifinder_report
    String? shigeifinder_docker = shigeifinder.shigeifinder_docker
    String? shigeifinder_version = shigeifinder.shigeifinder_version
    String? shigeifinder_ipaH_presence_absence = shigeifinder.shigeifinder_ipaH_presence_absence
    String? shigeifinder_num_virulence_plasmid_genes = shigeifinder.shigeifinder_num_virulence_plasmid_genes
    String? shigeifinder_cluster = shigeifinder.shigeifinder_cluster
    String? shigeifinder_serotype = shigeifinder.shigeifinder_serotype
    String? shigeifinder_O_antigen = shigeifinder.shigeifinder_O_antigen
    String? shigeifinder_H_antigen = shigeifinder.shigeifinder_H_antigen
    String? shigeifinder_notes = shigeifinder.shigeifinder_notes
    # ShigeiFinder outputs but for task that uses reads instead of assembly as input
    File? shigeifinder_report_reads = shigeifinder_reads.shigeifinder_report
    String? shigeifinder_docker_reads = shigeifinder_reads.shigeifinder_docker
    String? shigeifinder_version_reads = shigeifinder_reads.shigeifinder_version
    String? shigeifinder_ipaH_presence_absence_reads = shigeifinder_reads.shigeifinder_ipaH_presence_absence
    String? shigeifinder_num_virulence_plasmid_genes_reads = shigeifinder_reads.shigeifinder_num_virulence_plasmid_genes
    String? shigeifinder_cluster_reads = shigeifinder_reads.shigeifinder_cluster
    String? shigeifinder_serotype_reads = shigeifinder_reads.shigeifinder_serotype
    String? shigeifinder_O_antigen_reads = shigeifinder_reads.shigeifinder_O_antigen
    String? shigeifinder_H_antigen_reads = shigeifinder_reads.shigeifinder_H_antigen
    String? shigeifinder_notes_reads = shigeifinder_reads.shigeifinder_notes
    # E coli only typing
    File? virulencefinder_report_tsv = virulencefinder.virulencefinder_report_tsv
    String? virulencefinder_docker = virulencefinder.virulencefinder_docker
    String? virulencefinder_hits = virulencefinder.virulencefinder_hits
    # Shigella sonnei Typing
    File? sonneityping_mykrobe_report_csv = sonneityping.sonneityping_mykrobe_report_csv
    File? sonneityping_mykrobe_report_json = sonneityping.sonneityping_mykrobe_report_json
    File? sonneityping_final_report_tsv = sonneityping.sonneityping_final_report_tsv
    String? sonneityping_mykrobe_version = sonneityping.sonneityping_mykrobe_version
    String? sonneityping_mykrobe_docker = sonneityping.sonneityping_mykrobe_docker
    String? sonneityping_species = sonneityping.sonneityping_species
    String? sonneityping_final_genotype = sonneityping.sonneityping_final_genotype
    String? sonneityping_genotype_confidence = sonneityping.sonneityping_genotype_confidence
    String? sonneityping_genotype_name = sonneityping.sonneityping_genotype_name
    # Listeria Typing
    File? lissero_results = lissero.lissero_results
    String? lissero_version = lissero.lissero_version
    String? lissero_serotype = lissero.lissero_serotype
    # Pseudomonas Aeruginosa Typing
    String? pasty_serogroup = pasty.pasty_serogroup
    Float? pasty_serogroup_coverage = pasty.pasty_serogroup_coverage
    Int? pasty_serogroup_fragments = pasty.pasty_serogroup_fragments
    File? pasty_summary_tsv = pasty.pasty_summary_tsv
    File? pasty_blast_hits = pasty.pasty_blast_hits
    File? pasty_all_serogroups = pasty.pasty_all_serogroups
    String? pasty_version = pasty.pasty_version
    String? pasty_docker = pasty.pasty_docker
    String? pasty_comment = pasty.pasty_comment
    # Salmonella Typing
    File? sistr_results = sistr.sistr_results
    File? sistr_allele_json = sistr.sistr_allele_json
    File? sistr_allele_fasta = sistr.sistr_allele_fasta
    File? sistr_cgmlst = sistr.sistr_cgmlst
    String? sistr_version = sistr.sistr_version
    String? sistr_predicted_serotype = sistr.sistr_predicted_serotype
    String seqsero2_report = select_first([seqsero2.seqsero2_report, seqsero2_assembly.seqsero2_report, ""])
    String seqsero2_version = select_first([seqsero2.seqsero2_version, seqsero2_assembly.seqsero2_version, ""])
    String seqsero2_predicted_antigenic_profile = select_first([seqsero2.seqsero2_predicted_antigenic_profile, seqsero2_assembly.seqsero2_predicted_antigenic_profile, ""])
    String seqsero2_predicted_serotype = select_first([seqsero2.seqsero2_predicted_serotype, seqsero2_assembly.seqsero2_predicted_serotype, ""])
    String? seqsero2_predicted_contamination = seqsero2.seqsero2_predicted_contamination
    # Salmonella serotype Typhi typing
    File? genotyphi_report_tsv = genotyphi_task.genotyphi_report_tsv 
    File? genotyphi_mykrobe_json = genotyphi_task.genotyphi_mykrobe_json
    String? genotyphi_version = genotyphi_task.genotyphi_version
    String? genotyphi_species = genotyphi_task.genotyphi_species
    Float? genotyphi_st_probes_percent_coverage = genotyphi_task.genotyphi_st_probes_percent_coverage
    String? genotyphi_final_genotype = genotyphi_task.genotyphi_final_genotype
    String? genotyphi_genotype_confidence = genotyphi_task.genotyphi_genotype_confidence
    # Klebsiella Typing
    File? kleborate_output_file = kleborate.kleborate_output_file
    String? kleborate_version = kleborate.kleborate_version
    String? kleborate_docker = kleborate.kleborate_docker
    String? kleborate_key_resistance_genes = kleborate.kleborate_key_resistance_genes
    String? kleborate_genomic_resistance_mutations = kleborate.kleborate_genomic_resistance_mutations
    String? kleborate_mlst_sequence_type = kleborate.kleborate_mlst_sequence_type
    String? kleborate_klocus = kleborate.kleborate_klocus
    String? kleborate_ktype = kleborate.kleborate_ktype
    String? kleborate_olocus = kleborate.kleborate_olocus
    String? kleborate_otype = kleborate.kleborate_otype
    String? kleborate_klocus_confidence = kleborate.kleborate_klocus_confidence
    String? kleborate_olocus_confidence = kleborate.kleborate_olocus_confidence
    # Neisseria gonorrhoeae Typing
    File? ngmaster_tsv = ngmaster.ngmaster_tsv
    String? ngmaster_version = ngmaster.ngmaster_version
    String? ngmaster_ngmast_sequence_type = ngmaster.ngmaster_ngmast_sequence_type
    String? ngmaster_ngmast_porB_allele = ngmaster.ngmaster_ngmast_porB_allele
    String? ngmaster_ngmast_tbpB_allele = ngmaster.ngmaster_ngmast_tbpB_allele
    String? ngmaster_ngstar_sequence_type = ngmaster.ngmaster_ngstar_sequence_type
    String? ngmaster_ngstar_penA_allele = ngmaster.ngmaster_ngstar_penA_allele
    String? ngmaster_ngstar_mtrR_allele = ngmaster.ngmaster_ngstar_mtrR_allele
    String? ngmaster_ngstar_porB_allele = ngmaster.ngmaster_ngstar_porB_allele
    String? ngmaster_ngstar_ponA_allele = ngmaster.ngmaster_ngstar_ponA_allele
    String? ngmaster_ngstar_gyrA_allele = ngmaster.ngmaster_ngstar_gyrA_allele
    String? ngmaster_ngstar_parC_allele = ngmaster.ngmaster_ngstar_parC_allele
    String? ngmaster_ngstar_23S_allele = ngmaster.ngmaster_ngstar_23S_allele
    # Neisseria meningitidis Typing
    File? meningotype_tsv = meningotype.meningotype_tsv
    String? meningotype_version = meningotype.meningotype_version
    String? meningotype_serogroup = meningotype.meningotype_serogroup
    String? meningotype_PorA = meningotype.meningotype_PorA
    String? meningotype_FetA = meningotype.meningotype_FetA
    String? meningotype_PorB = meningotype.meningotype_PorB
    String? meningotype_fHbp = meningotype.meningotype_fHbp
    String? meningotype_NHBA = meningotype.meningotype_NHBA
    String? meningotype_NadA = meningotype.meningotype_NadA
    String? meningotype_BAST = meningotype.meningotype_BAST
    # Acinetobacter Typing
    File? kaptive_output_file_k = kaptive.kaptive_output_file_k
    File? kaptive_output_file_oc = kaptive.kaptive_output_file_oc
    String? kaptive_version = kaptive.kaptive_version
    String? kaptive_k_match = kaptive.kaptive_k_match
    String? kaptive_k_type = kaptive.kaptive_k_type
    String? kaptive_k_confidence = kaptive.kaptive_k_confidence
    String? kaptive_oc_match = kaptive.kaptive_oc_match
    String? kaptive_oc_type = kaptive.kaptive_oc_type
    String? kaptive_oc_confidence = kaptive.kaptive_oc_confidence
    File? abricate_results = abricate.abricate_results
    String? abricate_genes = abricate.abricate_genes
    String? abricate_database = abricate.abricate_database
    String? abricate_version = abricate.abricate_version
    String? abricate_docker = abricate.abricate_docker
    # Mycobacterium Typing
    File? tbprofiler_output_file = tbprofiler.tbprofiler_output_csv
    File? tbprofiler_output_bam = tbprofiler.tbprofiler_output_bam
    File? tbprofiler_output_bai = tbprofiler.tbprofiler_output_bai
    String? tbprofiler_version = tbprofiler.version
    String? tbprofiler_main_lineage = tbprofiler.tbprofiler_main_lineage
    String? tbprofiler_sub_lineage = tbprofiler.tbprofiler_sub_lineage
    String? tbprofiler_dr_type = tbprofiler.tbprofiler_dr_type
    String? tbprofiler_resistance_genes = tbprofiler.tbprofiler_resistance_genes
    File? tbprofiler_additional_outputs_csv = tbprofiler.tbprofiler_additional_outputs_csv
    File? tbprofiler_laboratorian_report_csv = tbprofiler.tbprofiler_laboratorian_report_csv
    String? tbprofiler_gene_name = tbprofiler.tbprofiler_gene_name
    String? tbprofiler_locus_tag = tbprofiler.tbprofiler_locus_tag
    String? tbprofiler_variant_substitutions = tbprofiler.tbprofiler_variant_substitutions
    String? tbprofiler_output_seq_method_type = tbprofiler.tbprofiler_output_seq_method_type
    # Legionella pneumophila Typing
    File? legsta_results = legsta.legsta_results
    String? legsta_predicted_sbt = legsta.legsta_predicted_sbt
    String? legsta_version = legsta.legsta_version
    # Staphylococcus aureus
    File? spatyper_tsv = spatyper.spatyper_tsv
    String? spatyper_docker = spatyper.spatyper_docker
    String? spatyper_repeats = spatyper.spatyper_repeats
    String? spatyper_type = spatyper.spatyper_type
    String? spatyper_version = spatyper.spatyper_version
    File? staphopiasccmec_results_tsv = staphopiasccmec.staphopiasccmec_results_tsv
    File? staphopiasccmec_hamming_distance_tsv = staphopiasccmec.staphopiasccmec_hamming_distance_tsv
    String? staphopiasccmec_types_and_mecA_presence = staphopiasccmec.staphopiasccmec_types_and_mecA_presence
    String? staphopiasccmec_version = staphopiasccmec.staphopiasccmec_version
    String? staphopiasccmec_docker = staphopiasccmec.staphopiasccmec_docker
    File? agrvate_summary = agrvate.agrvate_summary
    File? agrvate_results = agrvate.agrvate_results
    String? agrvate_agr_group = agrvate.agrvate_agr_group
    String? agrvate_agr_match_score = agrvate.agrvate_agr_match_score
    String? agrvate_agr_canonical = agrvate.agrvate_agr_canonical
    String? agrvate_agr_multiple = agrvate.agrvate_agr_multiple
    String? agrvate_agr_num_frameshifts = agrvate.agrvate_agr_num_frameshifts
    String? agrvate_version = agrvate.agrvate_version
    String? agrvate_docker = agrvate.agrvate_docker
    # Streptococcus pneumoniae Typing
    String? pbptyper_predicted_1A_2B_2X = pbptyper_task.pbptyper_predicted_1A_2B_2X
    File? pbptyper_pbptype_predicted_tsv = pbptyper_task.pbptyper_pbptype_predicted_tsv
    File? pbptyper_pbptype_1A_tsv = pbptyper_task.pbptyper_pbptype_1A_tsv
    File? pbptyper_pbptype_2B_tsv = pbptyper_task.pbptyper_pbptype_2B_tsv
    File? pbptyper_pbptype_2X_tsv = pbptyper_task.pbptyper_pbptype_2X_tsv
    String? pbptyper_version = pbptyper_task.pbptyper_version
    String? pbptyper_docker = pbptyper_task.pbptyper_docker
    String? poppunk_gps_cluster = poppunk_task.poppunk_gps_cluster
    File? poppunk_gps_external_cluster_csv = poppunk_task.poppunk_gps_external_cluster_csv
    String? poppunk_GPS_db_version = poppunk_task.poppunk_GPS_db_version
    String? poppunk_version = poppunk_task.poppunk_version
    String? poppunk_docker = poppunk_task.poppunk_docker
    String? seroba_version = seroba_task.seroba_version
    String? seroba_docker = seroba_task.seroba_docker
    String? seroba_serotype = seroba_task.seroba_serotype
    String? seroba_ariba_serotype = seroba_task.seroba_ariba_serotype
    String? seroba_ariba_identity = seroba_task.seroba_ariba_identity
    File? seroba_details = seroba_task.seroba_details
    # Streptococcus pyogenes Typing
    String? emmtypingtool_emm_type = emmtypingtool.emmtypingtool_emm_type
    File? emmtypingtool_results_xml = emmtypingtool.emmtypingtool_results_xml
    String? emmtypingtool_version = emmtypingtool.emmtypingtool_version
    String? emmtypingtool_docker = emmtypingtool.emmtypingtool_docker
    # Haemophilus influenzae Typing
    String? hicap_serotype = hicap.hicap_serotype
    String? hicap_genes = hicap.hicap_genes
    File? hicap_results_tsv = hicap.hicap_results_tsv
    String? hicap_version = hicap.hicap_version
    String? hicap_docker = hicap.hicap_docker
    # Vibrio
    File? srst2_vibrio_detailed_tsv = srst2_vibrio.srst2_detailed_tsv
    String? srst2_vibrio_version = srst2_vibrio.srst2_version
    String? srst2_vibrio_ctxA = srst2_vibrio.srst2_vibrio_ctxA
    String? srst2_vibrio_ompW = srst2_vibrio.srst2_vibrio_ompW
    String? srst2_vibrio_toxR = srst2_vibrio.srst2_vibrio_toxR
    String? srst2_vibrio_serogroup = srst2_vibrio.srst2_vibrio_serogroup
    String? srst2_vibrio_biotype = srst2_vibrio.srst2_vibrio_biotype
    # theiaeuk
    # c auris 
    String? clade_type = cladetyper.gambit_cladetype
    String? cladetyper_analysis_date = cladetyper.date
    String? cladetyper_version = cladetyper.version
    String? cladetyper_docker_image = cladetyper.gambit_cladetyper_docker_image
    String? cladetype_annotated_ref = cladetyper.clade_spec_ref
    # snippy variants
    String snippy_variants_version = select_first([snippy_cauris.snippy_variants_version, snippy_calbicans.snippy_variants_version, snippy_afumigatus.snippy_variants_version, snippy_crypto.snippy_variants_version, "No matching taxon detected"])
    String snippy_variants_query = select_first([snippy_gene_query_cauris.snippy_variants_query, snippy_gene_query_calbicans.snippy_variants_query, snippy_gene_query_afumigatus.snippy_variants_query, snippy_gene_query_crypto.snippy_variants_query, "No matching taxon detected"])
    String snippy_variants_query_check = select_first([snippy_gene_query_cauris.snippy_variants_query_check, snippy_gene_query_calbicans.snippy_variants_query_check, snippy_gene_query_afumigatus.snippy_variants_query_check, snippy_gene_query_crypto.snippy_variants_query_check, "No matching taxon detected"])
    String snippy_variants_hits = select_first([snippy_gene_query_cauris.snippy_variants_hits, snippy_gene_query_calbicans.snippy_variants_hits, snippy_gene_query_afumigatus.snippy_variants_hits, snippy_gene_query_crypto.snippy_variants_hits, "No matching taxon detected"])
    File snippy_variants_gene_query_results = select_first([snippy_gene_query_cauris.snippy_variants_gene_query_results, snippy_gene_query_calbicans.snippy_variants_gene_query_results, snippy_gene_query_afumigatus.snippy_variants_gene_query_results, snippy_gene_query_crypto.snippy_variants_gene_query_results, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    File snippy_variants_outdir_tarball = select_first([snippy_cauris.snippy_variants_outdir_tarball, snippy_calbicans.snippy_variants_outdir_tarball, snippy_afumigatus.snippy_variants_outdir_tarball, snippy_crypto.snippy_variants_outdir_tarball, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    File snippy_variants_results = select_first([snippy_cauris.snippy_variants_results, snippy_calbicans.snippy_variants_results, snippy_afumigatus.snippy_variants_results, snippy_crypto.snippy_variants_results, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    File snippy_variants_bam = select_first([snippy_cauris.snippy_variants_bam, snippy_calbicans.snippy_variants_bam, snippy_afumigatus.snippy_variants_bam, snippy_crypto.snippy_variants_bam, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    File snippy_variants_bai = select_first([snippy_cauris.snippy_variants_bai, snippy_calbicans.snippy_variants_bai, snippy_afumigatus.snippy_variants_bai, snippy_crypto.snippy_variants_bai, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    File snippy_variants_summary = select_first([snippy_cauris.snippy_variants_summary, snippy_calbicans.snippy_variants_summary, snippy_afumigatus.snippy_variants_summary, snippy_crypto.snippy_variants_summary, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
  }
}