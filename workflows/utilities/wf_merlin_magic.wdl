version 1.0

# theiaprok
import "../../tasks/gene_typing/drug_resistance/task_abricate.wdl" as abricate_task
import "../../tasks/species_typing/acinetobacter/task_kaptive.wdl" as kaptive_task
import "../../tasks/species_typing/escherichia_shigella/task_ectyper.wdl" as ectyper_task
import "../../tasks/species_typing/escherichia_shigella/task_serotypefinder.wdl" as serotypefinder_task
import "../../tasks/species_typing/escherichia_shigella/task_shigatyper.wdl" as shigatyper_task
import "../../tasks/species_typing/escherichia_shigella/task_shigeifinder.wdl" as shigeifinder_task
import "../../tasks/species_typing/escherichia_shigella/task_sonneityping.wdl" as sonneityping_task
import "../../tasks/species_typing/escherichia_shigella/task_virulencefinder.wdl" as virulencefinder_task
import "../../tasks/species_typing/haemophilus/task_hicap.wdl" as hicap_task
import "../../tasks/species_typing/klebsiella/task_kleborate.wdl" as kleborate_task
import "../../tasks/species_typing/legionella/task_legsta.wdl" as legsta_task
import "../../tasks/species_typing/listeria/task_lissero.wdl" as lissero_task
import "../../tasks/species_typing/mycobacterium/task_clockwork.wdl" as clockwork_task
import "../../tasks/species_typing/mycobacterium/task_tbp_parser.wdl" as tbp_parser_task
import "../../tasks/species_typing/mycobacterium/task_tbprofiler.wdl" as tbprofiler_task
import "../../tasks/species_typing/neisseria/task_meningotype.wdl" as meningotype_task
import "../../tasks/species_typing/neisseria/task_ngmaster.wdl" as ngmaster_task
import "../../tasks/species_typing/pseudomonas/task_pasty.wdl" as pasty_task
import "../../tasks/species_typing/salmonella/task_genotyphi.wdl" as genotyphi
import "../../tasks/species_typing/salmonella/task_seqsero2.wdl" as seqsero2_task
import "../../tasks/species_typing/salmonella/task_sistr.wdl" as sistr_task
import "../../tasks/species_typing/staphylococcus/task_agrvate.wdl" as agrvate_task
import "../../tasks/species_typing/staphylococcus/task_spatyper.wdl" as spatyper_task
import "../../tasks/species_typing/staphylococcus/task_staphopiasccmec.wdl" as staphopia_sccmec_task
import "../../tasks/species_typing/streptococcus/task_emmtyper.wdl" as emmtyper_task
import "../../tasks/species_typing/streptococcus/task_emmtypingtool.wdl" as emmtypingtool_task
import "../../tasks/species_typing/streptococcus/task_pbptyper.wdl" as pbptyper
import "../../tasks/species_typing/streptococcus/task_poppunk_streppneumo.wdl" as poppunk_spneumo
import "../../tasks/species_typing/streptococcus/task_seroba.wdl" as seroba
import "../../tasks/species_typing/vibrio/task_srst2_vibrio.wdl" as srst2_vibrio_task
import "../../tasks/species_typing/vibrio/task_abricate_vibrio.wdl" as abricate_vibrio_task

# theiaeuk
import "../../tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl" as snippy_gene_query
import "../../tasks/gene_typing/variant_detection/task_snippy_variants.wdl" as snippy
import "../../tasks/species_typing/candida/task_cauris_cladetyper.wdl" as cauris_cladetyper

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
    # subworkflow logic
    Boolean assembly_only = false
    Boolean ont_data = false
    Boolean paired_end = true
    Boolean theiaeuk = false
    # activating tool logic
    Boolean call_poppunk = true
    Boolean call_shigeifinder_reads_input = false
    Boolean tbprofiler_additional_outputs = false # set to true to run tbp-parser
    # docker options
    String? abricate_abaum_docker_image
    String? abricate_vibrio_docker_image
    String? agrvate_docker_image
    String? cauris_cladetyper_docker_image
    String? clockwork_docker_image
    String? ectyper_docker_image
    String? emmtyper_docker_image
    String? emmtypingtool_docker_image
    String? genotyphi_docker_image
    String? hicap_docker_image
    String? kaptive_docker_image
    String? kleborate_docker_image
    String? legsta_docker_image
    String? lissero_docker_image
    String? meningotype_docker_image
    String? ngmaster_docker_image
    String? pasty_docker_image
    String? pbptyper_docker_image
    String? poppunk_docker_image
    String? seqsero2_docker_image
    String? seroba_docker_image
    String? serotypefinder_docker_image
    String? shigatyper_docker_image
    String? shigeifinder_docker_image
    String? sistr_docker_image
    String? snippy_gene_query_docker_image
    String? snippy_variants_docker_image
    String? sonneityping_docker_image
    String? spatyper_docker_image
    String? srst2_docker_image
    String? staphopia_sccmec_docker_image
    String? tbprofiler_docker_image
    String? tbp_parser_docker_image
    String? virulencefinder_docker_image
    # abricate abaum options
    Int abricate_abaum_minid = 95 # strict threshold of 95% identity for typing purposes
    Int? abricate_abaum_mincov
    # abricate vibrio options
    Int abricate_vibrio_minid = 80
    Int abricate_vibrio_mincov = 80
    # agrvate options
    Boolean? agrvate_agr_typing_only
    # cladetyper options - primarily files we host
    Int? cladetyper_kmer_size
    File? cladetyper_ref_clade1
    File? cladetyper_ref_clade1_annotated
    File? cladetyper_ref_clade2
    File? cladetyper_ref_clade2_annotated
    File? cladetyper_ref_clade3
    File? cladetyper_ref_clade3_annotated
    File? cladetyper_ref_clade4
    File? cladetyper_ref_clade4_annotated
    File? cladetyper_ref_clade5
    File? cladetyper_ref_clade5_annotated
    # ectyper options
    Int? ectyper_opid
    Int? ectyper_hpid
    Int? ectyper_opcov
    Int? ectyper_hpcov
    Boolean? ectyper_verify
    Boolean? ectyper_print_alleles
    # emmtyper options
    String? emmtyper_wf
    Int? emmtyper_cluster_distance
    Int? emmtyper_percid
    Int? emmtyper_culling_limit
    Int? emmtyper_mismatch
    Int? emmtyper_align_diff
    Int? emmtyper_gap
    Int? emmtyper_min_perfect
    Int? emmtyper_min_good
    Int? emmtyper_max_size
    # kaptive options
    Int? kaptive_start_end_margin
    Float? kaptive_min_identity
    Float? kaptive_min_coverage
    Float? kaptive_low_gene_id
    # kleborate options
    Boolean? kleborate_skip_resistance
    Boolean? kleborate_skip_kaptive
    Float? kleborate_min_identity
    Float? kleborate_min_coverage
    Float? kleborate_min_spurious_identity
    Float? kleborate_min_spurious_coverage
    String? kleborate_min_kaptive_confidence
    # lissero options
    Float? lissero_min_id
    Float? lissero_min_cov
    # pasty options
    Int? pasty_min_pident
    Int? pasty_min_coverage      
    # pbptyper options 
    Int? pbptyper_min_pident
    Int? pbptyper_min_coverage
    # popppunk options - primarily files we host
    File? poppunk_gps_dists_npy
    File? poppunk_gps_dists_pkl
    File? poppunk_gps_h5
    File? poppunk_gps_refs
    File? poppunk_gps_refs_dists_npy
    File? poppunk_gps_refs_dists_pkl
    File? poppunk_gps_refs_h5
    File? poppunk_gps_clusters_csv
    File? poppunk_gps_fit_npz
    File? poppunk_gps_fit_pkl
    File? poppunk_gps_graph_gt
    File? poppunk_gps_qcreport_txt
    File? poppunk_gps_unword_clusters_csv
    File? poppunk_gps_refs_graph_gt
    File? poppunk_gps_external_clusters_csv
    # sistr options
    Boolean? sistr_use_full_cgmlst_db
    Int? sistr_cpu
    Int? sistr_memory
    Int? sistr_disk_size
    # snippy options - mostly files we host
    String? snippy_query_gene
    File snippy_reference_afumigatus = "gs://theiagen-public-files/terra/theiaeuk_files/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.gbff"
    File snippy_reference_cryptoneo = "gs://theiagen-public-files/terra/theiaeuk_files/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.gbff"
    Int? snippy_map_qual
    Int? snippy_base_quality
    Int? snippy_min_coverage
    Float? snippy_min_frac
    Int? snippy_min_quality
    Int? snippy_maxsoft
    #File snippy_reference_calbicans = "gs://theiagen-public-files/terra/theiaeuk_files/Candida_albicans_GCF_000182965.3_ASM18296v3_genomic.gbff"
    # sonneityping options
    String? sonneityping_mykrobe_opts
    # spatyper options
    Boolean? spatyper_do_enrich
    # srst2 options
    Int srst2_min_cov = 80
    Int srst2_max_divergence = 20
    Int srst2_min_depth = 5
    Int srst2_min_edge_depth = 2
    Int srst2_gene_max_mismatch = 2000
    # tbprofiler options
    Boolean tbprofiler_run_custom_db = false
    File? tbprofiler_custom_db
    Int? tbprofiler_cov_frac_threshold
    Float? tbprofiler_min_af
    Float? tbprofiler_min_af_pred
    Int? tbprofiler_min_depth
    String? tbprofiler_mapper
    String? tbprofiler_variant_caller
    String? tbprofiler_variant_calling_params
    # tbp-parser options
    String tbp_parser_output_seq_method_type = "WGS"
    String? tbp_parser_operator
    Int? tbp_parser_min_depth
    Int? tbp_parser_min_frequency
    Int? tbp_parser_min_read_support
    Int? tbp_parser_coverage_threshold
    File? tbp_parser_coverage_regions_bed
    Boolean? tbp_parser_debug
    Boolean? tbp_parser_add_cs_lims
    # virulencefinder options
    Float? virulencefinder_coverage_threshold
    Float? virulencefinder_identity_threshold
    String? virulencefinder_database
  }
  # theiaprok
  if (merlin_tag == "Acinetobacter baumannii") {
    call kaptive_task.kaptive {
      input:
        assembly = assembly,
        samplename = samplename,
        start_end_margin = kaptive_start_end_margin,
        min_identity = kaptive_min_identity,
        min_coverage = kaptive_min_coverage,
        low_gene_id = kaptive_low_gene_id,
        docker = kaptive_docker_image
    }
    call abricate_task.abricate as abricate_abaum {
      input:
        assembly = assembly,
        samplename = samplename,
        database = "AcinetobacterPlasmidTyping",
        minid = abricate_abaum_minid, 
        mincov = abricate_abaum_mincov,
        docker = abricate_abaum_docker_image
    }
  }
  if (merlin_tag == "Escherichia" || merlin_tag == "Shigella sonnei" ) {
    # tools specific to ALL Escherichia and Shigella species
    #
    # FYI see the GAMBIT task for all merlin_tag designations but all Escherichia and Shigella species are given "Escherichia" merlin_tag designation, except for Shigella sonnei which is given "Shigella sonnei" merlin_tag.
    # The reason being is that S. sonnei does have a species-specific tool (sonneityping) where the other Shigella species do not (flexneri, dysenteriae, boydii, etc.)
    call serotypefinder_task.serotypefinder {
      input:
        assembly = assembly,
        samplename = samplename,
        docker = serotypefinder_docker_image
    }
    call ectyper_task.ectyper {
      input:
        assembly = assembly,
        samplename = samplename,
        opid = ectyper_opid,
        hpid = ectyper_hpid,
        opcov = ectyper_opcov,
        hpcov = ectyper_hpcov,
        verify = ectyper_verify,
        print_alleles = ectyper_print_alleles,
        docker = ectyper_docker_image
    }
    if (!assembly_only) {
      call shigatyper_task.shigatyper {
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          read1_is_ont = ont_data,
          docker = shigatyper_docker_image
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
          paired_end = paired_end,
          docker = shigeifinder_docker_image
      }
    }
    call virulencefinder_task.virulencefinder {
      input:
      #  read1 = read1,
      #  read2 = read2,
        assembly = assembly,
        samplename = samplename,
      #  paired_end = paired_end,
      #  assembly_only = assembly_only,
      #  ont_data = ont_data,
        coverage_threshold = virulencefinder_coverage_threshold,
        identity_threshold = virulencefinder_identity_threshold,
        database = virulencefinder_database,
        docker = virulencefinder_docker_image
    }
  }
  if (merlin_tag == "Shigella sonnei") {
    if (!assembly_only) {
      call sonneityping_task.sonneityping { 
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          ont_data = ont_data,
          mykrobe_opts = sonneityping_mykrobe_opts,
          docker = sonneityping_docker_image
      }
    }
  }
  if (merlin_tag == "Listeria") {
    call lissero_task.lissero {
      input:
        assembly = assembly,
        samplename = samplename,
        min_id = lissero_min_id,
        min_cov = lissero_min_cov,
        docker = lissero_docker_image
    }
  }
  if (merlin_tag == "Salmonella") {
    call sistr_task.sistr {
      input: 
        assembly = assembly,
        samplename = samplename,
        use_full_cgmlst_db = sistr_use_full_cgmlst_db,
        docker = sistr_docker_image,
        cpu = sistr_cpu,
        memory = sistr_memory,
        disk_size = sistr_disk_size
    }
    if (!ont_data && !assembly_only) {
      call seqsero2_task.seqsero2 { 
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          paired_end = paired_end,
          docker = seqsero2_docker_image
      }
    }
    if (ont_data || assembly_only) {
      call seqsero2_task.seqsero2_assembly {
        input:
          assembly_fasta = assembly,
          samplename = samplename,
          docker = seqsero2_docker_image
      }
    }
    if ((select_first([seqsero2.seqsero2_predicted_serotype, seqsero2_assembly.seqsero2_predicted_serotype]) == "Typhi" || sistr.sistr_predicted_serotype == "Typhi") && !assembly_only) {
      call genotyphi.genotyphi as genotyphi_task {
        input: 
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          ont_data = ont_data,
          docker = genotyphi_docker_image
      }
    }
  }
  # see here for appropriate Klebsiella species & subspecies to be analyzed w Kleborate: https://github.com/klebgenomics/Kleborate/wiki/Species-detection
  if (merlin_tag == "Klebsiella" || merlin_tag == "Klebsiella pneumoniae" || merlin_tag == "Klebsiella variicola" || merlin_tag == "Klebsiella aerogenes" || merlin_tag == "Klebsiella oxytoca") {
    call kleborate_task.kleborate {
      input:
        assembly = assembly,
        samplename = samplename,
        skip_resistance = kleborate_skip_resistance,
        skip_kaptive = kleborate_skip_kaptive,
        min_identity = kleborate_min_identity,
        min_coverage = kleborate_min_coverage,
        min_spurious_identity = kleborate_min_spurious_identity,
        min_spurious_coverage = kleborate_min_spurious_coverage,
        min_kaptive_confidence = kleborate_min_kaptive_confidence,
        docker = kleborate_docker_image
    }
  }
  if (merlin_tag == "Neisseria gonorrhoeae") {
    call ngmaster_task.ngmaster {
      input:
        assembly = assembly,
        samplename = samplename,
        docker = ngmaster_docker_image
    }
  }
  if (merlin_tag == "Neisseria meningitidis") {
    call meningotype_task.meningotype {
      input:
        assembly = assembly,
        samplename = samplename,
        docker = meningotype_docker_image
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
      if (paired_end && !ont_data) {
        call clockwork_task.clockwork_decon_reads {
          input:
            read1 = select_first([read1]),
            read2 = read2,
            samplename = samplename,
            docker = clockwork_docker_image
        }
      }
      call tbprofiler_task.tbprofiler {
        input:
          read1 = select_first([clockwork_decon_reads.clockwork_cleaned_read1, read1]),
          read2 = select_first([clockwork_decon_reads.clockwork_cleaned_read2, read2, "gs://theiagen-public-files/terra/theiaprok-files/no-read2.txt"]),
          samplename = samplename,
          ont_data = ont_data,
          tbprofiler_run_custom_db = tbprofiler_run_custom_db,
          tbprofiler_custom_db = tbprofiler_custom_db,
          cov_frac_threshold = tbprofiler_cov_frac_threshold,
          min_af = tbprofiler_min_af,
          min_af_pred = tbprofiler_min_af_pred,
          min_depth = tbprofiler_min_depth,
          mapper = tbprofiler_mapper,
          variant_caller = tbprofiler_variant_caller,
          variant_calling_params = tbprofiler_variant_calling_params,
          docker = tbprofiler_docker_image
      }
      if (tbprofiler_additional_outputs) {
        call tbp_parser_task.tbp_parser {
          input:
            tbprofiler_json = tbprofiler.tbprofiler_output_json,
            tbprofiler_bam = tbprofiler.tbprofiler_output_bam,
            tbprofiler_bai = tbprofiler.tbprofiler_output_bai,
            samplename = samplename, 
            sequencing_method = tbp_parser_output_seq_method_type,
            operator = tbp_parser_operator,
            coverage_threshold = tbp_parser_coverage_threshold,
            coverage_regions_bed = tbp_parser_coverage_regions_bed,
            min_depth = tbp_parser_min_depth,
            min_frequency = tbp_parser_min_frequency,
            min_read_support = tbp_parser_min_read_support,
            tbp_parser_debug = tbp_parser_debug,
            add_cycloserine_lims = tbp_parser_add_cs_lims,
            docker = tbp_parser_docker_image
        }
      }
    }
  }
  if (merlin_tag == "Legionella pneumophila") {
    call legsta_task.legsta {
      input:
        assembly = assembly,
        samplename = samplename,
        docker = legsta_docker_image
    }
  }
  if (merlin_tag == "Staphylococcus aureus") {
    call spatyper_task.spatyper {
      input:
        assembly = assembly,
        samplename = samplename,
        do_enrich = spatyper_do_enrich,
        docker = spatyper_docker_image
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
        typing_only = agrvate_agr_typing_only,
        docker = agrvate_docker_image
     }
  }
  if (merlin_tag == "Streptococcus pneumoniae") {
    if (paired_end && !ont_data) {
      call seroba.seroba as seroba_task {
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          docker = seroba_docker_image
      }
    }
    call pbptyper.pbptyper as pbptyper_task {
      input:
        assembly = assembly,
        samplename = samplename,
        min_pident = pbptyper_min_pident,
        min_coverage = pbptyper_min_coverage,
        docker = pbptyper_docker_image
    }      
    if (call_poppunk) {
      call poppunk_spneumo.poppunk as poppunk_task {
        input:
          assembly = assembly,
          samplename = samplename,
          GPS_dists_npy = poppunk_gps_dists_npy,
          GPS_dists_pkl = poppunk_gps_dists_pkl,
          GPS_h5 = poppunk_gps_h5,
          GPS_refs = poppunk_gps_refs,
          GPS_refs_dists_npy = poppunk_gps_refs_dists_npy,
          GPS_refs_dists_pkl = poppunk_gps_refs_dists_pkl,
          GPS_refs_h5 = poppunk_gps_refs_h5,
          GPS_clusters_csv = poppunk_gps_clusters_csv,
          GPS_fit_npz = poppunk_gps_fit_npz,
          GPS_fit_pkl = poppunk_gps_fit_pkl,
          GPS_graph_gt = poppunk_gps_graph_gt,
          GPS_qcreport_txt = poppunk_gps_qcreport_txt,
          GPS_unword_clusters_csv = poppunk_gps_unword_clusters_csv,
          GPS_refs_graph_gt = poppunk_gps_refs_graph_gt,
          GPS_external_clusters_csv = poppunk_gps_external_clusters_csv,
          docker = poppunk_docker_image, 
      }  
    }
  }
  if (merlin_tag == "Streptococcus pyogenes") {
    call emmtyper_task.emmtyper {
      input:
        assembly = assembly,
        samplename = samplename,
        wf = emmtyper_wf,
        cluster_distance = emmtyper_cluster_distance,
        percid = emmtyper_percid,
        culling_limit = emmtyper_culling_limit,
        mismatch = emmtyper_mismatch,
        align_diff = emmtyper_align_diff,
        gap = emmtyper_gap,
        min_perfect = emmtyper_min_perfect,
        min_good = emmtyper_min_good,
        max_size = emmtyper_max_size,
        docker = emmtyper_docker_image
    }
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
  if (merlin_tag == "Vibrio" || merlin_tag == "Vibrio cholerae") {
    if (!assembly_only && !ont_data) {
      call srst2_vibrio_task.srst2_vibrio {
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          srst2_min_cov = srst2_min_cov,
          srst2_max_divergence = srst2_max_divergence,
          srst2_min_depth = srst2_min_depth,
          srst2_min_edge_depth = srst2_min_edge_depth,
          srst2_gene_max_mismatch = srst2_gene_max_mismatch,
          docker = srst2_docker_image
      }
    }
    call abricate_vibrio_task.abricate_vibrio {
      input:
        assembly = assembly,
        samplename = samplename,
        minid = abricate_vibrio_minid,
        mincov = abricate_vibrio_mincov,
        docker = abricate_vibrio_docker_image
    }
  }
  
  # theiaeuk
  if (theiaeuk) {
    if (merlin_tag == "Candida auris") {
      call cauris_cladetyper.cauris_cladetyper as cladetyper {
        input: 
          assembly_fasta = assembly,
          samplename = samplename,
          kmer_size = cladetyper_kmer_size,
          ref_clade1 = cladetyper_ref_clade1,
          ref_clade1_annotated = cladetyper_ref_clade1_annotated,
          ref_clade2 = cladetyper_ref_clade2,
          ref_clade2_annotated = cladetyper_ref_clade2_annotated,
          ref_clade3 = cladetyper_ref_clade3,
          ref_clade3_annotated = cladetyper_ref_clade3_annotated,
          ref_clade4 = cladetyper_ref_clade4,
          ref_clade4_annotated = cladetyper_ref_clade4_annotated,
          ref_clade5 = cladetyper_ref_clade5,
          ref_clade5_annotated = cladetyper_ref_clade5_annotated,
          docker = cauris_cladetyper_docker_image
      }
      if (!assembly_only && !ont_data) {
        call snippy.snippy_variants as snippy_cauris { # no ONT support right now
          input:
            reference_genome_file = cladetyper.clade_spec_ref,
            read1 = select_first([read1]),
            read2 = read2,
            samplename = samplename,
            map_qual = snippy_map_qual,
            base_quality = snippy_base_quality,
            min_coverage = snippy_min_coverage,
            min_frac = snippy_min_frac,
            min_quality = snippy_min_quality,
            maxsoft = snippy_maxsoft,
            docker = snippy_variants_docker_image
        }
        call snippy_gene_query.snippy_gene_query as snippy_gene_query_cauris {
          input:
            samplename = samplename,
            snippy_variants_results = snippy_cauris.snippy_variants_results,
            reference = cladetyper.clade_spec_ref,
            query_gene = select_first([snippy_query_gene, "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase"]),
            docker = snippy_gene_query_docker_image
        }
      }
    }
    # Removing C.albicans subworkflow for now as current workflows not designed for diploid assembly
    # if (merlin_tag == "Candida albicans") {
    #   if (!assembly_only && !ont_data) {
    #     call snippy.snippy_variants as snippy_calbicans {
    #       input:
    #         reference_genome_file = snippy_reference_calbicans,
    #         read1 = select_first([read1]),
    #         read2 = read2,
            # map_qual = snippy_map_qual,
            # base_quality = snippy_base_quality,
            # min_coverage = snippy_min_coverage,
            # min_frac = snippy_min_frac,
            # min_quality = snippy_min_quality,
            # maxsoft = snippy_maxsoft,
    #         samplename = samplename
    #     }
    #     call snippy_gene_query.snippy_gene_query as snippy_gene_query_calbicans {
    #       input:
    #         samplename = samplename,
    #         snippy_variants_results = snippy_calbicans.snippy_variants_results,
    #         reference = snippy_reference_calbicans,
    #         query_gene = select_first([snippy_query_gene,"GCS1,ERG11,FUR1,RTA2"]), # GCS1 is another name for FKS1
    #     }
    #   }
    # }
    if (merlin_tag == "Aspergillus fumigatus") {
      if (!assembly_only && !ont_data) {
        call snippy.snippy_variants as snippy_afumigatus {
          input:
            reference_genome_file = snippy_reference_afumigatus,
            read1 = select_first([read1]),
            read2 = read2,
            samplename = samplename,
            map_qual = snippy_map_qual,
            base_quality = snippy_base_quality,
            min_coverage = snippy_min_coverage,
            min_frac = snippy_min_frac,
            min_quality = snippy_min_quality,
            maxsoft = snippy_maxsoft,
            docker = snippy_variants_docker_image
        }
        call snippy_gene_query.snippy_gene_query as snippy_gene_query_afumigatus {
          input:
            samplename = samplename,
            snippy_variants_results = snippy_afumigatus.snippy_variants_results,
            reference = snippy_reference_afumigatus,
            query_gene = select_first([snippy_query_gene, "Cyp51A,HapE,AFUA_4G08340"]), # AFUA_4G08340 is COX10 according to MARDy
            docker = snippy_gene_query_docker_image
        }
      }
    }
    if (merlin_tag == "Cryptococcus neoformans") {
      if (!assembly_only && !ont_data) {
        call snippy.snippy_variants as snippy_crypto {
          input:
            reference_genome_file = snippy_reference_cryptoneo,
            read1 = select_first([read1]),
            read2 = read2,
            samplename = samplename,
            map_qual = snippy_map_qual,
            base_quality = snippy_base_quality,
            min_coverage = snippy_min_coverage,
            min_frac = snippy_min_frac,
            min_quality = snippy_min_quality,
            maxsoft = snippy_maxsoft,
            docker = snippy_variants_docker_image
        }
        call snippy_gene_query.snippy_gene_query as snippy_gene_query_crypto {
          input:
            samplename = samplename,
            snippy_variants_results = snippy_crypto.snippy_variants_results,
            reference = snippy_reference_cryptoneo,
            query_gene = select_first([snippy_query_gene, "CNA00300"]), # CNA00300 is ERG11 for this reference genome
            docker = snippy_gene_query_docker_image
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
    String? kleborate_virulence_score = kleborate.kleborate_virulence_score
    String? kleborate_resistance_score = kleborate.kleborate_resistance_score
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
    # Acinetobacter baumannii typing
    File? abricate_abaum_results = abricate_abaum.abricate_results
    String? abricate_abaum_genes = abricate_abaum.abricate_genes
    String? abricate_abaum_database = abricate_abaum.abricate_database
    String? abricate_abaum_version = abricate_abaum.abricate_version
    String? abricate_abaum_docker = abricate_abaum.abricate_docker
    # Mycobacterium Typing
    File? tbprofiler_output_file = tbprofiler.tbprofiler_output_csv
    File? tbprofiler_output_bam = tbprofiler.tbprofiler_output_bam
    File? tbprofiler_output_bai = tbprofiler.tbprofiler_output_bai
    File? tbprofiler_output_vcf = tbprofiler.tbprofiler_output_vcf
    String? tbprofiler_version = tbprofiler.version
    String? tbprofiler_main_lineage = tbprofiler.tbprofiler_main_lineage
    String? tbprofiler_sub_lineage = tbprofiler.tbprofiler_sub_lineage
    String? tbprofiler_dr_type = tbprofiler.tbprofiler_dr_type
    String? tbprofiler_resistance_genes = tbprofiler.tbprofiler_resistance_genes
    Int? tbprofiler_median_coverage = tbprofiler.tbprofiler_median_coverage
    Float? tbprofiler_pct_reads_mapped = tbprofiler.tbprofiler_pct_reads_mapped
    String? tbp_parser_version = tbp_parser.tbp_parser_version
    String? tbp_parser_docker = tbp_parser.tbp_parser_docker
    File? tbp_parser_lims_report_csv = tbp_parser.tbp_parser_lims_report_csv
    File? tbp_parser_laboratorian_report_csv = tbp_parser.tbp_parser_laboratorian_report_csv
    File? tbp_parser_looker_report_csv = tbp_parser.tbp_parser_looker_report_csv
    File? tbp_parser_coverage_report = tbp_parser.tbp_parser_coverage_report
    Float? tbp_parser_genome_percent_coverage = tbp_parser.tbp_parser_genome_percent_coverage
    Float? tbp_parser_average_genome_depth = tbp_parser.tbp_parser_average_genome_depth
    File? clockwork_cleaned_read1 = clockwork_decon_reads.clockwork_cleaned_read1
    File? clockwork_cleaned_read2 = clockwork_decon_reads.clockwork_cleaned_read2
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
    String? emmtyper_emm_type = emmtyper.emmtyper_emm_type
    File? emmtyper_results_tsv = emmtyper.emmtyper_results_tsv
    String? emmtyper_version = emmtyper.emmtyper_version
    String? emmtyper_docker = emmtyper.emmtyper_docker
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
    String? srst2_vibrio_docker = srst2_vibrio.srst2_docker
    String? srst2_vibrio_database = srst2_vibrio.srst2_database
    String? srst2_vibrio_version = srst2_vibrio.srst2_version
    String? srst2_vibrio_ctxA = srst2_vibrio.srst2_vibrio_ctxA
    String? srst2_vibrio_ompW = srst2_vibrio.srst2_vibrio_ompW
    String? srst2_vibrio_toxR = srst2_vibrio.srst2_vibrio_toxR
    String? srst2_vibrio_serogroup = srst2_vibrio.srst2_vibrio_serogroup
    String? srst2_vibrio_biotype = srst2_vibrio.srst2_vibrio_biotype
    File? abricate_vibrio_detailed_tsv = abricate_vibrio.abricate_vibrio_results
    String? abricate_vibrio_database = abricate_vibrio.abricate_vibrio_database
    String? abricate_vibrio_docker = abricate_vibrio.abricate_vibrio_docker
    String? abricate_vibrio_version = abricate_vibrio.abricate_vibrio_version
    String? abricate_vibrio_ctxA = abricate_vibrio.abricate_vibrio_ctxA
    String? abricate_vibrio_ompW = abricate_vibrio.abricate_vibrio_ompW
    String? abricate_vibrio_toxR = abricate_vibrio.abricate_vibrio_toxR
    String? abricate_vibrio_biotype = abricate_vibrio.abricate_vibrio_biotype
    String? abricate_vibrio_serogroup = abricate_vibrio.abricate_vibrio_serogroup
    
    # theiaeuk
    # c auris 
    String? clade_type = cladetyper.gambit_cladetype
    String? cladetyper_analysis_date = cladetyper.date
    String? cladetyper_version = cladetyper.version
    String? cladetyper_docker_image = cladetyper.gambit_cladetyper_docker_image
    String? cladetype_annotated_ref = cladetyper.clade_spec_ref
    # snippy variants
    String snippy_variants_reference_genome = select_first([snippy_cauris.snippy_variants_reference_genome, snippy_afumigatus.snippy_variants_reference_genome, snippy_crypto.snippy_variants_reference_genome, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    String snippy_variants_version = select_first([snippy_cauris.snippy_variants_version, snippy_afumigatus.snippy_variants_version, snippy_crypto.snippy_variants_version, "No matching taxon detected"])
    String snippy_variants_query = select_first([snippy_gene_query_cauris.snippy_variants_query, snippy_gene_query_afumigatus.snippy_variants_query, snippy_gene_query_crypto.snippy_variants_query, "No matching taxon detected"])
    String snippy_variants_query_check = select_first([snippy_gene_query_cauris.snippy_variants_query_check, snippy_gene_query_afumigatus.snippy_variants_query_check, snippy_gene_query_crypto.snippy_variants_query_check, "No matching taxon detected"])
    String snippy_variants_hits = select_first([snippy_gene_query_cauris.snippy_variants_hits, snippy_gene_query_afumigatus.snippy_variants_hits, snippy_gene_query_crypto.snippy_variants_hits, "No matching taxon detected"])
    String snippy_variants_gene_query_results = select_first([snippy_gene_query_cauris.snippy_variants_gene_query_results, snippy_gene_query_afumigatus.snippy_variants_gene_query_results, snippy_gene_query_crypto.snippy_variants_gene_query_results, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    String snippy_variants_outdir_tarball = select_first([snippy_cauris.snippy_variants_outdir_tarball, snippy_afumigatus.snippy_variants_outdir_tarball, snippy_crypto.snippy_variants_outdir_tarball, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    String snippy_variants_results = select_first([snippy_cauris.snippy_variants_results, snippy_afumigatus.snippy_variants_results, snippy_crypto.snippy_variants_results, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    String snippy_variants_bam = select_first([snippy_cauris.snippy_variants_bam, snippy_afumigatus.snippy_variants_bam, snippy_crypto.snippy_variants_bam, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    String snippy_variants_bai = select_first([snippy_cauris.snippy_variants_bai, snippy_afumigatus.snippy_variants_bai, snippy_crypto.snippy_variants_bai, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    String snippy_variants_summary = select_first([snippy_cauris.snippy_variants_summary, snippy_afumigatus.snippy_variants_summary, snippy_crypto.snippy_variants_summary, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    String snippy_variants_num_reads_aligned = select_first([snippy_cauris.snippy_variants_num_reads_aligned, snippy_afumigatus.snippy_variants_num_reads_aligned, snippy_crypto.snippy_variants_num_reads_aligned, "No matching taxon detected"])
    String snippy_variants_coverage_tsv = select_first([snippy_cauris.snippy_variants_coverage_tsv, snippy_afumigatus.snippy_variants_coverage_tsv, snippy_crypto.snippy_variants_coverage_tsv, "gs://theiagen-public-files/terra/theiaeuk_files/no_match_detected.txt"])
    String snippy_variants_num_variants = select_first([snippy_cauris.snippy_variants_num_variants, snippy_afumigatus.snippy_variants_num_variants, snippy_crypto.snippy_variants_num_reads_aligned, "No matching taxon detected"])
    String snippy_variants_percent_ref_coverage = select_first([snippy_cauris.snippy_variants_percent_ref_coverage, snippy_afumigatus.snippy_variants_percent_ref_coverage, snippy_crypto.snippy_variants_percent_ref_coverage, "No matching taxon detected"])
  }
}