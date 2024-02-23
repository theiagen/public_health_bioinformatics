version 1.0

task kleborate {
  # Inputs
  input {
    File assembly
    String samplename
    String kleborate_docker_image = "us-docker.pkg.dev/general-theiagen/staphb/kleborate:2.2.0"
    Int disk_size = 100
    Int cpu = 8
    Int memory = 16
    
    # Parameters
    Boolean skip_resistance = false # Turn on resistance genes screening (default: no resistance gene screening)
    Boolean skip_kaptive = false # Equivalent to --kaptive_k --kaptive_
    Float min_identity = 90.0 # Minimum alignment percent identity for main results (default: 90.0)
    Float min_coverage = 80.0 #  Minimum alignment percent coverage for main results (default: 80.0)
    Float min_spurious_identity = 80.0 # Minimum alignment percent identity for spurious results (default: 80.0)
    Float min_spurious_coverage = 40.0 #  Minimum alignment percent coverage for spurious results (default: 40.0)
    String min_kaptive_confidence = "Good" # {None,Low,Good,High,Very_high,Perfect} Minimum Kaptive confidence to call K/O loci - confidence levels below this will be reported as unknown (default: Good)
  }
  command <<<
    # Print and save date
    date | tee DATE
    
    # Print and save version
    kleborate --version | tee VERSION 

    # Run Kleborate on the input assembly with the --all flag and output with samplename prefix
    kleborate \
    ~{true="" false="--resistance" skip_resistance} \
    ~{true="" false="--kaptive" skip_kaptive} \
    ~{'--min_identity ' + min_identity} \
    ~{'--min_coverage ' + min_coverage} \
    ~{'--min_spurious_identity ' + min_spurious_identity} \
    ~{'--min_spurious_coverage ' + min_spurious_coverage} \
    ~{'--min_kaptive_confidence ' + min_kaptive_confidence} \
    --outfile ~{samplename}_kleborate_out.tsv \
    --assemblies ~{assembly} \
    --all
    
    # parse outputs
    python3 <<CODE
    import csv
    with open("./~{samplename}_kleborate_out.tsv",'r') as tsv_file:
      tsv_reader=csv.reader(tsv_file, delimiter="\t")
      tsv_data=list(tsv_reader)
      tsv_dict=dict(zip(tsv_data[0], tsv_data[1]))
      with open ("SPECIES", 'wt') as Species:
        kleb_species=tsv_dict['species']
        Species.write(kleb_species)
      with open ("MLST_SEQUENCE_TYPE", 'wt') as MLST_Sequence_Type:
        mlst=tsv_dict['ST']
        MLST_Sequence_Type.write(mlst)
      with open ("VIRULENCE_SCORE", 'wt') as Virulence_Score:
        virulence_level=tsv_dict['virulence_score']
        Virulence_Score.write(virulence_level)
      with open ("RESISTANCE_SCORE", 'wt') as Resistance_Score:
        resistance_level=tsv_dict['resistance_score']
        Resistance_Score.write(resistance_level)
      with open ("NUM_RESISTANCE_GENES", 'wt') as Num_Resistance_Genes:
        resistance_genes_count=tsv_dict['num_resistance_genes']
        Num_Resistance_Genes.write(resistance_genes_count)
      with open ("BLA_RESISTANCE_GENES", 'wt') as BLA_Resistance_Genes:
        bla_res_genes_list=['Bla_acquired', 'Bla_inhR_acquired', 'Bla_ESBL_acquired', 'Bla_ESBL_inhR_acquired', 'Bla_Carb_acquired']
        bla_res_genes=[]
        for i in bla_res_genes_list:
          if tsv_dict[i] != '-':
            bla_res_genes.append(tsv_dict[i])
        bla_res_genes_string=';'.join(bla_res_genes)
        BLA_Resistance_Genes.write(bla_res_genes_string)
      with open ("ESBL_RESISTANCE_GENES", 'wt') as ESBL_Resistance_Genes:
        esbl_res_genes_list=['Bla_ESBL_acquired', 'Bla_ESBL_inhR_acquired']
        esbl_res_genes=[]
        for i in esbl_res_genes_list:
          if tsv_dict[i] != '-':
            bla_res_genes.append(tsv_dict[i])
        esbl_res_genes_string=';'.join(esbl_res_genes)
        ESBL_Resistance_Genes.write(esbl_res_genes_string)
      with open ("KEY_RESISTANCE_GENES", 'wt') as Key_Resistance_Genes:
        key_res_genes_list= ['Col_acquired', 'Fcyn_acquired', 'Flq_acquired', 'Rif_acquired', 'Bla_acquired', 'Bla_inhR_acquired', 'Bla_ESBL_acquired', 'Bla_ESBL_inhR_acquired', 'Bla_Carb_acquired']
        key_res_genes=[]
        for i in key_res_genes_list:
          if tsv_dict[i] != '-':
            key_res_genes.append(tsv_dict[i])
        key_res_genes_string=';'.join(key_res_genes)
        Key_Resistance_Genes.write(key_res_genes_string)
      with open ("GENOMIC_RESISTANCE_MUTATIONS", 'wt') as Resistance_Mutations:
        res_mutations_list= ['Bla_chr', 'SHV_mutations', 'Omp_mutations', 'Col_mutations', 'Flq_mutations']
        res_mutations=[]
        for i in res_mutations_list:
          if tsv_dict[i] != '-':
            res_mutations.append(tsv_dict[i])
        res_mutations_string=';'.join(res_mutations)
        Resistance_Mutations.write(res_mutations_string)
      with open ("K_TYPE", 'wt') as Ktype:
        ktype=tsv_dict['K_type']
        Ktype.write(ktype)
      with open ("K_LOCUS", 'wt') as Klocus:
        klocus=tsv_dict['K_locus']
        Klocus.write(klocus)
      with open ("O_TYPE", 'wt') as Otype:
        otype=tsv_dict['O_type']
        Otype.write(otype)
      with open ("O_LOCUS", 'wt') as Olocus:
        olocus=tsv_dict['O_locus']
        Olocus.write(olocus)
      with open ("K_LOCUS_CONFIDENCE", 'wt') as K_locus_confidence:
        k_locus_confidence=tsv_dict['K_locus_confidence']
        K_locus_confidence.write(k_locus_confidence)
      with open ("O_LOCUS_CONFIDENCE", 'wt') as O_locus_confidence:
        o_locus_confidence=tsv_dict['O_locus_confidence']
        O_locus_confidence.write(o_locus_confidence)
    CODE
  >>>
  output {
    File kleborate_output_file = "~{samplename}_kleborate_out.tsv"
    String kleborate_version = read_string("VERSION")
    String kleborate_docker = kleborate_docker_image
    String kleborate_mlst_sequence_type = read_string("MLST_SEQUENCE_TYPE")
    String kleborate_virulence_score = read_string("VIRULENCE_SCORE")
    String kleborate_resistance_score = read_string("RESISTANCE_SCORE")
    String kleborate_num_resistance_genes = read_string("NUM_RESISTANCE_GENES")
    String kleborate_bla_resistance_genes = read_string("BLA_RESISTANCE_GENES")
    String kleborate_esbl_resistance_genes = read_string("ESBL_RESISTANCE_GENES")
    String kleborate_key_resistance_genes = read_string("KEY_RESISTANCE_GENES")
    String kleborate_genomic_resistance_mutations = read_string("GENOMIC_RESISTANCE_MUTATIONS")
    String kleborate_klocus = read_string("K_LOCUS")
    String kleborate_ktype = read_string("K_TYPE")
    String kleborate_olocus = read_string("O_LOCUS")
    String kleborate_otype = read_string("O_TYPE")
    String kleborate_klocus_confidence = read_string("K_LOCUS_CONFIDENCE")
    String kleborate_olocus_confidence = read_string("O_LOCUS_CONFIDENCE")
  }
  runtime {
    docker: "~{kleborate_docker_image}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}