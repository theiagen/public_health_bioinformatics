version 1.0

task kaptive {
  # Inputs
  input {
    File assembly
    String samplename
    String kaptive_docker_image = "quay.io/staphb/kaptive:2.0.3"
    Int cpu = 4
    # Parameters
    Int start_end_margin = 10 # determines flexibility in identifying the start and end of a locus - if this value is 10, a locus match that is missing the first 8 base pairs will still count as capturing the start of the locus (default: 10) 
    Float min_identity = 90.0 # minimum required percent coverage for the gene BLAST search via tBLASTn (default: 90.0)
    Float min_coverage = 80.0 # minimum required percent identity for the gene BLAST search via tBLASTn (default: 80.0)
    Float low_gene_id = 95.0 # percent identity threshold for what counts as a low identity match in the gene BLAST search (default: 95.0)
    #Int min_assembly_piece = 100 # smallest piece of the assembly (measured in bases) that will be included in the output FASTA files (default: 100)
    #Int gap_fill_size = 100 # size of assembly gaps to be filled in when producing the output FASTA files
  }

  command <<<
    #find absolute path of kaptive directory
    KAPTIVE_DIR=$(dirname "$(which kaptive.py)")
    # capture date and version
    # Print and save date
    date | tee DATE
    # Print and save version
    kaptive.py --version | tee VERSION 
    # Run Kaptive on the input assembly with the --all flag and output with samplename prefix
    kaptive.py \
    -t ~{cpu} \
    ~{'--start_end_margin ' + start_end_margin} \
    ~{'--min_gene_id ' + min_identity} \
    ~{'--min_gene_cov ' + min_coverage} \
    ~{'--low_gene_id ' + low_gene_id} \
    --no_seq_out \
    --no_json \
    --out ~{samplename}_kaptive_out_k \
    --assembly ~{assembly} \
    --k_refs ${KAPTIVE_DIR}/reference_database/Acinetobacter_baumannii_k_locus_primary_reference.gbk
    # parse outputs
    python3 <<CODE
    import csv
    with open("./~{samplename}_kaptive_out_k_table.txt",'r') as tsv_file:
      tsv_reader=csv.reader(tsv_file, delimiter="\t")
      tsv_data=list(tsv_reader)
      tsv_dict=dict(zip(tsv_data[0], tsv_data[1]))
      with open ("BEST_MATCH_LOCUS_K", 'wt') as Best_Match_Locus_K:
        kaptive_locus_k=tsv_dict['Best match locus']
        Best_Match_Locus_K.write(kaptive_locus_k)
      with open ("BEST_MATCH_TYPE_K", 'wt') as Best_Match_Type_K:
        kaptive_type_k=tsv_dict['Best match type']
        Best_Match_Type_K.write(kaptive_type_k)
      with open ("MATCH_CONFIDENCE_K", 'wt') as Match_Confidence_K:
        kaptive_confidence_k=tsv_dict['Match confidence']
        Match_Confidence_K.write(kaptive_confidence_k)
      with open ("NUM_EXPECTED_INSIDE_K", 'wt') as Num_Expected_Inside_K:
        expected_count_inside_k=tsv_dict['Expected genes in locus']
        Num_Expected_Inside_K.write(expected_count_inside_k)
      with open ("EXPECTED_GENES_IN_LOCUS_K", 'wt') as Expected_Inside_K:
        expected_in_k=tsv_dict['Expected genes in locus, details']
        Expected_Inside_K.write(expected_in_k)
      with open ("NUM_EXPECTED_OUTSIDE_K", 'wt') as Num_Expected_Outside_K:
        expected_count_outside_k=tsv_dict['Expected genes outside locus']
        Num_Expected_Outside_K.write(expected_count_outside_k)
      with open ("EXPECTED_GENES_OUT_LOCUS_K", 'wt') as Expected_Outside_K:
        expected_out_k=tsv_dict['Expected genes outside locus, details']
        Expected_Outside_K.write(expected_out_k)
      with open ("NUM_OTHER_INSIDE_K", 'wt') as Num_Other_Inside_K:
        other_count_inside_k=tsv_dict['Other genes in locus']
        Num_Other_Inside_K.write(other_count_inside_k)
      with open ("OTHER_GENES_IN_LOCUS_K", 'wt') as Other_Inside_K:
        other_in_k=tsv_dict['Other genes in locus, details']
        Other_Inside_K.write(other_in_k)
      with open ("NUM_OTHER_OUTSIDE_K", 'wt') as Num_Other_Outside_K:
        other_count_outside_k=tsv_dict['Other genes outside locus']
        Num_Other_Outside_K.write(other_count_outside_k)
      with open ("OTHER_GENES_OUT_LOCUS_K", 'wt') as Other_Outside_K:
        other_out_k=tsv_dict['Expected genes outside locus, details']
        Other_Outside_K.write(other_out_k)
    CODE
    kaptive.py \
    -t ~{cpu} \
    ~{'--start_end_margin ' + start_end_margin} \
    ~{'--min_gene_id ' + min_identity} \
    ~{'--min_gene_cov ' + min_coverage} \
    ~{'--low_gene_id ' + low_gene_id} \
    --no_seq_out \
    --no_json \
    --out ~{samplename}_kaptive_out_oc \
    --assembly ~{assembly} \
    --k_refs ${KAPTIVE_DIR}/reference_database/Acinetobacter_baumannii_OC_locus_primary_reference.gbk
    python3 <<CODE
    import csv
    with open("./~{samplename}_kaptive_out_oc_table.txt",'r') as tsv_file:
      tsv_reader=csv.reader(tsv_file, delimiter="\t")
      tsv_data=list(tsv_reader)
      tsv_dict=dict(zip(tsv_data[0], tsv_data[1]))
      with open ("BEST_MATCH_LOCUS_OC", 'wt') as Best_Match_Locus_OC:
        kaptive_locus_oc=tsv_dict['Best match locus']
        Best_Match_Locus_OC.write(kaptive_locus_oc)
      with open ("BEST_MATCH_TYPE_OC", 'wt') as Best_Match_Type_OC:
        kaptive_type_oc=tsv_dict['Best match type']
        Best_Match_Type_OC.write(kaptive_type_oc)
      with open ("MATCH_CONFIDENCE_OC", 'wt') as Match_Confidence_OC:
        kaptive_confidence_oc=tsv_dict['Match confidence']
        Match_Confidence_OC.write(kaptive_confidence_oc)
      with open ("NUM_EXPECTED_INSIDE_OC", 'wt') as Num_Expected_Inside_OC:
        expected_count_inside_oc=tsv_dict['Expected genes in locus']
        Num_Expected_Inside_OC.write(expected_count_inside_oc)
      with open ("EXPECTED_GENES_IN_LOCUS_OC", 'wt') as Expected_Inside_OC:
        expected_in_oc=tsv_dict['Expected genes in locus, details']
        Expected_Inside_OC.write(expected_in_oc)
      with open ("NUM_EXPECTED_OUTSIDE_OC", 'wt') as Num_Expected_Outside_OC:
        expected_count_outside_oc=tsv_dict['Expected genes outside locus']
        Num_Expected_Outside_OC.write(expected_count_outside_oc)
      with open ("EXPECTED_GENES_OUT_LOCUS_OC", 'wt') as Expected_Outside_OC:
        expected_out_oc=tsv_dict['Expected genes outside locus, details']
        Expected_Outside_OC.write(expected_out_oc)
      with open ("NUM_OTHER_INSIDE_OC", 'wt') as Num_Other_Inside_OC:
        other_count_inside_oc=tsv_dict['Other genes in locus']
        Num_Other_Inside_OC.write(other_count_inside_oc)
      with open ("OTHER_GENES_IN_LOCUS_OC", 'wt') as Other_Inside_OC:
        other_in_oc=tsv_dict['Other genes in locus, details']
        Other_Inside_OC.write(other_in_oc)
      with open ("NUM_OTHER_OUTSIDE_OC", 'wt') as Num_Other_Outside_OC:
        other_count_outside_oc=tsv_dict['Other genes outside locus']
        Num_Other_Outside_OC.write(other_count_outside_oc)
      with open ("OTHER_GENES_OUT_LOCUS_OC", 'wt') as Other_Outside_OC:
        other_out_oc=tsv_dict['Expected genes outside locus, details']
        Other_Outside_OC.write(other_out_oc)
    CODE
    mv -v ~{samplename}_kaptive_out_k_table.txt ~{samplename}_kaptive_out_k_table.tsv
    mv -v ~{samplename}_kaptive_out_oc_table.txt ~{samplename}_kaptive_out_oc_table.tsv
  >>>
  output {
    File kaptive_output_file_k = "~{samplename}_kaptive_out_k_table.tsv"
    File kaptive_output_file_oc = "~{samplename}_kaptive_out_oc_table.tsv"
    String kaptive_version = read_string("VERSION")
    String kaptive_k_match = read_string("BEST_MATCH_LOCUS_K")
    String kaptive_k_type = read_string("BEST_MATCH_TYPE_K")
    String kaptive_k_confidence = read_string("MATCH_CONFIDENCE_K")
    String kaptive_k_expected_inside_count = read_string("NUM_EXPECTED_INSIDE_K")
    String kaptive_k_expected_inside_genes = read_string("EXPECTED_GENES_IN_LOCUS_K")
    String kaptive_k_expected_outside_count = read_string("NUM_EXPECTED_OUTSIDE_K")
    String kaptive_k_expected_outside_genes = read_string("EXPECTED_GENES_OUT_LOCUS_K")
    String kaptive_k_other_inside_count = read_string("NUM_OTHER_INSIDE_K")
    String kaptive_k_other_inside_genes = read_string("OTHER_GENES_IN_LOCUS_K")
    String kaptive_k_other_outside_count = read_string("NUM_OTHER_OUTSIDE_K")
    String kaptive_k_other_outside_genes = read_string("OTHER_GENES_OUT_LOCUS_K")
    String kaptive_oc_match = read_string("BEST_MATCH_LOCUS_OC")
    String kaptive_oc_type = read_string("BEST_MATCH_TYPE_OC")
    String kaptive_oc_confidence = read_string("MATCH_CONFIDENCE_OC")
    String kaptive_oc_expected_inside_count = read_string("NUM_EXPECTED_INSIDE_OC")
    String kaptive_oc_expected_inside_genes = read_string("EXPECTED_GENES_IN_LOCUS_OC")
    String kaptive_oc_expected_outside_count = read_string("NUM_EXPECTED_OUTSIDE_OC")
    String kaptive_oc_expected_outside_genes = read_string("EXPECTED_GENES_OUT_LOCUS_OC")
    String kaptive_oc_other_inside_count = read_string("NUM_OTHER_INSIDE_OC")
    String kaptive_oc_other_inside_genes = read_string("OTHER_GENES_IN_LOCUS_OC")
    String kaptive_oc_other_outside_count = read_string("NUM_OTHER_OUTSIDE_OC")
    String kaptive_oc_other_outside_genes = read_string("OTHER_GENES_OUT_LOCUS_OC")
  }
  runtime {
    docker: "~{kaptive_docker_image}"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 SSD"
  }
}