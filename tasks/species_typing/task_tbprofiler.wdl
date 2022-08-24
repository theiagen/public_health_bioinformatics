version 1.0

task tbprofiler {
  # Inputs
  input {
    File read1
    File? read2
    String samplename
    String tbprofiler_docker_image = "quay.io/biocontainers/tb-profiler:3.0.8--pypyh5e36f6f_0"
    String? mapper = "bwa"
    String? caller = "bcftools"
    Int? min_depth = 10
    Float? min_af = 0.1
    Float? min_af_pred = 0.1
    Int? cov_frac_threshold = 1
  }
  command <<<
    # update TBDB
    # tb-profiler update_tbdb
    # Print and save date
    date | tee DATE
    # Print and save version
    tb-profiler --version > VERSION && sed -i -e 's/^/TBProfiler version /' VERSION
    
    if [ -z "~{read2}" ] ; then
      INPUT_READS="-1 ~{read1}"
    else
      INPUT_READS="-1 ~{read1} -2 ~{read2}"
    fi

    # Run Kleborate on the input assembly with the --all flag and output with samplename prefix
    tb-profiler profile \
      ${INPUT_READS} \
      --prefix ~{samplename} \
      --mapper ~{mapper} \
      --caller ~{caller} \
      --min_depth ~{min_depth} \
      --af ~{min_af} \
      --reporting_af \
      ~{min_af_pred} \
      --coverage_fraction_threshold ~{cov_frac_threshold} \
      --csv --txt

    #Collate results
    tb-profiler collate --prefix ~{samplename}

    python3 <<CODE
    import csv
    with open("./~{samplename}.txt",'r') as tsv_file:
      tsv_reader=csv.reader(tsv_file, delimiter="\t")
      tsv_data=list(tsv_reader)
      tsv_dict=dict(zip(tsv_data[0], tsv_data[1]))
      with open ("MAIN_LINEAGE", 'wt') as Main_Lineage:
        main_lin=tsv_dict['main_lineage']
        Main_Lineage.write(main_lin)
      with open ("SUB_LINEAGE", 'wt') as Sub_Lineage:
        sub_lin=tsv_dict['sub_lineage']
        Sub_Lineage.write(sub_lin)
      with open ("DR_TYPE", 'wt') as DR_Type:
        dr_type=tsv_dict['DR_type']
        DR_Type.write(dr_type)
      with open ("NUM_DR_VARIANTS", 'wt') as Num_DR_Variants:
        num_dr_vars=tsv_dict['num_dr_variants']
        Num_DR_Variants.write(num_dr_vars)
      with open ("NUM_OTHER_VARIANTS", 'wt') as Num_Other_Variants:
        num_other_vars=tsv_dict['num_other_variants']
        Num_Other_Variants.write(num_other_vars)
      with open ("RESISTANCE_GENES", 'wt') as Resistance_Genes:
        res_genes_list=['rifampicin', 'isoniazid', 'pyrazinamide', 'ethambutol', 'streptomycin', 'fluoroquinolones', 'moxifloxacin', 'ofloxacin', 'levofloxacin', 'ciprofloxacin', 'aminoglycosides', 'amikacin', 'kanamycin', 'capreomycin', 'ethionamide', 'para-aminosalicylic_acid', 'cycloserine', 'linezolid', 'bedaquiline', 'clofazimine', 'delamanid']
        res_genes=[]
        for i in res_genes_list:
          if tsv_dict[i] != '-':
            res_genes.append(tsv_dict[i])
        res_genes_string=';'.join(res_genes)
        Resistance_Genes.write(res_genes_string)
    CODE
  >>>
  output {
    File tbprofiler_output_csv = "./results/~{samplename}.results.csv"
    File tbprofiler_output_tsv = "./results/~{samplename}.results.txt"
    File tbprofiler_output_bam = "./bam/~{samplename}.bam"
    File tbprofiler_output_bai = "./bam/~{samplename}.bam.bai"
    String version = read_string("VERSION")
    String tbprofiler_main_lineage = read_string("MAIN_LINEAGE")
    String tbprofiler_sub_lineage = read_string("SUB_LINEAGE")
    String tbprofiler_dr_type = read_string("DR_TYPE")
    String tbprofiler_num_dr_variants = read_string("NUM_DR_VARIANTS")
    String tbprofiler_num_other_variants = read_string("NUM_OTHER_VARIANTS")
    String tbprofiler_resistance_genes = read_string("RESISTANCE_GENES")
  }
  runtime {
    docker: "~{tbprofiler_docker_image}"
    memory: "16 GB"
    cpu: 8
    disks: "local-disk 100 SSD"
  }
}

task tbprofiler_ont {
  # Inputs
  input {
    File reads
    String samplename
    String tbprofiler_docker_image = "quay.io/biocontainers/tb-profiler:3.0.8--pypyh5e36f6f_0"
    String? mapper = "bwa"
    String? caller = "bcftools"
    Int? min_depth = 10
    Float? min_af = 0.1
    Float? min_af_pred = 0.1
    Int? cov_frac_threshold = 1
  }
  command <<<
    # update TBDB
    # tb-profiler update_tbdb
    # Print and save date
    date | tee DATE
    # Print and save version
    tb-profiler --version > VERSION && sed -i -e 's/^/TBProfiler version /' VERSION
    # Run TBProfiler on the input sample
    tb-profiler profile --platform nanopore -1 ~{reads} --prefix ~{samplename} --mapper ~{mapper} --caller ~{caller} --min_depth ~{min_depth} --af ~{min_af} --reporting_af ~{min_af_pred} --coverage_fraction_threshold ~{cov_frac_threshold} --csv --txt

    #Collate results
    tb-profiler collate --prefix ~{samplename}

    python3 <<CODE
    import csv
    with open("./~{samplename}.txt",'r') as tsv_file:
      tsv_reader=csv.reader(tsv_file, delimiter="\t")
      tsv_data=list(tsv_reader)
      tsv_dict=dict(zip(tsv_data[0], tsv_data[1]))
      with open ("MAIN_LINEAGE", 'wt') as Main_Lineage:
        main_lin=tsv_dict['main_lineage']
        Main_Lineage.write(main_lin)
      with open ("SUB_LINEAGE", 'wt') as Sub_Lineage:
        sub_lin=tsv_dict['sub_lineage']
        Sub_Lineage.write(sub_lin)
      with open ("DR_TYPE", 'wt') as DR_Type:
        dr_type=tsv_dict['DR_type']
        DR_Type.write(dr_type)
      with open ("NUM_DR_VARIANTS", 'wt') as Num_DR_Variants:
        num_dr_vars=tsv_dict['num_dr_variants']
        Num_DR_Variants.write(num_dr_vars)
      with open ("NUM_OTHER_VARIANTS", 'wt') as Num_Other_Variants:
        num_other_vars=tsv_dict['num_other_variants']
        Num_Other_Variants.write(num_other_vars)
      with open ("RESISTANCE_GENES", 'wt') as Resistance_Genes:
        res_genes_list=['rifampicin', 'isoniazid', 'pyrazinamide', 'ethambutol', 'streptomycin', 'fluoroquinolones', 'moxifloxacin', 'ofloxacin', 'levofloxacin', 'ciprofloxacin', 'aminoglycosides', 'amikacin', 'kanamycin', 'capreomycin', 'ethionamide', 'para-aminosalicylic_acid', 'cycloserine', 'linezolid', 'bedaquiline', 'clofazimine', 'delamanid']
        res_genes=[]
        for i in res_genes_list:
          if tsv_dict[i] != '-':
            res_genes.append(tsv_dict[i])
        res_genes_string=';'.join(res_genes)
        Resistance_Genes.write(res_genes_string)
    CODE
  >>>
  output {
    File tbprofiler_output_csv = "./results/~{samplename}.results.csv"
    File tbprofiler_output_tsv = "./results/~{samplename}.results.txt"
    File tbprofiler_output_bam = "./bam/~{samplename}.bam"
    File tbprofiler_output_bai = "./bam/~{samplename}.bam.bai"
    String version = read_string("VERSION")
    String tbprofiler_main_lineage = read_string("MAIN_LINEAGE")
    String tbprofiler_sub_lineage = read_string("SUB_LINEAGE")
    String tbprofiler_dr_type = read_string("DR_TYPE")
    String tbprofiler_num_dr_variants = read_string("NUM_DR_VARIANTS")
    String tbprofiler_num_other_variants = read_string("NUM_OTHER_VARIANTS")
    String tbprofiler_resistance_genes = read_string("RESISTANCE_GENES")
  }
  runtime {
    docker: "~{tbprofiler_docker_image}"
    memory: "16 GB"
    cpu: 8
    disks: "local-disk 100 SSD"
    maxRetries: 3
  }
}
