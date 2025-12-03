version 1.0

task tbprofiler {
  input {
    File read1
    File? read2
    String samplename
    Boolean ont_data = false
    
    String mapper = "bwa"
    String variant_caller = "gatk"
    String? variant_calling_params
    
    String? additional_parameters # for tbprofiler
    
    Int min_depth = 10
    Float min_af = 0.1
    
    File? tbprofiler_custom_db
    String? tbdb_branch
    
    Int cpu = 8     
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/tbprofiler:6.6.3"
    Int memory = 16
  }
  command <<<
    # Print and save version
    tb-profiler version > VERSION && sed -i -e 's/TBProfiler version //' VERSION && sed -n -i '$p' VERSION
    
    # check if file is non existant or non empty
    if [ -z "~{read2}" ] || [ ! -s "~{read2}" ]; then
      INPUT_READS="-1 ~{read1}"
    else
      INPUT_READS="-1 ~{read1} -2 ~{read2}"
    fi

    # check if new database file is provided and not empty - if so, use that database preferentially
    if [ -s "~{tbprofiler_custom_db}" ]; then
      echo "Found new database file ~{tbprofiler_custom_db}"
      prefix=$(basename "~{tbprofiler_custom_db}" | sed 's/\.tar\.gz$//')
      tar xfv ~{tbprofiler_custom_db}
      
      tb-profiler load_library ./"$prefix"/"$prefix"
      TBDB="--db $prefix"

    elif [ -n "~{tbdb_branch}" ]; then
      echo "Using tbdb branch ~{tbdb_branch}"
      tb-profiler update_tbdb --branch ~{tbdb_branch}
      TBDB="--db ~{tbdb_branch}"
    
    else
      echo "Using default tbdb database"
      TBDB=""
    fi

    # Run tb-profiler on the input reads with samplename prefix
    tb-profiler profile \
      ${INPUT_READS} \
      --prefix ~{samplename} \
      --mapper ~{mapper} \
      --caller ~{variant_caller} \
      --calling_params "~{variant_calling_params}" \
      --depth ~{min_depth} \
      --af ~{min_af} \
      --threads ~{cpu} \
      --csv --txt \
      ~{true="--platform nanopore" false="" ont_data} \
      ~{additional_parameters} \
      ${TBDB}

    # Collate results
    tb-profiler collate --prefix ~{samplename}

    # convert any bcf files to vcf
    for bcf_file in ./vcf/*.bcf; do
      if [[ -f "$bcf_file" ]]; then
          bcftools view "$bcf_file" -Ov -o "${bcf_file%.bcf}.vcf"
      fi
    done

    # decompress any gzipped files. must be recompressed with 'bgzip' specifically for bcftools merge
    gunzip ./vcf/*.gz

    # bgzip compress and index all vcf files (must do this one file at a time)
    for vcf_file in ./vcf/*.vcf; do
      if [[ -f "$vcf_file" ]]; then
        bgzip "$vcf_file"
        bcftools index "${vcf_file}.gz"
      fi
    done

    # merge all vcf files (if we can) into a single vcf file
    vcf_count=$(ls ./vcf/*.vcf.gz 2>/dev/null | wc -l)
    if [ "$vcf_count" -eq 1 ]; then
      mv ./vcf/*.vcf.gz ./vcf/~{samplename}.targets.csq.merged.vcf.gz
    else
      bcftools merge --force-samples $(ls ./vcf/*.vcf.gz) > ./vcf/~{samplename}.targets.csq.merged.vcf.gz
    fi

    # decompress the final merged vcf file for output
    gunzip ./vcf/~{samplename}.targets.csq.merged.vcf.gz

    python3 <<CODE
    import csv
    import json
    import os

    with open("./~{samplename}.txt",'r') as tsv_file:
      tsv_reader=csv.reader(tsv_file, delimiter="\t")
      tsv_data=list(tsv_reader)
      tsv_dict=dict(zip(tsv_data[0], tsv_data[1]))

      with open ("MAIN_LINEAGE", 'wt') as main_lineage:
        main_lineage.write(tsv_dict['main_lineage'])
      with open ("SUB_LINEAGE", 'wt') as sublineage:
        sublineage.write(tsv_dict['sub_lineage'])
        
      with open ("DR_TYPE", 'wt') as dr_type:
        dr_type.write(tsv_dict['drtype'])
      with open ("NUM_DR_VARIANTS", 'wt') as num_dr_variants:
        num_dr_variants.write(tsv_dict['num_dr_variants'])
      with open ("NUM_OTHER_VARIANTS", 'wt') as num_other_variants:
        num_other_variants.write(tsv_dict['num_other_variants'])

      with open ("RESISTANCE_GENES", 'wt') as resistance_genes:
        res_genes_list=['rifampicin', 'isoniazid', 'ethambutol', 'pyrazinamide', 'moxifloxacin', 'levofloxacin', 'bedaquiline', 'delamanid', 'pretomanid', 'linezolid', 'streptomycin', 'amikacin', 'kanamycin', 'capreomycin', 'clofazimine', 'ethionamide', 'para-aminosalicylic_acid', 'cycloserine']
        res_genes=[]
        for i in res_genes_list:
          if tsv_dict[i] != '-':
            res_genes.append(tsv_dict[i])
        res_genes_string=';'.join(res_genes)
        resistance_genes.write(res_genes_string)

      with open ("MEDIAN_DEPTH", 'wt') as median_depth:
        median_depth.write(tsv_dict['target_median_depth'])
      with open ("PCT_READS_MAPPED", 'wt') as pct_reads_mapped:
        pct_reads_mapped.write(tsv_dict['pct_reads_mapped'])
    CODE
  >>>
  output {
    File tbprofiler_output_csv = "./results/~{samplename}.results.csv"
    File tbprofiler_output_tsv = "./results/~{samplename}.results.txt"
    File tbprofiler_output_json = "./results/~{samplename}.results.json"
    File tbprofiler_output_bam = "./bam/~{samplename}.bam"
    File tbprofiler_output_bai = "./bam/~{samplename}.bam.bai"
    File? tbprofiler_output_vcf = "./vcf/~{samplename}.targets.csq.merged.vcf"
    String version = read_string("VERSION")
    String tbprofiler_main_lineage = read_string("MAIN_LINEAGE")
    String tbprofiler_sub_lineage = read_string("SUB_LINEAGE")
    String tbprofiler_dr_type = read_string("DR_TYPE")
    String tbprofiler_num_dr_variants = read_string("NUM_DR_VARIANTS")
    String tbprofiler_num_other_variants = read_string("NUM_OTHER_VARIANTS")
    String tbprofiler_resistance_genes = read_string("RESISTANCE_GENES")
    Float tbprofiler_median_depth = read_float("MEDIAN_DEPTH")
    Float tbprofiler_pct_reads_mapped = read_float("PCT_READS_MAPPED")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 1
  }
}