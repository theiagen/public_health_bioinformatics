version 1.0

task tbprofiler {
  # Inputs
  input {
    File read1
    File? read2
    String samplename
    Boolean tbprofiler_additional_outputs
    String output_seq_method_type
    String tbprofiler_docker_image = "staphb/tbprofiler:4.4.2"
    Int disk_size = 100
    String mapper = "bwa"
    String caller = "bcftools"
    Int min_depth = 10
    Float min_af = 0.1
    Float min_af_pred = 0.1
    Int cov_frac_threshold = 1
    Int cpu = 8 
    Boolean ont_data = false
  }
  command <<<
    # Print and save date
    date | tee DATE

    # Print and save version
    tb-profiler --version > VERSION && sed -i -e 's/^/TBProfiler version /' VERSION
    
    if [ -z "~{read2}" ] ; then
      INPUT_READS="-1 ~{read1}"
    else
      INPUT_READS="-1 ~{read1} -2 ~{read2}"
    fi
    
    if [ "~{ont_data}" = true ]; then
      mode="--platform nanopore"
      export ont_data="true"
    else
      export ont_data="false"
    fi

    # Run tb-profiler on the input reads with samplename prefix
    tb-profiler profile \
      ${mode} \
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

    # Collate results
    tb-profiler collate --prefix ~{samplename}

    # touch optional output files because wdl
    touch GENE_NAME LOCUS_TAG VARIANT_SUBSTITUTIONS OUTPUT_SEQ_METHOD_TYPE

    # transform boolean tbprofiler_additional_outputs into string for python comparison
    if ~{tbprofiler_additional_outputs}; then
      export tbprofiler_additional_outputs="true"
    else 
      export tbprofiler_additional_outputs="false"
    fi

    python3 <<CODE
    import csv
    import json
    import os

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

    if (os.environ["tbprofiler_additional_outputs"] == "true"):

      gene_name = []
      locus_tag = []
      variant_substitutions = []
      confidence = []
      depth = []
      frequency = []

      with open("./results/~{samplename}.results.json") as results_json_fh:
        results_json = json.load(results_json_fh)
        for dr_variant in results_json["dr_variants"]:  # reported mutation by tb-profiler, all confering resistance
          gene_name.append(dr_variant["gene"])
          locus_tag.append(dr_variant["locus_tag"])  
          variant_substitutions.append(dr_variant["type"] + ":" + dr_variant["nucleotide_change"] + "(" + dr_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
          depth.append(dr_variant["depth"])
          frequency.append(dr_variant["freq"])
          if "annotation" in dr_variant:
            try:  # sometimes annotation is an empty list
              if dr_variant["annotation"][0]["who_confidence"] == "":
                confidence.append("No WHO annotation")
              else:
                confidence.append(dr_variant["annotation"][0]["who_confidence"])
            except:
              confidence.append("No WHO annotation")
          else:
            confidence.append("No WHO annotation")

        for other_variant in results_json["other_variants"]:  # mutations not reported by tb-profiler
          if other_variant["type"] != "synonymous_variant":
            if other_variant["gene"] == "katG" or other_variant["gene"] == "pncA" or other_variant["gene"] == "rpoB" or other_variant["gene"] == "ethA" or other_variant["gene"] == "gid":  # hardcoded for genes of interest that are reported to always confer resistance when mutated
              gene_name.append(other_variant["gene"])
              locus_tag.append(other_variant["locus_tag"])  
              variant_substitutions.append(other_variant["type"] + ":" + other_variant["nucleotide_change"] + "(" + other_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
              depth.append(other_variant["depth"])
              frequency.append(other_variant["freq"])
              if "annotation" in other_variant:
                try:  # sometimes annotation is an empty list
                  if other_variant["annotation"][0]["who_confidence"] == "":
                    confidence.append("No WHO annotation")
                  else:
                    confidence.append(other_variant["annotation"][0]["who_confidence"])
                except:
                  confidence.append("No WHO annotation")
              else:
                confidence.append("No WHO annotation")
            else:
              if "annotation" in other_variant:  # check if who annotation field is present
                for annotation in other_variant["annotation"]:
                  if annotation["who_confidence"] != "Not assoc w R" or annotation["who_confidence"] != "":
                    gene_name.append(other_variant["gene"])
                    locus_tag.append(other_variant["locus_tag"])  
                    variant_substitutions.append(other_variant["type"] + ":" + other_variant["nucleotide_change"] + "(" + other_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
                    depth.append(other_variant["depth"])
                    frequency.append(other_variant["freq"])
                    confidence.append(annotation["who_confidence"])

        with open("GENE_NAME", "wt") as gene_name_fh:
          gene_name_fh.write(','.join(gene_name))

        with open("LOCUS_TAG", "wt") as locus_tag_fh:
          locus_tag_fh.write(','.join(locus_tag))

        with open("VARIANT_SUBSTITUTIONS", "wt") as variant_substitutions_fh:
          variant_substitutions_fh.write(','.join(variant_substitutions))

        with open("OUTPUT_SEQ_METHOD_TYPE", "wt") as output_seq_method_type_fh:
          output_seq_method_type_fh.write("~{output_seq_method_type}")

        # file to be ingested into CDPH LIMS system
        with open("tbprofiler_additional_outputs.csv", "wt") as additional_outputs_csv:
          additional_outputs_csv.write("tbprofiler_gene_name,tbprofiler_locus_tag,tbprofiler_variant_substitutions,confidence,tbprofiler_output_seq_method_type\n")
          additional_outputs_csv.write(";".join(gene_name) + "," + ";".join(locus_tag) + "," + ";".join(variant_substitutions) + ',' + ";".join(confidence) + ',' + "~{output_seq_method_type}")

        # laboratorian report
        with open("tbprofiler_laboratorian_report.csv", "wt") as report_fh:
          report_fh.write("tbprofiler_gene_name,tbprofiler_locus_tag,tbprofiler_variant_substitutions,confidence,depth,frequency,warning\n")

          if (os.environ(ont_data) == "false"):
            for i in range(0, len(gene_name)):
              if not depth[i]:  # for cases when depth is null, it gets converted to 0
                depth[i] = 0
              warning = "Low depth coverage" if  depth[i] < int('~{min_depth}') else "" # warning when coverage is lower than the defined 'min_depth' times
              report_fh.write(gene_name[i] + ',' + locus_tag[i] + ',' + variant_substitutions[i] + ',' + confidence[i] + ',' + str(depth[i]) + ',' + str(frequency[i]) + ',' + warning + '\n')
          
          if (os.environ(ont_data) == "true"):
            for i in range(0, len(gene_name)):
              warning = []
              if not depth[i]:  # for cases when depth is null
                depth[i] = "NA"  # Deletion instead of value
                warning.append("Deletion")
              else:
                if depth[i] < int('~{min_depth}'): # warning when coverage is lower than the defined 'min_depth' times
                  warning.append("Low depth coverage") 
              if not frequency[i]:
                warning.append("No frequency")
              else:
                if frequency[i] < 0.05: # warning when frequency is lower than 5%
                  warning.append("Low frequency") 
              report_fh.write(gene_name[i] + ',' + locus_tag[i] + ',' + variant_substitutions[i] + ',' + confidence[i] + ',' + str(depth[i]) + ',' + str(frequency[i]) + ',' + ';'.join(warning) + '\n')
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
    File? tbprofiler_additional_outputs_csv = "tbprofiler_additional_outputs.csv"
    File? tbprofiler_laboratorian_report_csv = "tbprofiler_laboratorian_report.csv"
    String tbprofiler_gene_name = read_string("GENE_NAME")
    String tbprofiler_locus_tag = read_string("LOCUS_TAG")
    String tbprofiler_variant_substitutions = read_string("VARIANT_SUBSTITUTIONS")
    String tbprofiler_output_seq_method_type = read_string("OUTPUT_SEQ_METHOD_TYPE")
  }
  runtime {
    docker: "~{tbprofiler_docker_image}"
    memory: "16 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3 
  }
}