version 1.0

task tbprofiler_output_parsing {
  input {
    File json
    String output_seq_method_type
    String samplename
    Int min_depth = 10
  }
  command <<<
    python3 <<CODE
    import csv
    import json
    import pandas as pd

    ## Laboratorian report generation
    def parse_json_lab_report(json_file):

      gene_name = []
      locus_tag = []
      variant_substitutions = []
      confidence = []
      depth = []
      frequency = []

      with open(json_file) as results_json_fh:
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

        with open("tbprofiler_laboratorian_report.csv", "wt") as report_fh:
          report_fh.write("tbprofiler_gene_name,tbprofiler_locus_tag,tbprofiler_variant_substitutions,confidence,depth,frequency,read_support,warning\n")
          for i in range(0, len(gene_name)):
            if not depth[i]:  # for cases when depth is null, it gets converted to 0
              depth[i] = 0
            warning = "Low depth coverage" if  depth[i] < int('~{min_depth}') else "" # warning when coverage is lower than the defined 'min_depth' times
            report_fh.write(gene_name[i] + ',' + locus_tag[i] + ',' + variant_substitutions[i] + ',' + confidence[i] + ',' + str(depth[i]) + ',' + str(frequency[i]) + ',' + str(int(depth[i]*frequency[i])) + ',' + warning + '\n')

    parse_json_lab_report("~{json}")

    """
    ## LIMS ingestion report generation
    gene_dict={'Sample': '~{samplename}',
               'gyrB': '', 'gyrA': '', 'fgd1': '', 'mshA': '', 'ccsA': '', 'rpoB': '', 'rpoC': '', 'mmpL5': '', 'mmpS5': '', 'mmpR5': '', 'rpsL': '', 'rplC': '', 'fbiC': '', 'Rv1258c': '', 'embR': '',
               'atpE': '', 'rrs': '', 'rrl': '', 'fabG1': '', 'inhA': '', 'rpsA': '', 'tlyA': '', 'ndh': '', 'katG': '', 'PPE35': '', 'Rv1979c': '', 'pncA': '', 'kasA': '', 'eis': '', 'ahpC': '', 'folC': '',
               'pepQ': '', 'ribD': '', 'Rv2752c': '', 'thyX': '', 'thyA': '', 'ald': '', 'fbiD': '', 'Rv3083': '', 'fprA': '', 'whiB7': '', 'Rv3236c': '', 'fbiA': '', 'fbiB': '', 'alr': '', 'rpoA': '', 
               'ddn': '', 'clpC1': '', 'panD': '', 'embC': '', 'embA': '', 'embB': '', 'aftB': '', 'ubiA': '', 'ethA': '', 'ethR': '', 'whiB6': '', 'gid': '',
               'rifampicin': '', 'isoniazid': '', 'pyrazinamide': '', 'ethambutol': '', 'streptomycin': '', 'fluoroquinolones': '', 'moxifloxacin': '', 'ofloxacin': '', 
               'levofloxacin': '', 'ciprofloxacin': '', 'aminoglycosides': '', 'amikacin': '', 'kanamycin': '', 'capreomycin': '', 'ethionamide': '', 'para-aminosalicylic_acid': '', 
               'cycloserine': '', 'linezolid': '', 'bedaquiline': '', 'clofazimine': '', 'delamanid': ''}

    with open("~{json}") as results_json_fh:
      results_json = json.load(results_json_fh)
      for dr_variant in results_json["dr_variants"]:  # reported mutation by tb-profiler, all confering resistance
        name = dr_variant["gene"]
        substitutions =str(dr_variant["type"] + ":" + dr_variant["nucleotide_change"] + "(" + dr_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
        # check for cases where depth is zero because it's a deletion
        depth = other_variant["depth"]
        if 'del' not in substitutions and depth < int('~{min_depth}'):
          gene_dict[name] = gene_dict[name] + ';' + 'low coverage'
        else: 
          gene_dict[name] = gene_dict[name] + ';' + substitutions
        for annotation in dr_variant['annotation']:
          if 'drug' in annotation.keys():
            gene_dict[annotation['drug']] = annotation['who_confidence']

      for other_variant in results_json["other_variants"]:  # mutations not reported by tb-profiler
        name = other_variant["gene"]
        substitutions =str(other_variant["type"] + ":" + other_variant["nucleotide_change"] + "(" + other_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
        # check for cases where depth is zero because it's a deletion
        depth = other_variant["depth"]
        if 'del' not in substitutions and depth < int('~{min_depth}'):
          gene_dict[name] = gene_dict[name] + ';' + 'low coverage'
        else: 
          gene_dict[name] = gene_dict[name] + ';' + substitutions
        if 'annotation' in other_variant.keys():
          for annotation in other_variant['annotation']:
            if 'drug' in annotation.keys():
              gene_dict[annotation['drug']] = annotation['who_confidence']

    with open("tmp_tbprofiler_additional_outputs.csv", "w") as additional_outputs_fh:
      for key in gene_dict.keys():
        additional_outputs_fh.write(key)
        additional_outputs_fh.write(',')
      additional_outputs_fh.write('\n')
      for value in gene_dict.values():
        additional_outputs_fh.write(value)
        additional_outputs_fh.write(',')

    # file clean up with pandas
    df_additional_outputs = pd.read_csv("tmp_tbprofiler_additional_outputs.csv")
    df_additional_outputs = df_additional_outputs.loc[:, ~df_additional_outputs.columns.str.contains('^Unnamed')]
    df_additional_outputs = df_additional_outputs.fillna('')
    for column in df_additional_outputs:
      df_additional_outputs[column] = df_additional_outputs[column].str.lstrip(';')    
    df_additional_outputs.to_csv("tbprofiler_additional_outputs.csv", index=False)

    # Looker file - both concatenated side to side
    # file to be ingested into CDPH LIMS system
    with open("temp_looker.csv", "wt") as additional_outputs_csv:
      additional_outputs_csv.write("tbprofiler_gene_name,tbprofiler_locus_tag,tbprofiler_variant_substitutions,confidence,tbprofiler_output_seq_method_type\n")
      additional_outputs_csv.write(";".join(gene_name) + "," + ";".join(locus_tag) + "," + ";".join(variant_substitutions) + ',' + ";".join(confidence) + ',' + "~{output_seq_method_type}")

    df_looker = pd.read_csv("temp_looker.csv")
    result = pd.concat([df_looker, df_additional_outputs], axis=1)
    result.to_csv("tbprofiler_looker.csv", index=False)
    """

    def parse_json_mutations(json_file):
      mutations_dict = {}

      with open(json_file) as js_fh:
        results_json = json.load(js_fh)
        for dr_variant in results_json["dr_variants"]:  # reported mutation by tb-profiler, all confering resistance
          name = dr_variant["gene"]
          substitution =str(dr_variant["type"] + ":" + dr_variant["nucleotide_change"] + "(" + dr_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
          if name not in mutations_dict.keys():
            mutations_dict[name] = substitution
          else:
            mutations_dict[name] = mutations_dict[name] + ';' + substitution
        for other_variant in results_json["other_variants"]:  # mutations not reported by tb-profiler
          if other_variant["type"] != "synonymous_variant":  # report all non-synonymous mutations
              name = other_variant["gene"]
              substitution =str(other_variant["type"] + ":" + other_variant["nucleotide_change"] + "(" + other_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
              if name not in mutations_dict.keys():
                mutations_dict[name] = substitution
              else:
                mutations_dict[name] = mutations_dict[name] + ';' + substitution

      return mutations_dict


    def rank_annotation(annotation):
      if annotation == "Assoc w R":
        return 1
      elif annotation == "Assoc w R - interim":
        return 2
      elif annotation == "Uncertain significance":
        return 3
      else:
        return 4

    def translate(annotation, drug):
      if annotation == "Not assoc w R":
        return "No resistance to {} detected".format(drug)
      elif annotation == "Uncertain significance":
        return "The detected genetic determinant(s) have uncertain significance, resistance to {} cannot be rulled out".format(drug)
      elif annotation == "Assoc w R - interim":
        return "The detected genetic determinant(s) have uncertain significance, resistance to {} cannot be rulled out".format(drug)
      elif annotation == "Assoc w R":
        return "Genetic determinant(s) associated with resistance to {} detected".format(drug)
      else:
        return "No resistance to {} detected".format(drug)

    def parse_json_resistance(json_file):
      resistance_dict = {}

      with open(json_file) as js_fh:
        results_json = json.load(js_fh)
        for dr_variant in results_json["dr_variants"]:  # reported mutation by tb-profiler, all confering resistance
          for antimicrobial in dr_variant["gene_associated_drugs"]:
            if antimicrobial not in resistance_dict.keys():
              resistance_dict[antimicrobial] = "Assoc w R"
            
        for other_variant in results_json["other_variants"]:  # mutations not reported by tb-profiler
          if other_variant["type"] != "synonymous_variant":  # report all non-synonymous mutations
              if "annotation" in other_variant.keys():
                for annotation in other_variant["annotation"]:
                  drug = annotation["drug"]
                  resistance = annotation["who_confidence"]
                if drug not in resistance_dict.keys():
                  resistance_dict[antimicrobial] = resistance
                else:
                  if rank_annotation(resistance_dict[antimicrobial]) < rank_annotation(resistance):
                    resistance_dict[antimicrobial] = resistance
      return resistance_dict
    

    def get_lineage(json_file):
      with open(json_file) as js_fh:
        results_json = json.load(js_fh)
        if results_json["main_lin"] != "":
          return "DNA of M. tuberculosis complex detected (not M. bovis)"
        else:
          return "DNA of M. bovis detected"

    # lookup dictionary - antimicrobial code to name
    antimicrobial_dict = {"M_DST_B01_INH": "isoniazid", "M_DST_C01_ETO": "ethionamide",
                          "M_DST_D01_RIF": "rifampicin", "M_DST_E01_PZA": "pyrazinamide",
                          "M_DST_F01_EMB": "ethambutol","M_DST_G01_STM": "sulfamethazine",
                          "M_DST_H01_AMK": "amikacin", "M_DST_I01_KAN": "kanamycin",
                          "M_DST_J01_CAP": "capreomycin", "M_DST_K01_MFX": "moxifloxacin",
                          "M_DST_L01_LFX": "levofloxacin", "M_DST_M01_BDQ": "bedaquiline",
                          "M_DST_N01_CFZ": "clofazimine", "M_DST_o01_LZD": "linezolid" 
                         }

    # lookup dictionary - gene name to gene column name
    gene_dict = {"M_DST_B01_INH": {"katG": "M_DST_B02_katG", "fabG1": "M_DST_B03_fabG1", 
                                   "inhA": "M_DST_B04_inhA","ahpC": "M_DST_B05_ahpC"},
                 "M_DST_C01_ETO": {"ethA": "M_DST_C02_ethA", "fabG1": "M_DST_C03_fabG1",
                                   "inhA": "M_DST_C04_inhA"},
                 "M_DST_D01_RIF": {"rpoB": "M_DST_D02_rpoB"},
                 "M_DST_E01_PZA": {"pncA": "M_DST_E02_pncA", "panD": "M_DST_E03_panD",
                                   "clpC1": "M_DST_E04_clpC1"},
                 "M_DST_F01_EMB": {"embA": "M_DST_F02_embA", "embB": "M_DST_F03_embB",
                                   "embC": "M_DST_F04_embC"},
                 "M_DST_G01_STM": {"rrs": "M_DST_G02_rrs", "rpsL": "M_DST_G03_rpsL",
                                   "gid": "M_DST_G04_gid", "whiB7": "M_DST_G05_whiB7",
                                   "Rv1258c": "M_DST_G06_Rv1258c"},
                 "M_DST_H01_AMK": {"rrs": "M_DST_H02_rrs", "eis": "M_DST_H03_eis"},
                 "M_DST_I01_KAN": {"rrs": "M_DST_I02_rrs", "eis": "M_DST_I03_eis"},
                 "M_DST_J01_CAP": {"rrs": "M_DST_J02_rrs", "tlyA": "M_DST_J03_tlyA"},
                 "M_DST_K01_MFX": {"gyrA": "M_DST_K02_gyrA", "gyrB": "M_DST_K03_gyrB"},
                 "M_DST_L01_LFX": {"gyrA": "M_DST_L02_gyrA", "gyrB": "M_DST_L03_gyrB"},
                 "M_DST_M01_BDQ": {"Rv0678": "M_DST_M02_Rv0678", "atpE": "M_DST_M03_atpE",
                                   "pepQ": "M_DST_M04_pepQ", "mmpL5": "M_DST_M05_mmpL5",
                                   "mmpS5": "M_DST_M06_mmpS5"},
                 "M_DST_N01_CFZ": {"Rv0678":"M_DST_N02_Rv0678", "pepQ": "M_DST_N03_pepQ",
                                   "mmpL5":"M_DST_N04_mmpL5", "mmpS5": "M_DST_N05_mmpS5"},
                 "M_DST_o01_LZD": {"rrl": "M_DST_o02_rrl", "rplC": "M_DST_o03_rplC"}
                }
    
    lineage = get_lineage("~{json}")
    mutations = parse_json_mutations("~{json}")
    resistance = parse_json_resistance("~{json}")
    df_lims = pd.DataFrame({"MDL sample accession numbers":"~{samplename}", "M_DST_A01_ID": lineage},index=[0])

    for antimicrobial, genes in gene_dict.items():
      if antimicrobial_dict[antimicrobial] in resistance.keys():
        df_lims[antimicrobial] = translate(resistance[antimicrobial_dict[antimicrobial]], antimicrobial_dict[antimicrobial])
      else:
        df_lims[antimicrobial] = "No resistance to {} detected".format(antimicrobial_dict[antimicrobial]) # what to write when no resistance is reported?
      for gene_name, gene_id in genes.items():
        if gene_name in mutations.keys():
          df_lims[gene_id] = mutations[gene_name]
        else:
          df_lims[gene_id] = "No mutations detected"
    
    df_lims.to_csv("tbprofiler_lims_report.csv", index=False)
    CODE
  >>>
  output {
    #File tbprofiler_looker_csv = "tbprofiler_looker.csv"
    #File tbprofiler_additional_outputs_csv = "tbprofiler_additional_outputs.csv"
    File tbprofiler_laboratorian_report_csv = "tbprofiler_laboratorian_report.csv"
    File tbprofiler_lims_report_csv = "tbprofiler_lims_report.csv"
  }
  runtime {
    docker: "quay.io/theiagen/utility:1.2"
    memory: "4 GB"
    cpu: 1
    disks: "local-disk " + 10 + " SSD"
    disk: 10 + " GB"
    maxRetries: 0 
  }
}