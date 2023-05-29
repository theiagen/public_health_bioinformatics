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
    import re
    import pandas as pd

    ## Lookup Dictionaries ##

    # lookup dictionary - antimicrobial code to name
    antimicrobial_dict = {"M_DST_B01_INH": "isoniazid", "M_DST_C01_ETO": "ethionamide",
                          "M_DST_D01_RIF": "rifampicin", "M_DST_E01_PZA": "pyrazinamide",
                          "M_DST_F01_EMB": "ethambutol","M_DST_G01_STM": "sulfamethazine",
                          "M_DST_H01_AMK": "amikacin", "M_DST_I01_KAN": "kanamycin",
                          "M_DST_J01_CAP": "capreomycin", "M_DST_K01_MFX": "moxifloxacin",
                          "M_DST_L01_LFX": "levofloxacin", "M_DST_M01_BDQ": "bedaquiline",
                          "M_DST_N01_CFZ": "clofazimine", "M_DST_o01_LZD": "linezolid" 
                         }
    
    # Lookup list - antimicrobials
    antimicrobial_list = ["isoniazid", "ethionamide", "rifampicin", "pyrazinamide", "ethambutol",
                          "sulfamethazine", "amikacin", "kanamycin", "capreomycin", "moxifloxacin",
                          "levofloxacin", "bedaquiline", "clofazimine", "linezolid"
                         ]

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

    ## Auxiliary Functions ##

    def get_codon(protein_mut):
      """
      This function recieves a protein mutation (e.g. 'p.Met291Ile')
      and returns the codon (numerical part) where that mutation occurs
      in Int format.
      """
      pattern = r"\.\D+(\d+)"

      match = re.search(pattern, protein_mut)
      if match:
        position = match.group(1)
        return int(position)
      return 0


    def parse_json_lab_report(json_file):
      """
      This function recieves the tbprofiler output json file and
      writes the laboratorian report that includes the following information
      per mutation:
        - tbprofiler_gene_name
        - tbprofiler_locus_tag
        - tbprofiler_variant_substitutions
        - confidence
        - depth
        - frequency
        - read_support
        - rationale
        - warning
      """
      sample_id = []
      gene_name = []
      locus_tag = []
      variant_substitutions_type = []
      variant_substitutions_nt = []
      variant_substitutions_aa = []
      confidence = []
      resistance = []
      depth = []
      frequency = []
      rule = []

      with open(json_file) as results_json_fh:
        results_json = json.load(results_json_fh)
        for dr_variant in results_json["dr_variants"]:  # reported mutation by tb-profiler, all confering resistance
          sample_id.append("~{samplename}")
          gene_name.append(dr_variant["gene"])
          locus_tag.append(dr_variant["locus_tag"])  
          #variant_substitutions.append(dr_variant["type"] + ":" + dr_variant["nucleotide_change"] + "(" + dr_variant["protein_change"] + ")") 
          variant_substitutions_type.append(dr_variant["type"])
          variant_substitutions_nt.append(dr_variant["nucleotide_change"])
          variant_substitutions_aa.append(dr_variant["protein_change"] if dr_variant["protein_change"] != "" else "NA")
          depth.append(dr_variant["depth"])
          frequency.append(dr_variant["freq"])
          rule.append("WHO classification")
          resistance.append(dr_variant["gene_associated_drugs"][0])
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
            if other_variant["gene"] == "katG" or other_variant["gene"] == "pncA" or other_variant["gene"] == "ethA" or other_variant["gene"] == "gid":  # Expert rule: hardcoded for genes of interest that are reported to always confer resistance when mutated
              # report as uncertain significance based on expert rule
              sample_id.append("~{samplename}")
              gene_name.append(other_variant["gene"])
              locus_tag.append(other_variant["locus_tag"])  
              #variant_substitutions.append(other_variant["type"] + ":" + other_variant["nucleotide_change"] + "(" + other_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
              variant_substitutions_type.append(other_variant["type"])
              variant_substitutions_nt.append(other_variant["nucleotide_change"])
              variant_substitutions_aa.append(other_variant["protein_change"] if other_variant["protein_change"] != "" else "NA")
              depth.append(other_variant["depth"])
              frequency.append(other_variant["freq"])
              resistance.append(other_variant["gene_associated_drugs"][0])
              rule.append("Uncertain significance based on expert rule") # TODO: keep this? OR change to "Resistant based on expert rule"?
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
            if other_variant["gene"] == "rpoB":  # Expert rule: in case mutation occurs between codons 426 and 452 of rpoB gene, classify as resistant
              position = get_codon(other_variant["protein_change"])
              sample_id.append("~{samplename}")
              gene_name.append(other_variant["gene"])
              locus_tag.append(other_variant["locus_tag"])  
              #variant_substitutions.append(other_variant["type"] + ":" + other_variant["nucleotide_change"] + "(" + other_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
              variant_substitutions_type.append(other_variant["type"])
              variant_substitutions_nt.append(other_variant["nucleotide_change"])
              variant_substitutions_aa.append(other_variant["protein_change"] if other_variant["protein_change"] != "" else "NA")
              depth.append(other_variant["depth"])
              frequency.append(other_variant["freq"])
              resistance.append(other_variant["gene_associated_drugs"][0])
              if 426 <= position <= 452:  # considered resistant based on expert rule - TODO: Keep all or just the ones within the codons?
                rule.append("Resistant based on expert rule")
              else:
                rule.append("Uncertain significance based on expert rule")
            
        with open("tbprofiler_laboratorian_report.csv", "wt") as report_fh:
          report_fh.write("sample_id,tbprofiler_gene_name,tbprofiler_locus_tag,tbprofiler_variant_substitution_type,tbprofiler_variant_substitution_nt,tbprofiler_variant_substitution_aa,confidence,antimicrobial,depth,frequency,read_support,rationale,warning\n")
          for i in range(0, len(gene_name)):
            if not depth[i]:  # for cases when depth is null, it gets converted to 0
              depth[i] = 0
            warning = "Low depth coverage" if  depth[i] < int('~{min_depth}') else "" # warning when coverage is lower than the defined 'min_depth' times
            report_fh.write(sample_id[i] + ',' + gene_name[i] + ',' + locus_tag[i] + ',' + variant_substitutions_type[i] + ',' + variant_substitutions_nt[i] + ',' + variant_substitutions_aa[i] + ',' + confidence[i] + ',' + resistance[i] + ',' + str(depth[i]) + ',' + str(frequency[i]) + ',' + str(int(depth[i]*frequency[i])) + ',' + rule[i] + ',' + warning +'\n')


    def parse_json_mutations(json_file):
      """
      Function to parse the TBProfiler json file and store the found
      mutations into a dictionary
      """
      mutations_dict = {}

      with open(json_file) as js_fh:
        results_json = json.load(js_fh)
        for dr_variant in results_json["dr_variants"]:  # reported mutation by tb-profiler, all confering resistance
          name = dr_variant["gene"]
          substitution =str(dr_variant["nucleotide_change"] + "(" + dr_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
          if name not in mutations_dict.keys():
            mutations_dict[name] = substitution
          else:
            mutations_dict[name] = mutations_dict[name] + ';' + substitution
        for other_variant in results_json["other_variants"]:  # mutations not reported by tb-profiler
          if other_variant["type"] != "synonymous_variant":  # report all non-synonymous mutations
              name = other_variant["gene"]
              substitution =str(other_variant["nucleotide_change"] + "(" + other_variant["protein_change"] + ")")  # mutation_type:nt_sub(aa_sub)
              if name not in mutations_dict.keys():
                mutations_dict[name] = substitution
              else:
                mutations_dict[name] = mutations_dict[name] + ';' + substitution

      return mutations_dict


    def rank_annotation(annotation):
      """
      This function recieves tbprofiler WHO annotation and ranks it based on resistance,
      with 1 being the most resistant category and 4 the least.
      """
      if annotation == "Assoc w R":
        return 1
      elif annotation == "Assoc w R - interim":
        return 2
      elif annotation == "Uncertain significance":
        return 3
      else:
        return 4


    def translate(annotation, drug):
      """
      This function takes the annotation of resistance by tbprofiler and the target drug
      and returns the LIMS' report file appropriate annotation.
      """
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


    def decipher(annotation):
      """
      This function takes the annotation of resistance by TBProfiler and 
      returns simple R (resistant) and S (susceptible) notations
      """
      if annotation == "Not assoc w R":
        return "S"
      elif annotation == "Uncertain significance":
        return "S"
      elif annotation == "Assoc w R - interim":
        return "R"
      elif annotation == "Assoc w R":
        return "R"
      else:
        return "S"
    

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
                  resistance_dict[drug] = resistance
                else:
                  if rank_annotation(resistance_dict[drug]) < rank_annotation(resistance):
                    resistance_dict[drug] = resistance
      return resistance_dict
    

    def get_lineage(json_file):
      with open(json_file) as js_fh:
        results_json = json.load(js_fh)
        if results_json["main_lin"] == "":
          return "DNA of M. tuberculosis complex not detected"
        elif results_json["main_lin"] == "M.bovis":
          return "DNA of M. tuberculosis complex detected (M. bovis)"
        else:
          return "DNA of M. tuberculosis complex detected (not M. bovis)"

    def parse_json_lims_report(json_file):
      """
      This function recieves the tbprofiler output json file and
      writes the LIMS report that includes the following information
      per sample: 
        - MDL sample accession numbers: includes sample name
        - M_DST_A01_ID - includes lineage
        - The set of information in gene_dict dictionary with target drug resistance information
        in layman's terms, and the mutations responsible for the predicted phenotype 
      """
    
      lineage = get_lineage("~{json}")
      mutations = parse_json_mutations("~{json}")
      resistance = parse_json_resistance("~{json}")
      df_lims = pd.DataFrame({"MDL sample accession numbers":"~{samplename}", "M_DST_A01_ID": lineage},index=[0])

      for antimicrobial, genes in gene_dict.items():
        if antimicrobial_dict[antimicrobial] in resistance.keys():
          df_lims[antimicrobial] = translate(resistance[antimicrobial_dict[antimicrobial]], antimicrobial_dict[antimicrobial])
        else:
          df_lims[antimicrobial] = "No resistance to {} detected".format(antimicrobial_dict[antimicrobial])
        for gene_name, gene_id in genes.items():
          if gene_name in mutations.keys():
            df_lims[gene_id] = mutations[gene_name]
          else:
            df_lims[gene_id] = "No mutations detected"
    
      df_lims.to_csv("tbprofiler_lims_report.csv", index=False)
    
    def parse_laboratorian_report(mutations):
      """
      This function parses the laboratorian report for a given mutation and 
      returns the associated WHO annotation and frequency
      """
      df = pd.read_csv('tbprofiler_laboratorian_report.csv', dtype="str")
      confidences = []
      frequencies = []
      for mutation in mutations.split(';'):
        confidence = df.loc[df['tbprofiler_variant_substitutions'] == mutation, 'confidence']
        frequency = df.loc[df['tbprofiler_variant_substitutions'] == mutation, 'frequency']
        if len(confidence) > 0:
          confidences.append(confidence.array[0])
        else:
          confidences.append("No annotation")
        if len(frequency) > 0:
          frequencies.append(frequency.array[0])
        else:
          frequencies.append("1")
      return confidences, frequencies

    def parse_json_looker_report(json_file):
      """
      This function recieves the tbprofiler output json file and
      writes the Looker report that includes the following information
      per sample: 
        - sample_id: includes sample name
        - for each antimicrobial, indication if its resistant (R) or susceptible (S)
      """
      resistance = parse_json_resistance("~{json}")
      print(resistance)
      df_looker = pd.DataFrame({"sample_id":"~{samplename}", "output_seq_method_type": "~{output_seq_method_type}"},index=[0])

      for antimicrobial in antimicrobial_list:
        print(antimicrobial)
        if antimicrobial in resistance.keys():
          df_looker[antimicrobial] = decipher(resistance[antimicrobial])
        else:
          df_looker[antimicrobial] = "S"
    
      df_looker.to_csv("tbprofiler_looker.csv", index=False)


    ### Report Generation ###

    # Laboratorian report generation
    parse_json_lab_report("~{json}")

    # LIMS report generation
    parse_json_lims_report("~{json}")

    # LOOKER report generation
    parse_json_looker_report("~{json}")

    CODE
  >>>
  output {
    File tbprofiler_looker_csv = "tbprofiler_looker.csv"
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