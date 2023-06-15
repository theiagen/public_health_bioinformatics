version 1.0

task tbprofiler_output_parsing {
  input {
    File json
    String output_seq_method_type
    String operator
    String samplename
    Int min_depth = 10
  }
  command <<<
    python3 <<CODE
    import csv
    import json
    import re
    import pandas as pd
    import datetime

    ## Lookup Data Structures ##

    # lookup dictionary - antimicrobial code to name
    antimicrobial_dict = {"M_DST_B01_INH": "isoniazid", "M_DST_C01_ETO": "ethionamide",
                          "M_DST_D01_RIF": "rifampicin", "M_DST_E01_PZA": "pyrazinamide",
                          "M_DST_F01_EMB": "ethambutol","M_DST_H01_AMK": "amikacin", 
                          "M_DST_I01_KAN": "kanamycin","M_DST_J01_CAP": "capreomycin", 
                          "M_DST_K01_MFX": "moxifloxacin","M_DST_L01_LFX": "levofloxacin", 
                          "M_DST_M01_BDQ": "bedaquiline","M_DST_N01_CFZ": "clofazimine", 
                          "M_DST_o01_LZD": "linezolid" 
                         }
    
    # Lookup list - antimicrobials
    antimicrobial_list = ["isoniazid", "ethionamide", "rifampicin", "pyrazinamide", "ethambutol",
                          "streptomycin", "amikacin", "kanamycin", "capreomycin", "moxifloxacin",
                          "levofloxacin", "bedaquiline", "clofazimine", "linezolid"
                         ]

    # lookup dictionary - gene name to gene column name
    gene_dict = {"M_DST_B01_INH": {"katG": "M_DST_B02_katG", "fabG1": "M_DST_B03_fabG1", 
                                   "inhA": "M_DST_B04_inhA"},
                 "M_DST_C01_ETO": {"ethA": "M_DST_C02_ethA", "fabG1": "M_DST_C03_fabG1",
                                   "inhA": "M_DST_C04_inhA"},
                 "M_DST_D01_RIF": {"rpoB": "M_DST_D02_rpoB"},
                 "M_DST_E01_PZA": {"pncA": "M_DST_E02_pncA"},
                 "M_DST_F01_EMB": {"embA": "M_DST_F02_embA", "embB": "M_DST_F03_embB"},
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

    # lookup dictionary - gene to resistance (https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv)
    gene_to_resistance = {"ahpC":["isoniazid"], "ald":["cycloserine"], "alr": ["cycloserine"],
                          "ddn": ["delamanid"], "eis": ["amikacin", "kanamycin"], "embA": ["ethambutol"],
                          "embB": ["ethambutol"], "embC": ["ethambutol"], "embR": ["ethambutol"],
                          "ethA": ["ethionamide"], "ethR": ["ethionamide"], "fabG1": ["ethionamide", "isoniazid"],
                          "fbiA": ["delamanid"], "fgd1": ["delamanid"], "folC": ["para-aminosalicylic_acid"],
                          "gid": ["streptomycin"], "gyrA": ["ciprofloxacin", "fluoroquinolones", "levofloxacin",
                          "moxifloxacin", "ofloxacin"], "gyrB": ["moxifloxacin","ciprofloxacin","fluoroquinolones",
                          "levofloxacin", "ofloxacin"], "inhA": ["ethionamide", "isoniazid"],
                          "kasA": ["isoniazid"], "katG": ["isoniazid"], "panD": ["pyrazinamide"],
                          "pncA": ["pyrazinamide"], "ribD": ["para-aminosalicylic_acid"], "rplC": ["linezolid"],
                          "rpoB": ["rifampicin"], "rpoC": ["rifampicin"], "rpsA": ["pyrazinamide"],
                          "rpsL": ["streptomycin"], "rrl": ["linezolid"], "rrs": ["streptomycin", "amikacin",
                          "aminoglycosides", "capreomycin", "kanamycin"], "Rv0678": ["bedaquiline",
                          "clofazimine"], "thyA": ["para-aminosalicylic_acid"], "thyX": ["para-aminosalicylic_acid"],
                          "tlyA": ["capreomycin"]
                         }

    # lookup dictionary - gene to locus tag (https://github.com/jodyphelan/TBProfiler/blob/master/db/tbdb.bed)
    gene_to_locus_tag = {"ahpC":"Rv2428", "ald":"Rv2780", "alr": "Rv3423c",
                          "ddn": "Rv3547", "eis": "Rv2416c", "embA": "Rv3794",
                          "embB": "Rv3795", "embC": "Rv3793", "embR": "Rv1267c",
                          "ethA": "Rv3854c", "ethR": "Rv3855", "fabG1": "Rv1483",
                          "fbiA": "Rv3261", "fgd1": "Rv0407", "folC": "Rv2447c",
                          "gid": "Rv3919c", "gyrA": "Rv0006", "gyrB": "Rv0005",
                          "inhA": "Rv1484","kasA": "Rv2245", "katG": "Rv1908c", 
                          "panD": "Rv3601c", "pncA": "Rv2043c", "ribD": "Rv2671",
                          "rplC": "Rv0701", "rpoB": "Rv0667", "rpoC": "Rv0668", 
                          "rpsA": "Rv1630", "rpsL": "Rv0682", "rrl": "EBG00000313339", 
                          "rrs": "EBG00000313325", "Rv0678": "Rv0678", "thyA": "Rv2764c",
                          "thyX": "Rv2754c", "tlyA": "Rv1694"
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

    def decipher_looker(annotation):
      """
      This function takes the annotation of resistance by TBProfiler and 
      returns simple R (resistant), U (uncertain) and S (susceptible) notations
      """
      if annotation == "Not assoc w R":
        return "S"
      elif annotation == "Uncertain significance":
        return "U"
      elif annotation == "Assoc w R - interim":
        return "R-interim"
      elif annotation == "Assoc w R":
        return "R"
      else:
        return "S" # missing S-interim?
    
    def decipher_MDL(annotation):
      """
      This function takes the annotation of resistance by TBProfiler and 
      returns simple R (resistant), U (uncertain) and S (susceptible) notations
      """
      if annotation == "Assoc w R - interim":
        return "U"
      elif annotation == "Assoc w R":
        return "R"
      else:
        return "S"

    def get_lineage(json_file):
      with open(json_file) as js_fh:
        results_json = json.load(js_fh)
        if results_json["main_lin"] == "":
          return "DNA of M. tuberculosis complex not detected"
        elif results_json["main_lin"] == "M.bovis":
          return "DNA of M. tuberculosis complex detected (M. bovis)"
        else:
          return "DNA of M. tuberculosis complex detected (not M. bovis)"

    def parse_json_mutations(json_file):
      """
      Function to parse the TBProfiler json file and store the found
      mutations into a dictionary - LIMS Report
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
        return "The detected genetic determinant(s) have uncertain significance, resistance to {} cannot be ruled out".format(drug)
      elif annotation == "Assoc w R - interim":
        return "The detected genetic determinant(s) have uncertain significance, resistance to {} cannot be ruled out".format(drug)
      elif annotation == "Assoc w R":
        return "Genetic determinant(s) associated with resistance to {} detected".format(drug)
      else:
        return "No resistance to {} detected".format(drug)

    def parse_json_resistance(json_file):
      """
      TODO - Comment
      """
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

    ## Main Parsing Functions ## 

    def parse_json_lab_report(json_file):
      """
      This function recieved the tbprofiler output json file and
      writes the Laboratorian report that includes the following information
      per mutation:
        - sample_id: inclides sample name
        - tbprofiler_gene_name: gene name
        - tbprofiler_locus_tag: locus tag
        - tbprofiler_variant_substitution_type: variant substitution type (missense_variant, upstream_gene_variant...)
        - tbprofiler_variant_substitution_nt: nucleotide substitution (c.1349C>G)
        - tbprofiler_variant_substitution_aa: aminoacid substitution (p.Ser450Trp)
        - confidence: tbprofiler annotation regarding resistance (Not assoc w R, Uncertain significance...)
        - antimicrobial: antimicrobial the mutation is confering resistance to (streptomycin, rifampicin...)
        - looker_interpretation: interpretation of resistance for Looker report (R, S, U, R-interim)
        - mdl_interpretation: MDL interpretation of resistance (R,S,U)
        - depth: depth of coverage at the mutation site (100)
        - frequency: frequency of mutation at the site (1)
        - read_support: number of reads supporting the mutation (100, depth*frequency)
        - rationale: rationale for resistance calling (WHO classification, Expert rule)
        - warning: column reserved for warnings such as low depth of coverage
      """

      df_laboratorian = pd.DataFrame(columns = ["sample_id","tbprofiler_gene_name","tbprofiler_locus_tag",
                                                "tbprofiler_variant_substitution_type","tbprofiler_variant_substitution_nt",
                                                "tbprofiler_variant_substitution_aa","confidence","antimicrobial",
                                                "looker_interpretation","mdl_interpretation","depth","frequency",
                                                "read_support","rationale","warning"])
      
      row_list = []
      genes_reported = []

      with open(json_file) as results_json_fh:
        results_json = json.load(results_json_fh)

        # reported mutation by tb-profiler, all confering resistance by WHO criteria
        for dr_variant in results_json["dr_variants"]: 
          if "annotation" in dr_variant:
            for annotation in dr_variant["annotation"]:
              drug = annotation["drug"]
              annotation_who = "No WHO annotation" if annotation["who_confidence"] == "" else annotation["who_confidence"]
              row = {}
              row["sample_id"] = "~{samplename}"
              row["tbprofiler_gene_name"] = dr_variant["gene"]
              row["tbprofiler_locus_tag"] = dr_variant["locus_tag"]
              row["tbprofiler_variant_substitution_type"] = dr_variant["type"]
              row["tbprofiler_variant_substitution_nt"] = dr_variant["nucleotide_change"]
              row["tbprofiler_variant_substitution_aa"] = dr_variant["protein_change"] if dr_variant["protein_change"] != "" else "NA"
              row["confidence"] = annotation_who
              row["antimicrobial"] = drug
              row["looker_interpretation"] = decipher_looker(row["confidence"])
              row["mdl_interpretation"] = decipher_MDL(row["confidence"])
              row["depth"] = int(dr_variant["depth"] or 0)
              row["frequency"] = dr_variant["freq"]
              row["read_support"] = row["depth"]*row["frequency"] 
              row["rationale"] = "WHO classification"
              row["warning"] = "Low depth coverage" if row["depth"] < int('~{min_depth}') else ""
              genes_reported.append(dr_variant["gene"])
              row_list.append(row)
      
      # mutations not reported by tb-profiler - application of expert rules to determine resistance
      for other_variant in results_json["other_variants"]: 

        # report only mutations that are NOT synonymous
        if other_variant["type"] != "synonymous_variant":

          # Expert rule: mutations in Rv0678, atpE, pepQ, mmpL5, mmpS5, rrl amd rplC
          if other_variant["gene"] == "Rv0678" or other_variant["gene"] == "atpE" or other_variant["gene"] == "pepQ" or other_variant["gene"] == "mmpL5" or other_variant["gene"] == "mmpS5" or other_variant["gene"] == "rrl" or other_variant["gene"] == "rplC":
            if "annotation" in other_variant:
              try:  # sometimes annotation is an empty list
                if other_variant["annotation"][0]["who_confidence"] == "":
                  confidence = "No WHO annotation"
                else:
                  confidence = other_variant["annotation"][0]["who_confidence"]
              except:
                confidence = "No WHO annotation"
            else:
              confidence = "No WHO annotation"
            row = {}
            row["sample_id"] = "~{samplename}"
            row["tbprofiler_gene_name"] = other_variant["gene"]
            row["tbprofiler_locus_tag"] = other_variant["locus_tag"]
            row["tbprofiler_variant_substitution_type"] = other_variant["type"]
            row["tbprofiler_variant_substitution_nt"] = other_variant["nucleotide_change"]
            row["tbprofiler_variant_substitution_aa"] = other_variant["protein_change"] if other_variant["protein_change"] != "" else "NA"
            row["confidence"] = confidence
            row["antimicrobial"] = ",".join(other_variant["gene_associated_drugs"])
            row["looker_interpretation"] = decipher_looker(row["confidence"])
            row["mdl_interpretation"] = decipher_MDL(row["confidence"])
            row["depth"] = int(other_variant["depth"] or 0)
            row["frequency"] = other_variant["freq"]
            row["read_support"] = row["depth"]*row["frequency"]
            row["rationale"] = "Resistant based on expert rule"
            row["warning"] = "Low depth coverage" if row["depth"] < int('~{min_depth}') else ""
            genes_reported.append(other_variant["gene"])
            row_list.append(row)

          # Expert rule: mutations in katG, pncA, ethA or gid, classify as resistant
          if other_variant["gene"] == "katG" or other_variant["gene"] == "pncA" or other_variant["gene"] == "ethA" or other_variant["gene"] == "gid":
            if "annotation" in other_variant:
              try:  # sometimes annotation is an empty list
                if other_variant["annotation"][0]["who_confidence"] == "":
                  confidence = "No WHO annotation"
                else:
                  confidence = other_variant["annotation"][0]["who_confidence"]
              except:
                confidence = "No WHO annotation"
            else:
              confidence = "No WHO annotation"
            row = {}
            row["sample_id"] = "~{samplename}"
            row["tbprofiler_gene_name"] = other_variant["gene"]
            row["tbprofiler_locus_tag"] = other_variant["locus_tag"]
            row["tbprofiler_variant_substitution_type"] = other_variant["type"]
            row["tbprofiler_variant_substitution_nt"] = other_variant["nucleotide_change"]
            row["tbprofiler_variant_substitution_aa"] = other_variant["protein_change"] if other_variant["protein_change"] != "" else "NA"
            row["confidence"] = confidence
            row["antimicrobial"] = ",".join(other_variant["gene_associated_drugs"])
            row["looker_interpretation"] = decipher_looker(row["confidence"])
            row["mdl_interpretation"] = decipher_MDL(row["confidence"])
            row["depth"] = int(other_variant["depth"] or 0)
            row["frequency"] = other_variant["freq"]
            row["read_support"] = row["depth"]*row["frequency"]
            row["rationale"] = "Resistant based on expert rule"
            row["warning"] = "Low depth coverage" if row["depth"] < int('~{min_depth}') else ""
            genes_reported.append(other_variant["gene"])
            row_list.append(row)
          
          # Expert rule: in case mutation occurs between codons 426 and 452 of rpoB gene, classify as resistant
          if other_variant["gene"] == "rpoB": 
            position = get_codon(other_variant["protein_change"])
            row = {}
            row["sample_id"] = "~{samplename}"
            row["tbprofiler_gene_name"] = other_variant["gene"]
            row["tbprofiler_locus_tag"] = other_variant["locus_tag"]
            row["tbprofiler_variant_substitution_type"] = other_variant["type"]
            row["tbprofiler_variant_substitution_nt"] = other_variant["nucleotide_change"]
            row["tbprofiler_variant_substitution_aa"] = other_variant["protein_change"] if other_variant["protein_change"] != "" else "NA"
            row["confidence"] = "No WHO annotation"
            row["antimicrobial"] = ",".join(other_variant["gene_associated_drugs"])
            row["looker_interpretation"] = decipher_looker(row["confidence"])
            row["mdl_interpretation"] = decipher_MDL(row["confidence"])
            row["depth"] = int(other_variant["depth"] or 0)
            row["frequency"] = other_variant["freq"]
            row["read_support"] = row["depth"]*row["frequency"]
            row["rationale"] = "Resistant based on expert rule" if 426 <= position <= 452 else "Uncertain significance based on expert rule"
            row["warning"] = "Low depth coverage" if row["depth"] < int('~{min_depth}') else ""
            genes_reported.append(other_variant["gene"])
            row_list.append(row)
      
      for gene, resistance_list in gene_to_resistance.items():
        for resistance in resistance_list:
          if gene not in genes_reported:
            row = {}
            row["sample_id"] = "~{samplename}"
            row["tbprofiler_gene_name"] = gene
            row["tbprofiler_locus_tag"] = gene_to_locus_tag[gene]
            row["tbprofiler_variant_substitution_type"] = "WT"
            row["tbprofiler_variant_substitution_nt"] = "NA"
            row["tbprofiler_variant_substitution_aa"] = "NA"
            row["confidence"] = "NA"
            row["antimicrobial"] = resistance
            row["looker_interpretation"] = "NA"
            row["mdl_interpretation"] = "NA"
            row["depth"] = "NA"
            row["frequency"] = "NA"
            row["read_support"] = "NA"
            row["rationale"] = "NA"
            row["warning"] = "NA"
            row_list.append(row)
      
      df_laboratorian = df_laboratorian.append(row_list, ignore_index=True)
      df_laboratorian.to_csv("tbprofiler_laboratorian_report.csv", index=False)
    
    def parse_json_lims_report(json_file, formatted_time):
      """
      This function recieves the tbprofiler output json file and
      writes the LIMS report that includes the following information
      per sample: 
        - MDL sample accession numbers: includes sample name
        - M_DST_A01_ID - includes lineage
        - The set of information in gene_dict dictionary with target drug resistance information
        in layman's terms, and the mutations responsible for the predicted phenotype
        - Date of analysis in YYYY-MM-DD HH:SS format
        - Operator information
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

      df_lims["Analysis date"] = formatted_time
      df_lims["Operator"] = "~{operator}"
      df_lims.to_csv("tbprofiler_lims_report.csv", index=False)
    
    def parse_json_looker_report(json_file, current_time):
      """
      This function recieves the tbprofiler output json file and
      writes the Looker report that includes the following information
      per sample: 
        - sample_id: includes sample name
        - for each antimicrobial, indication if its resistant (R) or susceptible (S)
      """
      resistance = parse_json_resistance("~{json}")
      df_looker = pd.DataFrame({"sample_id":"~{samplename}", "output_seq_method_type": "~{output_seq_method_type}"},index=[0])

      for antimicrobial in antimicrobial_list:
        if antimicrobial in resistance.keys():
          df_looker[antimicrobial] = decipher_looker(resistance[antimicrobial])
        else:
          df_looker[antimicrobial] = "S"
      
      df_looker["analysis_date"] = current_time
      df_looker["operator"] = "~{operator}"
    
      df_looker.to_csv("tbprofiler_looker.csv", index=False)

    ### Report Generation ###

    # get timestamp in YYYY-MM-DD HH:MM format
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

    # Laboratorian report generation
    parse_json_lab_report("~{json}")

    # LIMS report generation
    parse_json_lims_report("~{json}", current_time)

    # LOOKER report generation
    parse_json_looker_report("~{json}", current_time)

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