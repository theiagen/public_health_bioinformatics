version 1.0

task tbprofiler_output_parsing {
  input {
    File json
    File gene_coverage
    String output_seq_method_type
    String operator
    String samplename
    Int min_depth = 10
    Int coverage_threshold = 100
  }
  command <<<
    python3 <<CODE
    import csv
    import json
    import re
    import pandas as pd
    import datetime

    ## Lookup Data Structures ##

    # gene coverage as a dictionary
    gene_coverage_dict = pd.read_csv("~{gene_coverage}", delimiter="\t", skip_blank_lines=True).to_dict()
    gene_coverage_dict = gene_coverage_dict["#NOTE: THE VALUES BELOW ASSUME TBPROFILER (H37Rv) REFERENCE GENOME"] # skip first line

    # lookup dictionary - antimicrobial code to drug name
    antimicrobial_code_to_drug_name = {"M_DST_B01_INH": "isoniazid", "M_DST_C01_ETO": "ethionamide",
                          "M_DST_D01_RIF": "rifampicin", "M_DST_E01_PZA": "pyrazinamide",
                          "M_DST_F01_EMB": "ethambutol","M_DST_H01_AMK": "amikacin", 
                          "M_DST_I01_KAN": "kanamycin","M_DST_J01_CAP": "capreomycin", 
                          "M_DST_K01_MFX": "moxifloxacin","M_DST_L01_LFX": "levofloxacin", 
                          "M_DST_M01_BDQ": "bedaquiline","M_DST_N01_CFZ": "clofazimine", 
                          "M_DST_o01_LZD": "linezolid" 
                         }
    
    # Lookup list - antimicrobial drug names
    antimicrobial_drug_name_list = ["isoniazid", "ethionamide", "rifampicin", "pyrazinamide", "ethambutol",
                          "streptomycin", "amikacin", "kanamycin", "capreomycin", "moxifloxacin",
                          "levofloxacin", "bedaquiline", "clofazimine", "linezolid"
                         ]

    # lookup dictionary - antimicrobial code to gene names to antimicrobial code column names
    antimicrobial_code_to_genes = {"M_DST_B01_INH": {"katG": "M_DST_B02_katG", "fabG1": "M_DST_B03_fabG1", 
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

    # lookup dictionary - gene to antimicrobial drug name (https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv)
    gene_to_antimicrobial_drug_name = {"ahpC":["isoniazid"], "ald":["cycloserine"], "alr": ["cycloserine"],
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
    
    # lookup dictionary - antimicrobial drug name to gene name  (https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv)
    antimicrobial_drug_name_to_gene = { "isoniazid": ["ahpC", "fabG1", "inhA", "kasA", "katG"], 
                                        "ethionamide": ["ethA", "ethR", "fabG1", "inhA"], 
                                        "rifampicin": ["rpoB", "rpoC"], 
                                        "pyrazinamide": ["panD", "pncA", "rpsA"], 
                                        "ethambutol": ["embA", "embB", "embC", "embR"], 
                                        "streptomycin": ["gid", "rpsL", "rrs"], 
                                        "amikacin": ["eis", "rrs"], 
                                        "kanamycin": ["eis", "rrs"], 
                                        "capreomycin": ["rrs", "tlyA"], 
                                        "moxifloxacin": ["gyrA", "gyrB"],
                                        "levofloxacin": ["gyrA", "gyrB"], 
                                        "bedaquiline": ["Rv0678"], 
                                        "clofazimine": ["Rv0678"], 
                                        "linezolid": ["rplC", "rrl"]
                                      }
    
    # lookup dictionary - gene to tier
    gene_to_tier = {"ahpC": "Tier 1", "inhA": "Tier 1", "katG": "Tier 1", "rpoB": "Tier 1", "embA": "Tier 1", 
                    "embB": "Tier 1", "embC": "Tier 1", "pncA": "Tier 1", "clpC1": "Tier 1", "panD": "Tier 1", 
                    "gyrA": "Tier 1", "gyrB": "Tier 1", "pepQ": "Tier 1", "Rv0678": "Tier 1", "mmpL5": "Tier 1", 
                    "mmpS5": "Tier 1", "atpE": "Tier 1", "rplC": "Tier 1", "rrl": "Tier 1", "fgd1": "Tier 1", 
                    "ddn": "Tier 1", "fbiA": "Tier 1", "fbiB": "Tier 1", "fbiC": "Tier 1", "Rv2983": "Tier 1", 
                    "rrs": "Tier 1", "eis": "Tier 1", "whiB7": "Tier 1", "rpsL": "Tier 1", "gid": "Tier 1", 
                    "Rv1258c": "Tier 1", "ethA": "Tier 1", "tlyA": "Tier 1", "mshA": "Tier 2", "ndh": "Tier 2", 
                    "Rv2752c": "Tier 2", "rpoA": "Tier 2", "rpoC": "Tier 2", "embR": "Tier 2", "ubiA": "Tier 2", 
                    "PPE35": "Tier 2", "Rv3236c": "Tier 2", "Rv1979c": "Tier 2", "whiB6": "Tier 2", "ccsA": "Tier 2", 
                    "fprA": "Tier 2", "aftB": "Tier 2", "ethR": "Tier 2", "Rv3083": "Tier 2"
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

    gene_list_option_1 = ["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5", "rrl", "rplC"] # Rv0678 is mmpR
    gene_list_option_2 = ["katG", "pncA", "ethA", "gid", "rpoB"]
    gene_list_combined = ["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5", "rrl", "rplC", "katG", "pncA", "ethA", "gid", "rpoB"]
    low_depth_of_coverage_list = []
    ## Auxiliary Functions ##

    def get_position(protein_mut):
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

    def annotation_to_looker(annotation):
      """
      1.1: This function takes the WHO-annotation of resistance by TBProfiler and 
      returns simple R (resistant), U (uncertain) and S (susceptible) notations for Looker
      """
      if annotation == "Assoc w R":
        return "R"
      elif annotation == "Assoc w R - interim":
        return "R-Interim"
      elif annotation == "Uncertain significance":
        return "U"
      elif annotation == "Not assoc w R":
        return "S"
      elif annotation == "Not assoc w R - Interim":
        return "S-Interim"
      else:
        return "S"
    
    def annotation_to_MDL(annotation, gene = ""):
      """
      This function takes the annotation of resistance by TBProfiler and 
      returns simple R (resistant), U (uncertain) and S (susceptible) notations for MDL
      """
      if annotation == "Assoc w R":
        return "R"
      elif annotation == "Assoc w R - interim":
        return "U"
      elif annotation == "Uncertain significance":
        if gene in gene_list_option_1:
          return "U"
        else:
          return "S"
      elif annotation == "Not assoc w R":
        return "S"
      elif annotation == "Not assoc w R - Interim":
        if gene in gene_list_option_1:
          return "U"
        else:
          return "S"
      else:
        return "S"

    def apply_expert_rules(nucleotide_change, protein_change, gene, substitution_type, interpretation_destination):
      """
      Apply rules 1-3
      """

      position_nt = get_position(nucleotide_change)
      position_aa = get_position(protein_change)

      if gene in ["Rv0678", "atpE", "pepQ", "rplC", "mmpL5", "mmpS5"]: # apply expert rules 1.2
        if gene == "Rv0678" and (position_nt >= -84 and position_nt <= -1): # promoter regions
          return "U"
        elif gene == "atpE" and (position_nt >= -48 and position_nt <= -1):
          return "U"
        elif gene == "pepQ" and (position_nt >= -33 and position_nt <= -1):
          return "U"
        elif gene == "rplC" and (position_nt >= -18 and position_nt <= -1):
          return "U"
        else: # apply expert rules 1.2
          if not any(non_ORF in nucleotide_change for non_ORF in ["+", "-", "*"]) or nucleotide_change.endswith("*"): 
          # if a position includes either +, *, or - it's not in the ORF 
          #  UNLESS the * is at the end which means its a premature stop codon
            if substitution_type != "synonymous_variant":
              return "U"
            else:
              return "S"

      elif gene == "rrl": # apply expert rules 1.2
        if (position_nt >= 2003 and position_nt <= 2367) or (position_nt >= 2449 and position_nt <= 3056):
          return "U"
        else:
          return "S" if interpretation_destination == "MDL" else "U"

      elif gene in ["katG", "pncA", "ethA", "gid"]: # apply expert rules 2.2.1
        if any(indel_or_stop in nucleotide_change for indel_or_stop in ["del", "ins", "fs", "delins", "_"]) or nucleotide_change.endswith("*"):
          return "U"
        else:
            if substitution_type != "synonymous_variant":
              return "S" if interpretation_destination == "MDL" else "U"
            else:
              return "S"

      elif gene == "rpoB": # apply expert rules 2.2.2
        if (position_aa >= 426 and position_aa <= 452):
            if substitution_type != "synonymous_variant":
              return "R"
            else:
              return "S"   
        else:
            if substitution_type != "synonymous_variant":
              return "S" if interpretation_destination == "MDL" else "U"
            else:
              return "S"

      elif gene not in gene_list_combined: # NOT AN EXPERT RULE: 3.2
        if substitution_type != "synonymous_variant":
          return "Snoexpert" if interpretation_destination == "MDL" else "Unoexpert"
        else:
          return "Snoexpert"

      return ""

    def remove_no_expert(row):
      if "noexpert" in row["looker_interpretation"]:
        interpretation = row["looker_interpretation"]
        row["looker_interpretation"] = interpretation.replace("noexpert", "")
        row["rationale"] = "No WHO annotation or expert rule"
      
      if "noexpert" in row["mdl_interpretation"]:
        interpretation = row["mdl_interpretation"]
        row["mdl_interpretation"] = interpretation.replace("noexpert", "")
        row["rationale"] = "No WHO annotation or expert rule"
      
      return row

    def get_lineage_LIMS(json_file):
      """
      This function returns the lineage in English for the LIMS report
      """
      with open(json_file) as js_fh:
        results_json = json.load(js_fh)
        if results_json["main_lin"] == "":
          return "No main lineage could be identified"
        elif "La1" in results_json["main_lin"]:
          return "DNA of M. tuberculosis complex detected (M. bovis)"
        elif "lineage" in results_json["main_lin"]:
          return "DNA of M. tuberculosis complex detected (M. tb)"
        else:
          return "DNA of M. tuberculosis complex detected (not M. bovis and not M. tb)"
    
    def get_lineage_and_ID_Looker(json_file):
      """
      This function returns the Lineage and ID for Looker from the sublineage field
      """
      with open(json_file) as js_fh:
        results_json = json.load(js_fh)

        sublineage = results_json["sublin"]
        lineage = "NA"

        if "lineage" in sublineage:
          lineage = sublineage
          ID = "MTBC, not M. bovis"
        elif "BCG" in sublineage:
          ID = "M. bovis BCG"
        elif "bovis" in sublineage and "BCG" not in sublineage:
          ID = "M. bovis, not BCG"  
        elif sublineage == "":
          ID = "NA"
        else:
          ID = sublineage

      return lineage, ID

    def parse_json_mutations_for_LIMS(json_file):
      """
      Function to parse the TBProfiler json file and store the found
      mutations into a dictionary for the LIMS Report
      """
      mutations_dict = {}

      with open(json_file) as js_fh:
        results_json = json.load(js_fh)
        for dr_variant in results_json["dr_variants"]:  # reported mutation by tb-profiler, all confering resistance
          name = dr_variant["gene"]
          substitution = str(dr_variant["nucleotide_change"] + " (" + dr_variant["protein_change"] + ")")  # mutation_type:nt_sub (aa_sub)
          if name not in mutations_dict.keys():
            mutations_dict[name] = substitution
          else:
            mutations_dict[name] = mutations_dict[name] + '; ' + substitution
        for other_variant in results_json["other_variants"]:  # mutations not reported by tb-profiler
          aa_position = get_position(other_variant["protein_change"])
          if other_variant["type"] != "synonymous_variant" or (other_variant["gene"] == "rpoB" and (aa_position >= 426 and aa_position <= 452)):  # report all non-synonymous mutations, unless rpoB RRDR
            name = other_variant["gene"]
            substitution = str(other_variant["nucleotide_change"] + " (" + other_variant["protein_change"] + ")")  # mutation_type:nt_sub (aa_sub)
            if other_variant["type"] == "synonymous_variant": # rpoB RRDR
              substituion = subtitution + " [synonymous]"
            if name not in mutations_dict.keys():
              mutations_dict[name] = substitution
            else:
              mutations_dict[name] = mutations_dict[name] + '; ' + substitution
      return mutations_dict

    def rank_annotation(annotation):
      """
      This function recieves tbprofiler WHO annotation and ranks it based on resistance,
      with 4 being the most resistant category and 1 the least.
      """
      if annotation == "Assoc w R":
        return 4
      elif annotation == "Assoc w R - interim":
        return 3
      elif annotation == "Uncertain significance":
        return 2
      else:
        return 1

    def annotation_to_LIMS(annotation, drug):
      """
      This function takes the annotation of resistance by tbprofiler and the target drug
      and returns the LIMS' report file appropriate annotation.
      """
      if annotation == "Assoc w R":
        return "Genetic determinant(s) associated with resistance to {} detected".format(drug)
      elif (annotation == "Assoc w R - interim") or (annotation == "Uncertain significance"):
        return "The detected genetic determinant(s) have uncertain significance, resistance to {} cannot be ruled out".format(drug)
      else: # "Not assoc w R" and "Not assoc w R - Interim" and anything else
        return "No genetic determinants associated with resistance to {} detected".format(drug)

    def parse_json_resistance(json_file):
      """
      This function parses the tbprofiler json report and returns a resistance dictionary
      containing the WHO resistance annotation for each antimicrobial drug for LIMS and Looker. 
      The annotation corresponds to the highest ranked one regarding severity (R > I > S)
      """
      resistance_dict = {}

      with open(json_file) as js_fh:
        results_json = json.load(js_fh)

        for dr_variant in results_json["dr_variants"]: # mutation 
          if dr_variant["type"] != "synonymous_variant":  # report all non-synonymous mutations
            if "annotation" in dr_variant.keys(): # if an annotation is present,
              for annotation in dr_variant["annotation"]: # iterate through them
                drug = annotation["drug"]
                who_annotation = annotation["who_confidence"]
                if drug not in resistance_dict.keys():
                  resistance_dict[drug] = who_annotation
                else: # if the drug has already been seen in either the same variant or a different one,
                  if rank_annotation(resistance_dict[drug]) < rank_annotation(who_annotation): # if current annotation indicates higher severity than any previous annotation,
                    resistance_dict[drug] = who_annotation # overwrite with the who_annotation
          
        for other_variant in results_json["other_variants"]:
          if other_variant["type"] != "synonymous_variant":  # report all non-synonymous mutations
            if "annotation" in other_variant.keys(): # if an annotation is present,
              for annotation in other_variant["annotation"]: # iterate through them
                drug = annotation["drug"]
                who_annotation = annotation["who_confidence"]
                if drug not in resistance_dict.keys():
                  resistance_dict[drug] = who_annotation
                else: # if the drug has already been seen in either the same variant or a different one,
                  if rank_annotation(resistance_dict[drug]) < rank_annotation(who_annotation): # if current annotation indicates higher severity than any previous annotation,
                    resistance_dict[drug] = who_annotation # overwrite with the who_annotation

      return resistance_dict

    def variant_to_row(variant):
      """
      This function recieved a variants dictionary and returns a row dictionary containing the
      basic information for it to be filled in the laboratorian report
      """
      print(variant)
      row = {}
      row["sample_id"] = "~{samplename}"
      row["tbprofiler_gene_name"] = variant["gene"]
      if variant["gene"] in gene_to_tier.keys():
        row["gene_tier"] = gene_to_tier[variant["gene"]]
      else:
        row["gene_tier"] = "NA"
      row["tbprofiler_locus_tag"] = variant["locus_tag"]
      row["tbprofiler_variant_substitution_type"] = variant["type"]
      row["tbprofiler_variant_substitution_nt"] = variant["nucleotide_change"]
      row["tbprofiler_variant_substitution_aa"] = variant["protein_change"] if variant["protein_change"] != "" else "NA"
      row["depth"] = int(variant["depth"] or 0)
      row["frequency"] = variant["freq"]
      row["read_support"] = row["depth"]*row["frequency"]
      row["warning"] = ""
      print(gene_coverage_dict[variant["gene"]])
      if row["depth"] < int('~{min_depth}') or float(gene_coverage_dict[variant["gene"]]) < ~{coverage_threshold}:
        row["warning"] = "Insufficient coverage in locus"
        if "del" in variant["nucleotide_change"]:
          row["warning"] = "Insufficient coverage in locus (deletion identified)"
        else:
          low_depth_of_coverage_list.append(variant["gene"])

      return row

    ## Main Parsing Functions ## 

    def parse_json_lab_report(json_file):
      """
      This function recieves the tbprofiler output json file and
      writes the Laboratorian report that includes the following information
      per mutation:
        - sample_id: includes sample name
        - tbprofiler_gene_name: gene name
        - tbprofiler_locus_tag: locus tag
        - tbprofiler_variant_substitution_type: variant substitution type (missense_variant, upstream_gene_variant...)
        - tbprofiler_variant_substitution_nt: nucleotide substitution (c.1349C>G)
        - tbprofiler_variant_substitution_aa: aminoacid substitution (p.Ser450Trp)
        - confidence: tbprofiler annotation regarding resistance (Not assoc w R, Uncertain significance...)
        - antimicrobial: antimicrobial the mutation is confering resistance to (streptomycin, rifampicin...)
        - looker_interpretation: interpretation of resistance for Looker report (R, R-interim, U, S, S-interim)
        - mdl_interpretation: MDL interpretation of resistance (R, S, U)
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
          genes_reported.append(dr_variant["gene"])
          if "annotation" in dr_variant:
            # empty dictionary in the case where multiple drugs are identified
            drugs_to_row = {}
            # case: annotation field is present but it's an empty list, expert rule is applied directly 
            if len(dr_variant["annotation"]) == 0:
              row = variant_to_row(dr_variant)
              row["confidence"] = "No WHO annotation"
              row["looker_interpretation"] = apply_expert_rules(dr_variant["nucleotide_change"], dr_variant["protein_change"], dr_variant["gene"], dr_variant["type"], "looker")
              row["mdl_interpretation"] = apply_expert_rules(dr_variant["nucleotide_change"], dr_variant["protein_change"], dr_variant["gene"], dr_variant["type"], "MDL")
              row["rationale"] = "Expert rule applied"

              # iterate through any gene-associated drugs
              for identified_drug in dr_variant["gene_associated_drugs"]:
                if identified_drug not in drugs_to_row:
                  drugs_to_row[identified_drug] = {"other_variant": dr_variant, "who_confidence": row["confidence"], "drug": identified_drug, "nucleotide_change": dr_variant["nucleotide_change"], "protein_change": dr_variant["protein_change"], "gene": dr_variant["gene"], "type": dr_variant["type"]}
                elif rank_annotation(drugs_to_row[identified_drug]["who_confidence"]) > rank_annotation(row["confidence"]): # overwrite entry with the more severe annotation (higher value) if multiple drugs are present
                  drugs_to_row[identified_drug] = {"other_variant": dr_variant, "who_confidence": row["confidence"], "drug": identified_drug, "nucleotide_change": dr_variant["nucleotide_change"], "protein_change": dr_variant["protein_change"], "gene": dr_variant["gene"], "type": dr_variant["type"]}

            # case: drug confers resistance to multiple drugs - if the same drug shows multiple times in a single mutation, save only the most severe annotation
            for annotation in dr_variant["annotation"]: # iterate thorugh all possible annotations for the variant
              if annotation["drug"] not in drugs_to_row: # if this is the first time a drug is seen, add to dictionary
                drugs_to_row[annotation["drug"]] = {"other_variant": dr_variant, "who_confidence": annotation["who_confidence"], "drug": annotation["drug"], "nucleotide_change": dr_variant["nucleotide_change"], "protein_change": dr_variant["protein_change"], "gene": dr_variant["gene"], "type": dr_variant["type"]}
              elif rank_annotation(drugs_to_row[annotation["drug"]]["who_confidence"]) > rank_annotation(annotation["who_confidence"]): # overwrite entry with the more severe annotation (higher value) if multiple drugs are present
                drugs_to_row[annotation["drug"]] = {"other_variant": dr_variant, "who_confidence": annotation["who_confidence"], "drug": annotation["drug"], "nucleotide_change": dr_variant["nucleotide_change"], "protein_change": dr_variant["protein_change"], "gene": dr_variant["gene"], "type": dr_variant["type"]}
            
            for drug in drugs_to_row:
              print(drug)
              row = variant_to_row(drugs_to_row[drug]["other_variant"])
              print("HERE")
              row["confidence"] = "No WHO annotation" if drugs_to_row[drug]["who_confidence"] == "" else drugs_to_row[drug]["who_confidence"]
              row["antimicrobial"] = drugs_to_row[drug]["drug"]
              row["looker_interpretation"] = annotation_to_looker(row["confidence"])  if row["confidence"] != "No WHO annotation" else apply_expert_rules(drugs_to_row[drug]["nucleotide_change"], drugs_to_row[drug]["protein_change"], drugs_to_row[drug]["gene"], drugs_to_row[drug]["type"], "looker")
              row["mdl_interpretation"] = annotation_to_MDL(row["confidence"], drugs_to_row[drug]["gene"]) if row["confidence"] != "No WHO annotation" else apply_expert_rules(drugs_to_row[drug]["nucleotide_change"], drugs_to_row[drug]["protein_change"], drugs_to_row[drug]["gene"], drugs_to_row[drug]["type"], "MDL")
              row["rationale"] = "WHO classification"  if row["confidence"] != "No WHO annotation" else "Expert rule applied"
              row = remove_no_expert(row)
              row_list.append(row)

          # case: annotation field is not present, expert rule is applied directly
          else:
            for drug in dr_variant["gene_associated_drugs"]:
              row = variant_to_row(dr_variant)
              row["confidence"] = "No WHO annotation"
              row["looker_interpretation"] = apply_expert_rules(dr_variant["nucleotide_change"], dr_variant["protein_change"], dr_variant["gene"], dr_variant["type"], "looker")
              row["mdl_interpretation"] = apply_expert_rules(dr_variant["nucleotide_change"], dr_variant["protein_change"], dr_variant["gene"], dr_variant["type"], "MDL")
              row["rationale"] = "Expert rule applied"
              row["antimicrobial"] = drug
              row = remove_no_expert(row)
              row_list.append(row)
      
        # mutations not reported by tb-profiler
        for other_variant in results_json["other_variants"]: 
          genes_reported.append(other_variant["gene"])
          if "annotation" in other_variant:
            # case: annotation field is present but it's an empty list 
            if len(other_variant["annotation"]) == 0:
              row = variant_to_row(other_variant)
              row["confidence"] = "No WHO annotation"
              row["looker_interpretation"] = apply_expert_rules(other_variant["nucleotide_change"], other_variant["protein_change"], other_variant["gene"], other_variant["type"], "looker")
              row["mdl_interpretation"] = apply_expert_rules(other_variant["nucleotide_change"], other_variant["protein_change"], other_variant["gene"], other_variant["type"], "MDL")
              row["rationale"] = "Expert rule applied"
              row = remove_no_expert(row)
              row_list.append(row)
            # case: drug confers resistance to multiple drugs - if the same drug shows multiple times, save only the most severe annotation
            drugs_to_row = {}
            for annotation in other_variant["annotation"]:
              if annotation["drug"] not in drugs_to_row:
                drugs_to_row[annotation["drug"]] = {"other_variant": other_variant, "who_confidence": annotation["who_confidence"], "drug": annotation["drug"], "nucleotide_change": other_variant["nucleotide_change"],"protein_change": other_variant["protein_change"], "gene": other_variant["gene"], "type": other_variant["type"]}
              elif rank_annotation(drugs_to_row[annotation["drug"]]["who_confidence"]) > rank_annotation(annotation["who_confidence"]):
                drugs_to_row[annotation["drug"]] = {"other_variant": other_variant, "who_confidence": annotation["who_confidence"], "drug": annotation["drug"], "nucleotide_change": other_variant["nucleotide_change"],"protein_change": other_variant["protein_change"], "gene": other_variant["gene"], "type": other_variant["type"]}

            for drug in drugs_to_row:
              row = variant_to_row(drugs_to_row[drug]["other_variant"])
              row["confidence"] = "No WHO annotation" if drugs_to_row[drug]["who_confidence"] == "" else drugs_to_row[drug]["who_confidence"]
              row["antimicrobial"] = drugs_to_row[drug]["drug"]
              row["looker_interpretation"] = annotation_to_looker(row["confidence"])  if row["confidence"] != "No WHO annotation" else apply_expert_rules(drugs_to_row[drug]["nucleotide_change"], drugs_to_row[drug]["protein_change"], drugs_to_row[drug]["gene"], drugs_to_row[drug]["type"], "looker")
              row["mdl_interpretation"] = annotation_to_MDL(row["confidence"], drugs_to_row[drug]["gene"]) if row["confidence"] != "No WHO annotation" else apply_expert_rules(drugs_to_row[drug]["nucleotide_change"],drugs_to_row[drug]["protein_change"], drugs_to_row[drug]["gene"], drugs_to_row[drug]["type"], "MDL")
              row["rationale"] = "WHO classification"  if row["confidence"] != "No WHO annotation" else "Expert rule applied"
              row = remove_no_expert(row)
              row_list.append(row)
          else:
            for drug in other_variant["gene_associated_drugs"]:
              row = variant_to_row(other_variant)
              row["confidence"] = "No WHO annotation"
              row["looker_interpretation"] = apply_expert_rules(other_variant["nucleotide_change"], other_variant["protein_change"], other_variant["gene"], other_variant["type"], "looker")
              row["mdl_interpretation"] = apply_expert_rules(other_variant["nucleotide_change"], other_variant["protein_change"], other_variant["gene"], other_variant["type"], "MDL")
              row["rationale"] = "Expert rule applied"
              row["antimicrobial"] = drug
              row = remove_no_expert(row)
              row_list.append(row)
              
      for gene, antimicrobial_drug_names in gene_to_antimicrobial_drug_name.items():
        for drug_name in antimicrobial_drug_names:
          if gene not in genes_reported:

              row = {}
              row["sample_id"] = "~{samplename}"
              row["tbprofiler_gene_name"] = gene
              row["tbprofiler_locus_tag"] = gene_to_locus_tag[gene]
              row["tbprofiler_variant_substitution_nt"] = "NA"
              row["tbprofiler_variant_substitution_aa"] = "NA"
              if float(gene_coverage_dict[gene]) >= ~{coverage_threshold}:
                row["tbprofiler_variant_substitution_type"] = "WT"
                row["looker_interpretation"] = "S"
                row["mdl_interpretation"] = "WT"
                row["tbprofiler_variant_substitution_nt"] = "WT"
                row["tbprofiler_variant_substitution_aa"] = "WT"
              else:
                row["tbprofiler_variant_substitution_type"] = "Insufficient Coverage"
                row["looker_interpretation"] = "Insufficient Coverage"
                row["mdl_interpretation"] = "Insufficient Coverage"
              if gene in gene_to_tier.keys():
                row["gene_tier"] = gene_to_tier[gene]
              else:
                row["gene_tier"] = "NA"
              row["confidence"] = "NA"
              row["antimicrobial"] = drug_name
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
        - The set of information in antimicrobial_code_to_genes dictionary with target drug resistance information in layman's terms, and the mutations responsible for the predicted phenotype
        - Date of analysis in YYYY-MM-DD HH:SS format
        - Operator information
      """
    
      lineage = get_lineage_LIMS("~{json}")
      mutations = parse_json_mutations_for_LIMS("~{json}")
      resistance_annotation = parse_json_resistance("~{json}")
      df_lims = pd.DataFrame({"MDL sample accession numbers":"~{samplename}", "M_DST_A01_ID": lineage}, index=[0])

      for antimicrobial_code, genes in antimicrobial_code_to_genes.items():
        drug_name = antimicrobial_code_to_drug_name[antimicrobial_code]
        if drug_name in resistance_annotation.keys(): # if the drug has been mentioned in the results file
          df_lims[antimicrobial_code] = annotation_to_LIMS(resistance_annotation[drug_name], drug_name)
        else: # the drug is not in the results file
          df_lims[antimicrobial_code] = "No genetic determinants associated with resistance to {} detected".format(drug_name)
        
        for gene_name, gene_column_code in genes.items(): # iterate through the genes that are associated with resistance to each drug
          if gene_name in mutations.keys(): # if the gene is mentioned in the mutations
            df_lims[gene_column_code] = mutations[gene_name] # put the formatted mutation as the content for the column
            if gene_name == "rpoB": # rule 5.2.1.2
              if df_lims[antimicrobial_code][0] == "No genetic determinants associated with resistance to {} detected".format(drug_name):
                if len(mutations) > 0: # if any mutations are present
                  non_synomynous_count = 0
                  for mutation in mutations: # check if that mutation is synonymous (only output if rpoB RRDR -- see the parse_json_mutations_for_LIMS function)
                    if "synonymous" not in mutation: # if any nonsynymous mutations were identified.
                      non_synomynous_count += 1 # keep the original output for the antimicrobial code
                  if non_synomynous_count == 0: # otherwise, the only synonymous mutations were identified in rpoB RRDR
                    df_lims[antimicrobial_code] = "No genetic determinants associated with resistance to rifampin detected. The detected synonymous mutation(s) do not confer resistance but may result in false-resistance in PCR-based assays targeting the rpoB RRDR."
            
            if gene_name == "rrl": # Rule 5.2.2
              df_lims[gene_column_code] = "" # do not report mutations for rrl
            if df_lims[antimicrobial_code][0] == "No genetic determinants associated with resistance to {} detected".format(drug_name): # if the mutations detected were only "S", 
              df_lims[gene_column_code] = "No high confidence mutations detected"
          if float(gene_coverage_dict[gene_name]) < ~{coverage_threshold}: # HOWEVER, if the coverage is less than the indicated threshold
            df_lims[gene_column_code] = "Insufficient Coverage"
            if "del" in mutations[gene_name]:
              df_lims[gene_column_code] = "Insufficient Coverage (deletion identified)"
            if int(rank_annotation(resistance_annotation[drug_name])) < 4: # in addition, if the indicated annotation for the drug is not resistant (less than 4)
              df_lims[antimicrobial_code] = "Pending Retest"
          else: # the gene is not in the mutations list but has decent coverage
            df_lims[gene_column_code] = "No mutations detected"

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

      lineage, ID = get_lineage_and_ID_Looker("~{json}")
      resistance_annotation = parse_json_resistance("~{json}")
      df_looker = pd.DataFrame({"sample_id":"~{samplename}", "output_seq_method_type": "~{output_seq_method_type}"}, index=[0])

      for antimicrobial_drug in antimicrobial_drug_name_list:
        if antimicrobial_drug in resistance_annotation.keys():
          who_annotation = resistance_annotation[antimicrobial_drug]
          df_looker[antimicrobial_drug] = annotation_to_looker(who_annotation)
          for gene in antimicrobial_drug_name_to_gene[antimicrobial_drug]:
            # indicate warning if any genes failed to achieve 100% coverage_threshold and/or minimum depth  (10x) 
            if annotation_to_looker(who_annotation) != "R" and gene in low_depth_of_coverage_list:
              df_looker[antimicrobial_drug] = "Insufficient coverage for the locus"
          
        else: # the antimicrobial drug was not present in the results
          df_looker[antimicrobial_drug] = "S"
      
      df_looker["lineage"] = lineage 
      df_looker["ID"] = ID
      df_looker["analysis_date"] = current_time
      df_looker["operator"] = "~{operator}"
    
      df_looker.to_csv("tbprofiler_looker.csv", index=False)

    def regenerate_coverage_report(gene_coverage_dict, lab_report):
      """
      This function takes the information generated previously regarding
      breadth of coverage for each gene associated with resistance and 
      regenerates it including information on deletions that occured.
      """
      df_lab_report = pd.read_csv(lab_report, skip_blank_lines=True)
      df_coverage = pd.DataFrame(columns=["Gene", "Percent_Coverage","Warning"])
      
      for gene, percent_coverage in gene_coverage_dict.items():
        if gene == "Gene":
          continue
        warning = ""
        for mutation_type_nucleotide in df_lab_report["tbprofiler_variant_substitution_nt"][df_lab_report["tbprofiler_gene_name"] == gene]:
          if "del" in mutation_type_nucleotide:
            warning = "Deletion identified"
            if float(percent_coverage) == 100:
              warning = "Deletion identified (upstream)"
        df_coverage = df_coverage.append({"Gene": gene, "Percent_Coverage": percent_coverage, "Warning": warning}, ignore_index=True)
      
      df_coverage.to_csv("tbprofiler_coverage_report.csv", index=False)


    ### Report Generation ###

    # get timestamp in YYYY-MM-DD HH:MM format
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

    # Laboratorian report generation
    parse_json_lab_report("~{json}")

    # LIMS report generation
    parse_json_lims_report("~{json}", current_time)

    # LOOKER report generation
    parse_json_looker_report("~{json}", current_time)

    # Coverage report generation
    regenerate_coverage_report(gene_coverage_dict, "tbprofiler_laboratorian_report.csv")

    CODE
  >>>
  output {
    File tbprofiler_looker_csv = "tbprofiler_looker.csv"
    File tbprofiler_laboratorian_report_csv = "tbprofiler_laboratorian_report.csv"
    File tbprofiler_lims_report_csv = "tbprofiler_lims_report.csv"
    File tbprofiler_coverage_report = "tbprofiler_coverage_report.csv"
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