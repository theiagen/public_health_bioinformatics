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

    ## Lookup Data Structures - GLOBAL VARIABLES ##

    # gene coverage as a dictionary
    GENE_COVERAGE_DICT = pd.read_csv("~{gene_coverage}", delimiter="\t", skip_blank_lines=True).to_dict()
    GENE_COVERAGE_DICT = GENE_COVERAGE_DICT["#NOTE: THE VALUES BELOW ASSUME TBPROFILER (H37Rv) REFERENCE GENOME"] # skip first line

    # lookup dictionary - antimicrobial code to drug name
    ANTIMICROBIAL_CODE_TO_DRUG_NAME = {"M_DST_B01_INH": "isoniazid", 
                          "M_DST_C01_ETO": "ethionamide",
                          "M_DST_D01_RIF": "rifampicin", 
                          "M_DST_E01_PZA": "pyrazinamide",
                          "M_DST_F01_EMB": "ethambutol",
                          "M_DST_G01_AMK": "amikacin", 
                          "M_DST_H01_KAN": "kanamycin",
                          "M_DST_I01_CAP": "capreomycin", 
                          "M_DST_J01_MFX": "moxifloxacin",
                          "M_DST_K01_LFX": "levofloxacin", 
                          "M_DST_L01_BDQ": "bedaquiline",
                          "M_DST_M01_CFZ": "clofazimine", 
                          "M_DST_N01_LZD": "linezolid" 
                         }
    
    # Lookup list - antimicrobial drug names
    ANTIMICROBIAL_DRUG_NAME_LIST = ["isoniazid", "ethionamide", "rifampicin", "pyrazinamide", "ethambutol",
                          "streptomycin", "amikacin", "kanamycin", "capreomycin", "moxifloxacin",
                          "levofloxacin", "bedaquiline", "clofazimine", "linezolid"
                         ]

    # lookup dictionary - antimicrobial code to gene names to antimicrobial code column names
    ANTIMICROBIAL_CODE_TO_GENES = {
                  "M_DST_B01_INH": {"katG": "M_DST_B02_katG", 
                                    "fabG1": "M_DST_B03_fabG1",
                                    "inhA": "M_DST_B04_inhA"},
                  "M_DST_C01_ETO": {"ethA": "M_DST_C02_ethA", 
                                    "fabG1": "M_DST_C03_fabG1",
                                    "inhA": "M_DST_C04_inhA"},
                  "M_DST_D01_RIF": {"rpoB": "M_DST_D02_rpoB"},
                  "M_DST_E01_PZA": {"pncA": "M_DST_E02_pncA"},
                  "M_DST_F01_EMB": {"embA": "M_DST_F02_embA", 
                                    "embB": "M_DST_F03_embB"},
                  "M_DST_G01_AMK": {"rrs": "M_DST_G02_rrs", 
                                    "eis": "M_DST_G03_eis"},
                  "M_DST_H01_KAN": {"rrs": "M_DST_H02_rrs", 
                                    "eis": "M_DST_H03_eis"},
                  "M_DST_I01_CAP": {"rrs": "M_DST_I02_rrs", 
                                    "tlyA": "M_DST_I03_tlyA"},
                  "M_DST_J01_MFX": {"gyrA": "M_DST_J02_gyrA", 
                                    "gyrB": "M_DST_J03_gyrB"},
                  "M_DST_K01_LFX": {"gyrA": "M_DST_K02_gyrA", 
                                    "gyrB": "M_DST_K03_gyrB"},
                  "M_DST_L01_BDQ": {"Rv0678": "M_DST_L02_Rv0678", 
                                    "atpE": "M_DST_L03_atpE",
                                    "pepQ": "M_DST_L04_pepQ", 
                                    "mmpL5": "M_DST_L05_mmpL5",
                                    "mmpS5": "M_DST_L06_mmpS5"},
                  "M_DST_M01_CFZ": {"Rv0678":"M_DST_M02_Rv0678", 
                                    "pepQ": "M_DST_M03_pepQ",
                                    "mmpL5":"M_DST_M04_mmpL5", 
                                    "mmpS5": "M_DST_M05_mmpS5"},
                  "M_DST_N01_LZD": {"rrl": "M_DST_N02_rrl", 
                                    "rplC": "M_DST_N03_rplC"}
                }

    # lookup list - the genes to be considered for the LIMS report
    GENES_FOR_LIMS = ["katG", "fabG1", "inhA", "ethA", "rpoB", "pncA", "embA", "embB", "rrs", "eis", 
                    "tlyA", "gyrA", "gyrB", "Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5", "rrl", "rplC"
                    ]

    # lookup dictionary - gene to antimicrobial drug name, including genes in watchlist
    # (https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv)
    # (https://github.com/jodyphelan/tbdb/blob/master/tbdb.watchlist.csv)
    GENE_TO_ANTIMICROBIAL_DRUG_NAME = {"ahpC":["isoniazid"],
                                        "ald":["cycloserine"],
                                        "alr": ["cycloserine"],
                                        "atpE": ["bedaquiline"],
                                        "ddn": ["delamanid"],
                                        "ddn": ["delamanid"],
                                        "eis": ["amikacin", "kanamycin"],
                                        "embA": ["ethambutol"],
                                        "embB": ["ethambutol"],
                                        "embC": ["ethambutol"],
                                        "embR": ["ethambutol"],
                                        "ethA": ["ethionamide"],
                                        "ethR": ["ethionamide"],
                                        "fabG1": ["ethionamide", "isoniazid"],
                                        "fbiA": ["delamanid"],
                                        "fbiB": ["delamanid"],
                                        "fbiC": ["delamanid"],
                                        "fbiD": ["delamanid"],
                                        "fgd1": ["delamanid"],
                                        "fgd1": ["delamanid"],
                                        "folC": ["para-aminosalicylic_acid"],
                                        "gid": ["streptomycin"],
                                        "gyrA": ["ciprofloxacin", "fluoroquinolones", "levofloxacin", "moxifloxacin", "ofloxacin"],
                                        "gyrB": ["moxifloxacin","ciprofloxacin","fluoroquinolones", "levofloxacin", "ofloxacin"],
                                        "inhA": ["ethionamide", "isoniazid"],
                                        "kasA": ["isoniazid"],
                                        "katG": ["isoniazid"],
                                        "mshA": ["ethionamide"],
                                        "panD": ["pyrazinamide"],
                                        "pepQ": ["bedaquiline", "clofazimine"],
                                        "pncA": ["pyrazinamide"],
                                        "ribD": ["para-aminosalicylic_acid"],
                                        "rplC": ["linezolid"],
                                        "rpoA": ["rifampicin"],
                                        "rpoB": ["rifampicin"],
                                        "rpoC": ["rifampicin"],
                                        "rpsA": ["pyrazinamide"],
                                        "rpsL": ["streptomycin"],
                                        "rrl": ["linezolid"],
                                        "rrs": ["streptomycin", "amikacin", "aminoglycosides", "capreomycin", "kanamycin"],
                                        "Rv0678": ["bedaquiline", "clofazimine"],
                                        "thyA": ["para-aminosalicylic_acid"],
                                        "thyX": ["para-aminosalicylic_acid"],
                                        "tlyA": ["capreomycin"],
                                        "ubiA": ["ethambutol"]
                         }
    
    # lookup dictionary - antimicrobial drug name to gene name  (https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv)
    ANTIMICROBIAL_DRUG_NAME_TO_GENE_NAME = { "amikacin": ["eis", "rrs"], 
                                        "bedaquiline": ["Rv0678"], 
                                        "capreomycin": ["rrs", "tlyA"], 
                                        "clofazimine": ["Rv0678"], 
                                        "ethambutol": ["embA", "embB", "embC", "embR"], 
                                        "ethionamide": ["ethA", "ethR", "fabG1", "inhA"], 
                                        "isoniazid": ["ahpC", "fabG1", "inhA", "kasA", "katG"], 
                                        "kanamycin": ["eis", "rrs"], 
                                        "levofloxacin": ["gyrA", "gyrB"], 
                                        "linezolid": ["rplC", "rrl"],
                                        "moxifloxacin": ["gyrA", "gyrB"],
                                        "pyrazinamide": ["panD", "pncA", "rpsA"], 
                                        "rifampicin": ["rpoA", "rpoB", "rpoC"], 
                                        "streptomycin": ["gid", "rpsL", "rrs"] 
                                      }
    
    # lookup dictionary - gene to tier
    GENE_TO_TIER = {"ahpC": "Tier 1", "inhA": "Tier 1", "katG": "Tier 1", "rpoB": "Tier 1", "embA": "Tier 1", 
                    "embB": "Tier 1", "embC": "Tier 1", "pncA": "Tier 1", "clpC1": "Tier 1", "panD": "Tier 1", 
                    "gyrA": "Tier 1", "gyrB": "Tier 1", "pepQ": "Tier 1", "Rv0678": "Tier 1", "mmpL5": "Tier 1", 
                    "mmpS5": "Tier 1", "atpE": "Tier 1", "rplC": "Tier 1", "rrl": "Tier 1", "fgd1": "Tier 1", 
                    "ddn": "Tier 1", "fbiA": "Tier 1", "fbiB": "Tier 1", "fbiC": "Tier 1", "Rv2983": "Tier 1", 
                    "rrs": "Tier 1", "eis": "Tier 1", "whiB7": "Tier 1", "rpsL": "Tier 1", "gid": "Tier 1", 
                    "Rv1258c": "Tier 1", "ethA": "Tier 1", "tlyA": "Tier 1", "mshA": "Tier 2", "ndh": "Tier 2", 
                    "Rv2752c": "Tier 2", "rpoA": "Tier 2", "rpoC": "Tier 2", "embR": "Tier 2", "ubiA": "Tier 2", 
                    "PPE35": "Tier 2", "Rv3236c": "Tier 2", "Rv1979c": "Tier 2", "whiB6": "Tier 2", "ccsA": "Tier 2", 
                    "fprA": "Tier 2", "aftB": "Tier 2", "ethR": "Tier 2", "Rv3083": "Tier 2", "fabG1": "Tier 1"
                  }
                  
    # lookup dictionary - gene to locus tag (https://github.com/jodyphelan/TBProfiler/blob/master/db/tbdb.bed)
    GENE_TO_LOCUS_TAG = {"ahpC":"Rv2428", "ald":"Rv2780", "alr": "Rv3423c",
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
                          "thyX": "Rv2754c", "tlyA": "Rv1694", "rpoA": "Rv3457c",
                          "fgd1": "Rv0407", "fbiB": "Rv3262", "fbiC": "Rv1173",
                          "ddn": "Rv3547", "ubiA": "Rv3806c", "atpE": "Rv1305",
                          "mshA": "Rv0486", "pepQ": "Rv2535c", "fbiD": "Rv2983"
                         }

    GENE_LIST_OPTION_1 = ["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5", "rrl", "rplC"] # Rv0678 is mmpR
    GENE_LIST_OPTION_2 = ["katG", "pncA", "ethA", "gid", "rpoB"]
    
    # Turning TBProfiler annotations into Looker or MDL interpretations
    ANNOTATION_TO_INTERPRETATION = {
      "Assoc w R": {"looker": "R", 
                    "MDL": "R"},
      "Assoc w R - interim": {"looker": "R-Interim", 
                              "MDL": "U"},
      "Uncertain significance": {"looker": "U", 
                                 "MDL": "S", 
                                 "MDL-ingenelist1" : "U"},
      "Not assoc w R": {"looker": "S", 
                        "MDL": "S"},
      "Not assoc w R - Interim": {"looker": "S-Interim", 
                                  "MDL": "S", 
                                  "MDL-ingenelist1" : "U"}                              
    }

    # genes with promoter regions to consider
    PROMOTER_REGIONS = {"Rv0678": [-84, -1],
                        "atpE": [-48, -1],
                        "pepQ": [-33, -1],
                        "rplC": [-18, -1]
                        }

    # genes that have positions with special consideration
    SPECIAL_POSITIONS = {"rrl": [[2003, 2367], [2449, 3056]],
                         "rpoB": [426, 452],
                         "rss": [1401, 1402, 1484]
                        }

    LOW_DEPTH_OF_COVERAGE_LIST = []
    
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

    def apply_expert_rules(nucleotide_change, protein_change, gene, substitution_type, interpretation_destination):
      """
      Apply rules 1-3
      """

      position_nt = get_position(nucleotide_change)
      position_aa = get_position(protein_change)

      if gene in ["Rv0678", "atpE", "pepQ", "rplC", "mmpL5", "mmpS5"]: # apply expert rules 1.2         
        # check if position within promoter regions
        if PROMOTER_REGIONS[gene][1] <= position_nt <= PROMOTER_REGIONS[gene][2]: 
          return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
        elif "upstream_gene_variant" in substitution_type: # otherwise, check if it is an upstream gene variant
          return "S" if interpretation_destination == "MDL" else "U"
        else: # otherwise, apply expert rules 1.2
          if not any(non_ORF in nucleotide_change for non_ORF in ["+", "-", "*"]) or nucleotide_change.endswith("*"): 
          # if a position includes either +, *, or - it's not in the ORF 
          # UNLESS the * is at the end which means its a premature stop codon
            if substitution_type != "synonymous_variant":
              return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
            else:
              return "S"

      elif gene == "rrl": # apply expert rules 1.2
        if (SPECIAL_POSITIONS[gene][1][1] <= position_nt <= SPECIAL_POSITIONS[gene][1][2]) or (SPECIAL_POSITIONS[gene][2][1] <= position_nt <= SPECIAL_POSITIONS[gene][2][2]):
          return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
        else:
          return "S" if interpretation_destination == "MDL" else "U"

      elif gene in ["katG", "pncA", "ethA", "gid"]: # apply expert rules 2.2.1
        if any(indel_or_stop in nucleotide_change for indel_or_stop in ["del", "ins", "fs", "delins", "_"]) or nucleotide_change.endswith("*"):
          return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
        else:
            if substitution_type != "synonymous_variant" or "upstream_gene_variant" in substitution_type:
              return "S" if interpretation_destination == "MDL" else "U"
            else:
              return "S"

      elif gene == "rpoB": # apply expert rules 2.2.2
        if SPECIAL_POSITIONS[gene][1] <= position_nt <= SPECIAL_POSITIONS[gene][2]:
            if substitution_type != "synonymous_variant":
              return "Assoc with R" if interpretation_destination == "LIMS" else "R"
            else:
              return "S"   
        else:
            if substitution_type != "synonymous_variant" or "upstream_gene_variant" in substitution_type:
              return "S" if interpretation_destination == "MDL" else "U"
            else:
              return "S"

      elif gene not in GENE_LIST_OPTION_1 or gene not in GENE_LIST_OPTION_2: # NOT AN EXPERT RULE: 3.2
        if gene == "rrs": # apply rule 3.2.1
          if position_nt in SPECIAL_POSITIONS[gene]:
            return "Unoexpert"
          else:
            return "Snoexpert" if interpretation_destination == "MDL" else "Unoexpert"
        elif substitution_type != "synonymous_variant" or "upstream_gene_variant" in substitution_type:
          return "Snoexpert" if interpretation_destination == "MDL" else "Unoexpert"
        else:
          return "Snoexpert"

      return ""

    def remove_no_expert(row):
      """
      This function removes the 'noexpert' suffix in the case where 
      the logic applied is not considered an expert rule.
      """
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
          # if 90% of genes have coverage above threshold:
            # return "DNA of Mycobacterium tuberculosis complex detected"
          # else:
            # return "DNA of Mycobacterium tuberculosis complex NOT detected"
        elif "La1" in results_json["main_lin"]:
          return "DNA of M. tuberculosis complex detected (M. bovis)"
        elif "lineage" in results_json["main_lin"]:
          return "DNA of Mycobacterium tuberculosis species detected"
        elif "BCG" in results_json["main_lin"]:
          return "DNA of Mycobacterium bovis BCG detected"
        elif "bovis" in results_json["main_lin"]:
          return "DNA of non-BCG Mycobacterium bovis detected"
        else:
          return "DNA of Mycobacterium tuberculosis complex detected (not M. bovis and not M. tb)"
    
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
        
        for dr_variant in results_json["dr_variants"]:
          name = dr_variant["gene"]
          if name in GENES_FOR_LIMS:
            substitution = "{} ({})".format(dr_variant["nucleotide_change"], dr_variant["protein_change"])

            if name not in mutations_dict.keys():
              mutations_dict[name] = substitution
            else:
              mutations_dict[name] = "{}; {}".format(mutations_dict[name], substitution)
        
        for other_variant in results_json["other_variants"]:
          name = other_variant["gene"]
          if name in GENES_FOR_LIMS:
            aa_position = get_position(other_variant["protein_change"])
            
            # report all non-synonymous mutations, unless rpoB RRDR
            if other_variant["type"] != "synonymous_variant" or (other_variant["gene"] == "rpoB" and (SPECIAL_POSITIONS["rpoB"][1] <= position_nt <= SPECIAL_POSITIONS["rpoB"][2])): 
              name = other_variant["gene"]
              substitution = "{} ({})".format(other_variant["nucleotide_change"], other_variant["protein_change"])

              if other_variant["type"] == "synonymous_variant": # only catches rpoB RRDR
                substitution = "{} [synonymous]".format(substitution)
                
              if name not in mutations_dict.keys():
                mutations_dict[name] = substitution
              else:
                mutations_dict[name] = "{}; {}".format(mutations_dict[name], substitution)
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
        return "The detected genetic determinant(s) have uncertain significance. Resistance to {} cannot be ruled out".format(drug)
      else: # "Not assoc w R" and "Not assoc w R - Interim" and anything else
        return "No mutations associated with resistance to {} detected".format(drug)

    def parse_json_resistance(json_file, destination):
      """
      This function parses the tbprofiler json report and returns a resistance dictionary
      containing the WHO resistance annotation for each antimicrobial drug for LIMS and Looker. 
      The annotation corresponds to the highest ranked one regarding severity (R > I > S)
      """
      resistance_dict = {}


      with open(json_file) as js_fh:
        results_json = json.load(js_fh)

        for dr_variant in results_json["dr_variants"]: # mutation 
          name = dr_variant["gene"]

          if dr_variant["type"] != "synonymous_variant":  # report all non-synonymous mutations

            if "annotation" in dr_variant.keys(): # if an annotation is present,
              if len(dr_variant["annotation"]) == 0:
                for identified_drug in dr_variant["gene_associated_drugs"]:
                  who_annotation = apply_expert_rules(dr_variant["nucleotide_change"], dr_variant["protein_change"], dr_variant["gene"], dr_variant["type"], "LIMS")

                  if identified_drug not in resistance_dict.keys():
                    if (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                      resistance_dict[identified_drug] = who_annotation # overwrite with more severe annotation
                  else:
                    # if current annotation indicates higher severity than any previous annotation,
                    if rank_annotation(resistance_dict[identified_drug]) < rank_annotation(who_annotation): 
                      if (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                        resistance_dict[identified_drug] = who_annotation # overwrite with more severe annotation

              for annotation in dr_variant["annotation"]: # iterate through them
                drug = annotation["drug"]
                who_annotation = annotation["who_confidence"]
                if drug not in resistance_dict.keys():
                  if destination == "Looker":
                    resistance_dict[drug] = who_annotation
                  elif (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                    resistance_dict[drug] = who_annotation

                else: # if the drug has already been seen in either the same variant or a different one,
                  # if current annotation indicates higher severity than any previous annotation,
                  if rank_annotation(resistance_dict[drug]) < rank_annotation(who_annotation): 
                    if destination == "Looker": 
                      resistance_dict[drug] = who_annotation # overwrite with more severe annotation
                    elif (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                      resistance_dict[drug] = who_annotation # overwrite with more severe annotation
            else:
              for identified_drug in dr_variant["gene_associated_drugs"]:
                who_annotation = apply_expert_rules(dr_variant["nucleotide_change"], dr_variant["protein_change"], dr_variant["gene"], dr_variant["type"], "LIMS")
        
                if identified_drug not in resistance_dict.keys():
                  if (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                    resistance_dict[identified_drug] = who_annotation # overwrite with more severe annotation
                else:
                  # if current annotation indicates higher severity than any previous annotation,
                  if rank_annotation(resistance_dict[identified_drug]) < rank_annotation(who_annotation): 
                    if (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                      resistance_dict[identified_drug] = who_annotation # overwrite with more severe annotation

        for other_variant in results_json["other_variants"]:
          name = other_variant["gene"]
          if other_variant["type"] != "synonymous_variant":  # report all non-synonymous mutations

            if "annotation" in other_variant.keys(): # if an annotation is present,
              if len(other_variant["annotation"]) == 0:
                for identified_drug in other_variant["gene_associated_drugs"]:
                  who_annotation = apply_expert_rules(other_variant["nucleotide_change"], other_variant["protein_change"], other_variant["gene"], other_variant["type"], "LIMS")

                  if identified_drug not in resistance_dict.keys():
                    if (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                      resistance_dict[identified_drug] = who_annotation # overwrite with more severe annotation
                  else:
                    # if current annotation indicates higher severity than any previous annotation,
                    if rank_annotation(resistance_dict[identified_drug]) < rank_annotation(who_annotation): 
                      if (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                        resistance_dict[identified_drug] = who_annotation # overwrite with more severe annotation
              
              for annotation in other_variant["annotation"]: # iterate through them
                drug = annotation["drug"]
                who_annotation = annotation["who_confidence"]
                
                if drug not in resistance_dict.keys():
                  if destination == "Looker":
                    resistance_dict[drug] = who_annotation
                  elif (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                    resistance_dict[drug] = who_annotation # only 

                else: # if the drug has already been seen in either the same variant or a different one,
                  # if current annotation indicates higher severity than any previous annotation,
                  if rank_annotation(resistance_dict[drug]) < rank_annotation(who_annotation): 
                    if destination == "Looker": 
                      resistance_dict[drug] = who_annotation # overwrite with more severe annotation
                    elif (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                      resistance_dict[drug] = who_annotation # overwrite with more severe annotation
            else:
              for identified_drug in other_variant["gene_associated_drugs"]:
                who_annotation = apply_expert_rules(other_variant["nucleotide_change"], other_variant["protein_change"], other_variant["gene"], other_variant["type"], "LIMS")

                if identified_drug not in resistance_dict.keys():
                  if (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                    resistance_dict[identified_drug] = who_annotation # overwrite with more severe annotation
                else:
                  # if current annotation indicates higher severity than any previous annotation,
                  if rank_annotation(resistance_dict[identified_drug]) < rank_annotation(who_annotation): 
                    if (name in GENES_FOR_LIMS) and (destination == "LIMS"):
                      resistance_dict[identified_drug] = who_annotation # overwrite with more severe annotation

      return resistance_dict

    def variant_to_row(variant):
      """
      This function recieved a variants dictionary and returns a row dictionary containing the
      basic information for it to be filled in the laboratorian report
      """
      row = {}
      row["sample_id"] = "~{samplename}"
      row["tbprofiler_gene_name"] = variant["gene"]
      if variant["gene"] in GENE_TO_TIER.keys():
        row["gene_tier"] = GENE_TO_TIER[variant["gene"]]
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
      if row["depth"] < int('~{min_depth}') or float(GENE_COVERAGE_DICT[variant["gene"]]) < ~{coverage_threshold}:
        row["warning"] = "Insufficient coverage in locus"
        if "del" in variant["nucleotide_change"]:
          row["warning"] = "Insufficient coverage in locus (deletion identified)"
        else:
          LOW_DEPTH_OF_COVERAGE_LIST.append(variant["gene"])

      return row

    ## Main Parsing Functions ## 

    def parse_json_lab_report(json_file):
      """
      This function receives the tbprofiler output json file and
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
                  drugs_to_row[identified_drug] = {"other_variant": dr_variant, 
                                                   "who_confidence": row["confidence"], 
                                                   "drug": identified_drug, 
                                                   "nucleotide_change": dr_variant["nucleotide_change"], 
                                                   "protein_change": dr_variant["protein_change"], 
                                                   "gene": dr_variant["gene"], 
                                                   "type": dr_variant["type"]}
                # overwrite entry with the more severe annotation (higher value) if multiple drugs are present
                elif rank_annotation(drugs_to_row[identified_drug]["who_confidence"]) < rank_annotation(row["confidence"]): 
                  drugs_to_row[identified_drug] = {"other_variant": dr_variant, 
                                                   "who_confidence": row["confidence"], 
                                                   "drug": identified_drug, 
                                                   "nucleotide_change": dr_variant["nucleotide_change"], 
                                                   "protein_change": dr_variant["protein_change"], 
                                                   "gene": dr_variant["gene"], 
                                                   "type": dr_variant["type"]}

            # case: drug confers resistance to multiple drugs - if the same drug shows multiple times in a single mutation, save only the most severe annotation
            # iterate thorugh all possible annotations for the variant
            for annotation in dr_variant["annotation"]: 
              # if this is the first time a drug is seen, add to dictionary
              if annotation["drug"] not in drugs_to_row: 
                drugs_to_row[annotation["drug"]] = {"other_variant": dr_variant, 
                                                    "who_confidence": annotation["who_confidence"], 
                                                    "drug": annotation["drug"], 
                                                    "nucleotide_change": dr_variant["nucleotide_change"], 
                                                    "protein_change": dr_variant["protein_change"], 
                                                    "gene": dr_variant["gene"], 
                                                    "type": dr_variant["type"]}
              # overwrite entry with the more severe annotation (higher value) if multiple drugs are present
              elif rank_annotation(drugs_to_row[annotation["drug"]]["who_confidence"]) < rank_annotation(annotation["who_confidence"]): 
                drugs_to_row[annotation["drug"]] = {"other_variant": dr_variant, 
                                                    "who_confidence": annotation["who_confidence"], 
                                                    "drug": annotation["drug"], 
                                                    "nucleotide_change": dr_variant["nucleotide_change"], 
                                                    "protein_change": dr_variant["protein_change"], 
                                                    "gene": dr_variant["gene"], 
                                                    "type": dr_variant["type"]}
            
            for drug in drugs_to_row:
              gene = drugs_to_row[drug]["gene"]

              row = variant_to_row(drugs_to_row[drug]["other_variant"])

              row["confidence"] = "No WHO annotation" if drugs_to_row[drug]["who_confidence"] == "" else drugs_to_row[drug]["who_confidence"]
              
              row["antimicrobial"] = drugs_to_row[drug]["drug"]

              # annotation to interpretation logic
              if row["confidence"] != "No WHO annotation":
                row["looker_interpretation"] = ANNOTATION_TO_INTERPRETATION[row["confidence"]]["looker"]
                row["mdl_interpretation"] = ANNOTATION_TO_INTERPRETATION[row["confidence"]]["MDL" if gene not in GENE_LIST_OPTION_1 else "MDL-ingenelist1"]
                row["rationale"] = "WHO classification"
              else:
                row["looker_interpretation"] = apply_expert_rules(drugs_to_row[drug]["nucleotide_change"], drugs_to_row[drug]["protein_change"], gene, drugs_to_row[drug]["type"], "looker")
                row["mdl_interpretation"] = apply_expert_rules(drugs_to_row[drug]["nucleotide_change"], drugs_to_row[drug]["protein_change"], gene, drugs_to_row[drug]["type"], "MDL")
                row["rationale"] = "Expert rule applied"
              
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
                drugs_to_row[annotation["drug"]] = {"other_variant": other_variant, 
                                                    "who_confidence": annotation["who_confidence"], 
                                                    "drug": annotation["drug"], 
                                                    "nucleotide_change": other_variant["nucleotide_change"],
                                                    "protein_change": other_variant["protein_change"], 
                                                    "gene": other_variant["gene"], 
                                                    "type": other_variant["type"]}
              elif rank_annotation(drugs_to_row[annotation["drug"]]["who_confidence"]) < rank_annotation(annotation["who_confidence"]):
                drugs_to_row[annotation["drug"]] = {"other_variant": other_variant, 
                                                    "who_confidence": annotation["who_confidence"], 
                                                    "drug": annotation["drug"], 
                                                    "nucleotide_change": other_variant["nucleotide_change"],
                                                    "protein_change": other_variant["protein_change"], 
                                                    "gene": other_variant["gene"], 
                                                    "type": other_variant["type"]}

            for drug in drugs_to_row:
              row = variant_to_row(drugs_to_row[drug]["other_variant"])
              
              row["confidence"] = "No WHO annotation" if drugs_to_row[drug]["who_confidence"] == "" else drugs_to_row[drug]["who_confidence"]
              row["antimicrobial"] = drugs_to_row[drug]["drug"]
              
              # annotation to interpretation logic
              if row["confidence"] != "No WHO annotation":
                row["looker_interpretation"] = ANNOTATION_TO_INTERPRETATION[row["confidence"]]["looker"]
                row["mdl_interpretation"] = ANNOTATION_TO_INTERPRETATION[row["confidence"]]["MDL" if gene not in GENE_LIST_OPTION_1 else "MDL-ingenelist1"]
                row["rationale"] = "WHO classification"
              else:
                row["looker_interpretation"] = apply_expert_rules(drugs_to_row[drug]["nucleotide_change"], drugs_to_row[drug]["protein_change"], gene, drugs_to_row[drug]["type"], "looker")
                row["mdl_interpretation"] = apply_expert_rules(drugs_to_row[drug]["nucleotide_change"], drugs_to_row[drug]["protein_change"], gene, drugs_to_row[drug]["type"], "MDL")
                row["rationale"] = "Expert rule applied"
              
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
              
      for gene, antimicrobial_drug_names in GENE_TO_ANTIMICROBIAL_DRUG_NAME.items():
        for drug_name in antimicrobial_drug_names:
          if gene not in genes_reported:

              row = {}
              row["sample_id"] = "~{samplename}"
              row["tbprofiler_gene_name"] = gene
              row["tbprofiler_locus_tag"] = GENE_TO_LOCUS_TAG[gene]
              row["tbprofiler_variant_substitution_nt"] = "NA"
              row["tbprofiler_variant_substitution_aa"] = "NA"
              if float(GENE_COVERAGE_DICT[gene]) >= ~{coverage_threshold}:
                row["tbprofiler_variant_substitution_type"] = "WT"
                row["looker_interpretation"] = "S"
                row["mdl_interpretation"] = "WT"
                row["tbprofiler_variant_substitution_nt"] = "WT"
                row["tbprofiler_variant_substitution_aa"] = "WT"
              else:
                row["tbprofiler_variant_substitution_type"] = "Insufficient Coverage"
                row["looker_interpretation"] = "Insufficient Coverage"
                row["mdl_interpretation"] = "Insufficient Coverage"
              if gene in GENE_TO_TIER.keys():
                row["gene_tier"] = GENE_TO_TIER[gene]
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
      df_laboratorian.to_csv("~{samplename}_tbprofiler_laboratorian_report.csv", index=False)
    
    def parse_json_lims_report(json_file, formatted_time):
      """
      This function recieves the tbprofiler output json file and
      writes the LIMS report that includes the following information
      per sample: 
        - MDL sample accession numbers: includes sample name
        - M_DST_A01_ID - includes lineage
        - The set of information in ANTIMICROBIAL_CODE_TO_GENES dictionary with target drug resistance information in layman's terms, and the mutations responsible for the predicted phenotype
        - Date of analysis in YYYY-MM-DD HH:SS format
        - Operator information
      """
    
      lineage = get_lineage_LIMS("~{json}")
      mutations = parse_json_mutations_for_LIMS("~{json}")
      resistance_annotation = parse_json_resistance("~{json}", "LIMS")
      df_lims = pd.DataFrame({"MDL sample accession numbers":"~{samplename}", "M_DST_A01_ID": lineage}, index=[0])

      for antimicrobial_code, genes in ANTIMICROBIAL_CODE_TO_GENES.items():
        drug_name = ANTIMICROBIAL_CODE_TO_DRUG_NAME[antimicrobial_code]
        # if the drug has been mentioned in the results file
        if drug_name in resistance_annotation.keys(): 
          df_lims[antimicrobial_code] = annotation_to_LIMS(resistance_annotation[drug_name], drug_name)
        else: # the drug is not in the results file
          df_lims[antimicrobial_code] = "No mutations associated with resistance to {} detected".format(drug_name)

        # iterate through the genes that are associated with resistance to each drug
        for gene_name, gene_column_code in genes.items():
          # if the gene is mentioned in the mutations
          if gene_name in mutations.keys():
            # put the formatted mutation as the content for the column
            df_lims[gene_column_code] = mutations[gene_name]
            if gene_name == "rpoB": # rule 5.2.1.2
              if df_lims[antimicrobial_code][0] == "No mutations associated with resistance to {} detected".format(drug_name):
                # if any mutations are present
                if len(mutations) > 0: 
                  non_synomynous_count = 0
                  # check if that mutation is synonymous 
                  # (only output if rpoB RRDR -- see the parse_json_mutations_for_LIMS function)
                  for mutation in mutations: 
                    # if any nonsynymous mutations were identified.
                    if "synonymous" not in mutation: 
                      # keep the original output for the antimicrobial code
                      non_synomynous_count += 1 
                  # otherwise, the only synonymous mutations were identified in rpoB RRDR
                  if non_synomynous_count == 0: 
                    df_lims[antimicrobial_code] = "No mutations associated with resistance to rifampin detected. The detected synonymous mutation(s) do not confer resistance but may result in false-resistance in PCR-based assays targeting the rpoB RRDR."
            
            if gene_name == "rrl": # Rule 5.2.2
              # do not report mutations for rrl
              df_lims[gene_column_code] = "" 
              
            # if the mutations detected were only "S", 
            if df_lims[antimicrobial_code][0] == "No mutations associated with resistance to {} detected".format(drug_name): 
              df_lims[gene_column_code] = "No high confidence mutations detected"
          else: # the gene is not in the mutations list but has decent coverage
            if float(GENE_COVERAGE_DICT[gene_name]) < ~{coverage_threshold}: 
              df_lims[gene_column_code] = "No sequence"
            else:
              df_lims[gene_column_code] = "No mutations detected"

          # HOWEVER, if the coverage is less than the indicated threshold
          if float(GENE_COVERAGE_DICT[gene_name]) < ~{coverage_threshold}: 
            df_lims[gene_column_code] = "Insufficient Coverage"
            # catch for when there's no mutation on gene name but coverage is below threshold
            try: 
              if "del" in mutations[gene_name]:
                df_lims[gene_column_code] = "Insufficient Coverage (deletion identified)"
                # if this case is part of rule 2.2.1, output the mutation
                if gene_name in ["katG", "pncA", "ethA", "gid"]:
                  df_lims[gene_column_code] = mutations[gene_name]
            except:
              df_lims[antimicrobial_code] = "Insufficient Coverage" 
            
            # in addition, if the indicated annotation for the drug is not resistant (less than 4)
            if drug_name in resistance_annotation.keys() and int(rank_annotation(resistance_annotation[drug_name])) < 4:
              # catch for when there's no mutation on gene name but coverage is below threshold
              try: 
                if gene_name not in ["katG", "pncA", "ethA", "gid"] and "del" not in mutations[gene_name]:
                  df_lims[antimicrobial_code] = "Pending Retest"
              except: # there are no mutations in the gene
                df_lims[antimicrobial_code] = "Insufficient Coverage"
              


      df_lims["Analysis date"] = formatted_time
      df_lims["Operator"] = "~{operator}"
      df_lims.to_csv("~{samplename}_tbprofiler_lims_report.csv", index=False)
    
    def parse_json_looker_report(json_file, current_time):
      """
      This function recieves the tbprofiler output json file and
      writes the Looker report that includes the following information
      per sample: 
        - sample_id: includes sample name
        - for each antimicrobial, indication if its resistant (R) or susceptible (S)
      """

      lineage, ID = get_lineage_and_ID_Looker("~{json}")
      resistance_annotation = parse_json_resistance("~{json}", "Looker")
      df_looker = pd.DataFrame({"sample_id":"~{samplename}", "output_seq_method_type": "~{output_seq_method_type}"}, index=[0])

      for antimicrobial_drug in ANTIMICROBIAL_DRUG_NAME_LIST:
        if antimicrobial_drug in resistance_annotation.keys():
          who_annotation = resistance_annotation[antimicrobial_drug]
          df_looker[antimicrobial_drug] = ANNOTATION_TO_INTERPRETATION[who_annotation]["looker"]
          for gene in ANTIMICROBIAL_DRUG_NAME_TO_GENE_NAME[antimicrobial_drug]:
            # indicate warning if any genes failed to achieve 100% coverage_threshold and/or minimum depth  (10x) 
            if df_looker[antimicrobial_drug][0] != "R" and gene in LOW_DEPTH_OF_COVERAGE_LIST:
              df_looker[antimicrobial_drug] = "Insufficient coverage for the locus"
        else: # the antimicrobial drug was not present in the results
          df_looker[antimicrobial_drug] = "S"
      
      df_looker["lineage"] = lineage 
      df_looker["ID"] = ID
      df_looker["analysis_date"] = current_time
      df_looker["operator"] = "~{operator}"
    
      df_looker.to_csv("~{samplename}_tbprofiler_looker.csv", index=False)

    def regenerate_coverage_report(GENE_COVERAGE_DICT, lab_report):
      """
      This function takes the information generated previously regarding
      breadth of coverage for each gene associated with resistance and 
      regenerates it including information on deletions that occured.
      """
      df_lab_report = pd.read_csv(lab_report, skip_blank_lines=True)
      df_coverage = pd.DataFrame(columns=["Gene", "Percent_Coverage","Warning"])
      
      for gene, percent_coverage in GENE_COVERAGE_DICT.items():
        if gene == "Gene":
          continue
        warning = ""
        # sometimes no mutations were identified for a given gene, in those cases, ignore them
        # and save them to stdout
        try: 
          for mutation_type_nucleotide in df_lab_report["tbprofiler_variant_substitution_nt"][df_lab_report["tbprofiler_gene_name"] == gene]:
            if "del" in mutation_type_nucleotide:
              warning = "Deletion identified"
              if float(percent_coverage) == 100:
                warning = "Deletion identified (upstream)"
        except:
          print(gene, df_lab_report["tbprofiler_variant_substitution_nt"][df_lab_report["tbprofiler_gene_name"] == gene])
        df_coverage = df_coverage.append({"Gene": gene, "Percent_Coverage": percent_coverage, "Warning": warning}, ignore_index=True)
      
      df_coverage.to_csv("~{samplename}_tbprofiler_coverage_report.csv", index=False)

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
    regenerate_coverage_report(GENE_COVERAGE_DICT, "~{samplename}_tbprofiler_laboratorian_report.csv")

    CODE
  >>>
  output {
    File tbprofiler_looker_csv = "~{samplename}_tbprofiler_looker.csv"
    File tbprofiler_laboratorian_report_csv = "~{samplename}_tbprofiler_laboratorian_report.csv"
    File tbprofiler_lims_report_csv = "~{samplename}_tbprofiler_lims_report.csv"
    File tbprofiler_coverage_report = "~{samplename}_tbprofiler_coverage_report.csv"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.2"
    memory: "4 GB"
    cpu: 1
    disks: "local-disk " + 10 + " SSD"
    disk: 10 + " GB"
    maxRetries: 0 
  }
}