version 1.0

workflow allele_calling_parameters {
  meta {
    description: "Organizes all PulseNet 2.0 Allele Calling parameters into a single location for easier maintainence; similarity thresholds were extracted from <https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/main/utils/utils.nf#L42>"
  }
  input {
    String merlin_tag
    String gambit_predicted_taxon
  }
  if (merlin_tag == "Campylobacter") {
    File campy_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/CAMPY.tar.gz"
    String campy_scheme = "CAMPY"
    String campy_loci_path = "CAMPY/loci.tsv"
    Float campy_similarity = 70
    String campy_qc_genus = "CAMPY"
    # qc species
    if (gambit_predicted_taxon == "Campylobacter coli") {
      String c_coli_qc_species = "--organism.species COLI"
    }
    if (gambit_predicted_taxon == "Campylobacter fetus") {
      String c_fetus_qc_species = "--organism.species FETUS"
    }
    if (gambit_predicted_taxon == "Campylobacter jejuni") {
      String c_jejuni_qc_species = "--organism.species JEJUNI"
    }
    if (gambit_predicted_taxon == "Campylobacter lari") {
      String c_lari_qc_species = "--organism.species LARI"
    }
    if (gambit_predicted_taxon == "Campylobacter upsaliensis") {
      String c_upsaliensis_qc_species = "--organism.species UPSALIENSIS"
    }
  }
  if (merlin_tag == "Clostridium botulinum") {
    File cbot_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/CBOT.tar.gz"
    String cbot_scheme = "CBOT"
    String cbot_loci_path = "CBOT/loci.tsv"
    Float cbot_similarity = 85
    String cbot_qc_genus = "CBOT"
  }
  if (merlin_tag == "Cronobacter") {
    File crono_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/CRONO.tar.gz"
    String crono_scheme = "CRONO"
    String crono_loci_path = "CRONO/loci.tsv"
    Float crono_similarity = 80
    String crono_qc_genus = "CRONO"
  }
  if (merlin_tag == "Escherichia" || merlin_tag == "Shigella sonnei") {
    File stec_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/STEC.tar.gz"
    String stec_scheme = "STEC"
    String stec_loci_path = "STEC/loci/chromosomal_plasmid.tsv"
    Float stec_similarity = 85
    String stec_qc_genus = "STEC"
  }
  if (merlin_tag == "Listeria") {
    File listeria_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/LISTERIA.tar.gz"
    String listeria_scheme = "LISTERIA"
    String listeria_loci_path = "LISTERIA/loci.tsv"
    Float listeria_similarity = 85
    String listeria_qc_genus = "LISTERIA_" # not sure if we should be using the I, II, III, or IV ones
  }
  if (merlin_tag == "Salmonella") {
    String salm_scheme = "SALM"
    File salm_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/SALM.tar.gz"
    String salm_loci_path = "SALM/loci.tsv"
    Float salm_similarity = 75
    String salm_qc_genus = "SALM"
  }
  if (merlin_tag == "Vibrio" || merlin_tag == "Vibrio cholerae") { # do not run vibrio vulnificus ??
    String vibrio_scheme = "VIBR"
    File vibrio_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/VIBR.tar.gz"
    String vibrio_loci_path = "VIBR/loci/VIBRIO_loci.tsv"
    Float vibrio_similarity = 85
    String vibrio_qc_genus = "VIBRIO"
    # species-specific loci paths
    if (gambit_predicted_taxon == "Vibrio cholerae") {
      String v_cholerae_loci_path = "VIBR/loci/VIBRIO_cholerae_loci.tsv"
    }
    if (gambit_predicted_taxon == "Vibrio parahaemolyticus") {
      String v_parahaemolyticus_loci_path = "VIBR/loci/VIBRIO_parahaemloyticus_loci.tsv"
    }
    # qc species
    if (gambit_predicted_taxon == "Vibrio vulnificus") {
      String v_vulnificus_qc_species = "--organism.species VULNIFICUS"
    }
  }
  if (merlin_tag == "Yersinia") {
    String yersinia_scheme = "YERSINIA"
    File yersinia_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/YERSINIA.tar.gz"
    String yersinia_loci_path = "YERSINIA/loci.csv"
    Float yersinia_similarity = 85
    String yersinia_qc_genus = "YERSINIA"
  }
  output {
    String scheme = select_first([campy_scheme, cbot_scheme, crono_scheme, stec_scheme, listeria_scheme, salm_scheme, vibrio_scheme, yersinia_scheme, ""])
    File db = select_first([campy_db, cbot_db, crono_db, stec_db, listeria_db, salm_db, vibrio_db, yersinia_db, ""])
    String loci_path = select_first([campy_loci_path, cbot_loci_path, crono_loci_path, stec_loci_path, listeria_loci_path, salm_loci_path, v_cholerae_loci_path, v_parahaemolyticus_loci_path, vibrio_loci_path, yersinia_loci_path, ""])
    Float similarity = select_first([campy_similarity, cbot_similarity, crono_similarity, stec_similarity, listeria_similarity, salm_similarity, vibrio_similarity, yersinia_similarity, 0])
    String qc_genus = select_first([campy_qc_genus, cbot_qc_genus, crono_qc_genus, stec_qc_genus, listeria_qc_genus, salm_qc_genus, vibrio_qc_genus, yersinia_qc_genus, ""])
    String qc_species = select_first([c_jejuni_qc_species, c_coli_qc_species, c_fetus_qc_species, c_lari_qc_species, c_upsaliensis_qc_species, v_vulnificus_qc_species, ""])
  }
}
