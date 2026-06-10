version 1.0

workflow allele_caller_parameters {
  meta {
    description: "Organizes all PulseNet 2.0 Allele Calling parameters into a single location for easier maintainence; similarity thresholds were extracted from <https://github.com/ncezid-biome/pulsenet2.0-bfx/blob/47644186f2df27e9f01a000d47c451135a75f65d/main/utils/utils.nf#L42>"
  }
  input {
    String merlin_tag
  }
  if (merlin_tag == "Campylobacter") {
    Float campy_similarity = 70
    File campy_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/CAMPY.tar.gz"
  }
  if (merlin_tag == "Clostridium botulinum") {
    Float cbot_similarity = 85
    File cbot_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/CBOT.tar.gz"
  }
  if (merlin_tag == "Cronobacter") {
    Float crono_similarity = 80
    File crono_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/CRONO.tar.gz"
  }
  if (merlin_tag == "Escherichia" || merlin_tag == "Shigella sonnei") {
    Float stec_similarity = 85
    File stec_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/STEC.tar.gz"
  }
  if (merlin_tag == "Listeria") {
    Float listeria_similarity = 85
    File listeria_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/LISTERIA.tar.gz"
  }
  if (merlin_tag == "Salmonella") {
    Float salm_similarity = 75
    File salm_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/SALM.tar.gz"
  }
  if (merlin_tag == "Vibrio" || merlin_tag == "Vibrio cholerae") {
    Float vibrio_similarity = 85
    File vibrio_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/VIBR.tar.gz"
  }
  if (merlin_tag == "Yersinia") {
    Float yersinia_similarity = 85
    File yersinia_db = "gs://theiagen-public-resources-rp/reference_data/bacterial/pn2.0-mlst-databases/YERSINIA.tar.gz"
  }
  output {
    Float similarity = select_first([campy_similarity, cbot_similarity, crono_similarity, stec_similarity, listeria_similarity, salm_similarity, vibrio_similarity, yersinia_similarity, 0])
    File db = select_first([campy_db, cbot_db, crono_db, stec_db, listeria_db, salm_db, vibrio_db, yersinia_db, ""])
  }
}
