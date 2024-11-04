version 1.0

workflow convert_taxon_ids {
  input {
    Int taxon_id
  }
  String unsupported_organism = "unsupported"
  if (taxon_id == "2697049") {
    String sars_cov_2 = "sars-cov-2"
  }
  if (taxon_id == "10244") {
    String mpox = "MPXV"
  }
  if (taxon_id == "11082") {
    String wnv = "WNV"
  }
  if (taxon_id == "11320") {
    String flu_a = "flu" # flu A
  }
  if (taxon_id == "11520") {
    String flu_b = "flu" # flu B
  }
  if (taxon_id == "12814") {
    String rsv_a = "rsv_a"
  }
  if (taxon_id == "12815") {
    String rsv_b = "rsv_b"
  }
  if (taxon_id == "11676") {
    String hiv = "HIV"
  }
  output {
    organism = select_first([sars_cov_2, mpox, wnv, flu_a, flu_b, rsv_a, rsv_b, hiv, unsupported_organism])
  }
}