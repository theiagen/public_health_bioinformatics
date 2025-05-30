version 1.0

import "../../tasks/alignment/task_mafft.wdl" as mafft_task
import "../../tasks/gene_typing/drug_resistance/task_flu_antiviral_subs.wdl" as flu_antiviral

workflow flu_antiviral_substitutions {
  input {
    File? ha_segment_assembly
    File? na_segment_assembly
    File? pa_segment_assembly
    File? pb1_segment_assembly
    File? pb2_segment_assembly
    File? mp_segment_assembly
    String abricate_flu_subtype
    String irma_flu_subtype
     # Amino acid sunstitutions
    File flu_h1_ha_ref = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1_ha.fasta"
    File flu_h3_ha_ref = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3_ha.fasta"
    File flu_n1_na_ref = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_n1_na.fasta"
    File flu_n2_na_ref = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_n2_na.fasta"
    File flu_pa_ref = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_pa.fasta"
    File flu_pb1_ref = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_pb1.fasta"
    File flu_pb2_ref = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_pb2.fasta"
    File flu_h1n1_m2_ref = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h1n1_m2.fasta"
    File flu_h3n2_m2_ref = "gs://theiagen-public-resources-rp/reference_data/viral/flu/reference_h3n2_m2.fasta"
    String? antiviral_aa_subs #user input for antiviral aa subs to be reported
  }
# Identify AA changes for Influenza NA, PA, PB1 and PB2 segments, and associated antiviral mutations
  if (defined(ha_segment_assembly) && ((abricate_flu_subtype == "H1N1") || (irma_flu_subtype == "H1N1"))) {
    call mafft_task.mafft as mafft_h1_ha {
      input:
        genomes = select_all([flu_h1_ha_ref, ha_segment_assembly])
    }
    call flu_antiviral.aa_subs as aa_subs_h1_ha {
      input:
        alignment = mafft_h1_ha.msa,
        protein_name = "HA",
        subtype = "H1N1"
    }
    call flu_antiviral.antiviral_mutations_parser as flu_antiviral_parser_h1_ha {
      input:
        mutations_tsv = aa_subs_h1_ha.aa_changes_tsv,
        antiviral_aa_subs = antiviral_aa_subs
    }
  }
  if (defined(na_segment_assembly) && ((abricate_flu_subtype == "H1N1") || (irma_flu_subtype == "H1N1"))) {
    call mafft_task.mafft as mafft_n1_na {
      input:
        genomes = select_all([flu_n1_na_ref, na_segment_assembly])
    }
    call flu_antiviral.aa_subs as aa_subs_n1_na {
      input:
        alignment = mafft_n1_na.msa,
        protein_name = "NA",
        subtype = "H1N1"
    }
    call flu_antiviral.antiviral_mutations_parser as flu_antiviral_parser_n1_na {
      input:
        mutations_tsv = aa_subs_n1_na.aa_changes_tsv,
        antiviral_aa_subs = antiviral_aa_subs
    }
  }
  if (defined(mp_segment_assembly) && ((abricate_flu_subtype == "H1N1") || (irma_flu_subtype == "H1N1"))) {
    call mafft_task.mafft as mafft_h1n1_mp {
      input:
        genomes = select_all([flu_h1n1_m2_ref, mp_segment_assembly])
    }
    call flu_antiviral.aa_subs as aa_subs_h1n1_mp {
      input:
        alignment = mafft_h1n1_mp.msa,
        protein_name = "MP",
        subtype = "H1N1"
    }
    call flu_antiviral.antiviral_mutations_parser as flu_antiviral_parser_h1n1_mp {
      input:
        mutations_tsv = aa_subs_h1n1_mp.aa_changes_tsv,
        antiviral_aa_subs = antiviral_aa_subs
    }
  }
  if (defined(ha_segment_assembly) && ((abricate_flu_subtype == "H3N2") || (irma_flu_subtype == "H3N2"))) {
    call mafft_task.mafft as mafft_h3_ha {
      input:
        genomes = select_all([flu_h3_ha_ref, ha_segment_assembly])
    }
    call flu_antiviral.aa_subs as aa_subs_h3_ha {
      input:
        alignment = mafft_h3_ha.msa,
        protein_name = "HA",
        subtype = "H3N2"
    }
    call flu_antiviral.antiviral_mutations_parser as flu_antiviral_parser_h3_ha {
      input:
        mutations_tsv = aa_subs_h3_ha.aa_changes_tsv,
        antiviral_aa_subs = antiviral_aa_subs
    }
  }
  if (defined(na_segment_assembly) && ((abricate_flu_subtype == "H3N2") || (irma_flu_subtype == "H3N2"))) {
    call mafft_task.mafft as mafft_n2_na {
      input:
        genomes = select_all([flu_n2_na_ref,na_segment_assembly])
    }
    call flu_antiviral.aa_subs as aa_subs_n2_na {
      input:
        alignment = mafft_n2_na.msa,
        protein_name = "NA",
        subtype = "H3N2"
    }
    call flu_antiviral.antiviral_mutations_parser as flu_antiviral_parser_n2_na {
      input:
        mutations_tsv = aa_subs_n2_na.aa_changes_tsv,
        antiviral_aa_subs = antiviral_aa_subs
    }
  }
  if (defined(mp_segment_assembly) && ((abricate_flu_subtype == "H3N2") || (irma_flu_subtype == "H3N2"))) {
    call mafft_task.mafft as mafft_h3n2_mp {
      input:
        genomes = select_all([flu_h3n2_m2_ref, mp_segment_assembly])
    }
    call flu_antiviral.aa_subs as aa_subs_h3n2_mp {
      input:
        alignment = mafft_h3n2_mp.msa,
        protein_name = "MP",
        subtype = "H3N2"
    }
    call flu_antiviral.antiviral_mutations_parser as flu_antiviral_parser_h3n2_mp {
      input:
        mutations_tsv = aa_subs_h3n2_mp.aa_changes_tsv,
        antiviral_aa_subs = antiviral_aa_subs
    }
  }
  if (defined(pa_segment_assembly)) {
    call mafft_task.mafft as mafft_pa {
      input:
        genomes = select_all([flu_pa_ref, pa_segment_assembly])
    }
    call flu_antiviral.aa_subs as aa_subs_pa {
      input:
        alignment = mafft_pa.msa,
        protein_name = "PA"
    }
    call flu_antiviral.antiviral_mutations_parser as flu_antiviral_parser_pa {
      input:
        mutations_tsv = aa_subs_pa.aa_changes_tsv,
        antiviral_aa_subs = antiviral_aa_subs
    }
  }
  if (defined(pb1_segment_assembly)) {
    call mafft_task.mafft as mafft_pb1 {
      input:
        genomes = select_all([flu_pb1_ref, pb1_segment_assembly])
    }
    call flu_antiviral.aa_subs as aa_subs_pb1 {
      input:
        alignment = mafft_pb1.msa,
        protein_name = "PB1"
    }
    call flu_antiviral.antiviral_mutations_parser as flu_antiviral_parser_pb1 {
      input:
        mutations_tsv = aa_subs_pb1.aa_changes_tsv,
        antiviral_aa_subs = antiviral_aa_subs
    }
  }
  if (defined(pb2_segment_assembly)) {
    call mafft_task.mafft as mafft_pb2 {
      input:
        genomes = select_all([flu_pb2_ref, pb2_segment_assembly])
    }
    call flu_antiviral.aa_subs as aa_subs_pb2 {
      input:
        alignment = mafft_pb2.msa,
        protein_name = "PB2"
    }
    call flu_antiviral.antiviral_mutations_parser as flu_antiviral_parser_pb2 {
      input:
        mutations_tsv = aa_subs_pb2.aa_changes_tsv,
        antiviral_aa_subs = antiviral_aa_subs
    }
  }
  call flu_antiviral.serialization {
    input:
      flu_A_315675_resistance_array = select_all([flu_antiviral_parser_h1_ha.A_315675_aa_subs, flu_antiviral_parser_h3_ha.A_315675_aa_subs, flu_antiviral_parser_n1_na.A_315675_aa_subs, flu_antiviral_parser_n2_na.A_315675_aa_subs, flu_antiviral_parser_pa.A_315675_aa_subs, flu_antiviral_parser_pb1.A_315675_aa_subs, flu_antiviral_parser_pb2.A_315675_aa_subs]),
      flu_amantadine_resistance_array = select_all([flu_antiviral_parser_h1_ha.amantadine_aa_subs, flu_antiviral_parser_h3_ha.amantadine_aa_subs, flu_antiviral_parser_n1_na.amantadine_aa_subs, flu_antiviral_parser_n2_na.amantadine_aa_subs, flu_antiviral_parser_pa.amantadine_aa_subs, flu_antiviral_parser_pb1.amantadine_aa_subs, flu_antiviral_parser_pb2.amantadine_aa_subs, flu_antiviral_parser_h1n1_mp.amantadine_aa_subs, flu_antiviral_parser_h3n2_mp.amantadine_aa_subs]),
      flu_compound_367_resistance_array = select_all([flu_antiviral_parser_h1_ha.compound_367_aa_subs, flu_antiviral_parser_h3_ha.compound_367_aa_subs, flu_antiviral_parser_n1_na.compound_367_aa_subs, flu_antiviral_parser_n2_na.compound_367_aa_subs, flu_antiviral_parser_pa.compound_367_aa_subs, flu_antiviral_parser_pb1.compound_367_aa_subs, flu_antiviral_parser_pb2.compound_367_aa_subs]),
      flu_favipiravir_resistance_array = select_all([flu_antiviral_parser_h1_ha.favipiravir_aa_subs, flu_antiviral_parser_h3_ha.favipiravir_aa_subs, flu_antiviral_parser_n1_na.favipiravir_aa_subs, flu_antiviral_parser_n2_na.favipiravir_aa_subs, flu_antiviral_parser_pa.favipiravir_aa_subs, flu_antiviral_parser_pb1.favipiravir_aa_subs, flu_antiviral_parser_pb2.favipiravir_aa_subs]),
      flu_fludase_resistance_array = select_all([flu_antiviral_parser_h1_ha.fludase_aa_subs, flu_antiviral_parser_h3_ha.fludase_aa_subs, flu_antiviral_parser_n1_na.fludase_aa_subs, flu_antiviral_parser_n2_na.fludase_aa_subs, flu_antiviral_parser_pa.fludase_aa_subs, flu_antiviral_parser_pb1.fludase_aa_subs, flu_antiviral_parser_pb2.fludase_aa_subs]),
      flu_L_742_001_resistance_array = select_all([flu_antiviral_parser_h1_ha.L_742_001_aa_subs, flu_antiviral_parser_h3_ha.L_742_001_aa_subs, flu_antiviral_parser_n1_na.L_742_001_aa_subs, flu_antiviral_parser_n2_na.L_742_001_aa_subs, flu_antiviral_parser_pa.L_742_001_aa_subs, flu_antiviral_parser_pb1.L_742_001_aa_subs, flu_antiviral_parser_pb2.L_742_001_aa_subs]),
      flu_laninamivir_resistance_array = select_all([flu_antiviral_parser_h1_ha.laninamivir_aa_subs, flu_antiviral_parser_h3_ha.laninamivir_aa_subs, flu_antiviral_parser_n1_na.laninamivir_aa_subs, flu_antiviral_parser_n2_na.laninamivir_aa_subs, flu_antiviral_parser_pa.laninamivir_aa_subs, flu_antiviral_parser_pb1.laninamivir_aa_subs, flu_antiviral_parser_pb2.laninamivir_aa_subs]),
      flu_peramivir_resistance_array = select_all([flu_antiviral_parser_h1_ha.peramivir_aa_subs, flu_antiviral_parser_h3_ha.peramivir_aa_subs, flu_antiviral_parser_n1_na.peramivir_aa_subs, flu_antiviral_parser_n2_na.peramivir_aa_subs, flu_antiviral_parser_pa.peramivir_aa_subs, flu_antiviral_parser_pb1.peramivir_aa_subs, flu_antiviral_parser_pb2.peramivir_aa_subs]),
      flu_pimodivir_resistance_array = select_all([flu_antiviral_parser_h1_ha.pimodivir_aa_subs, flu_antiviral_parser_h3_ha.pimodivir_aa_subs, flu_antiviral_parser_n1_na.pimodivir_aa_subs, flu_antiviral_parser_n2_na.pimodivir_aa_subs, flu_antiviral_parser_pa.pimodivir_aa_subs, flu_antiviral_parser_pb1.pimodivir_aa_subs, flu_antiviral_parser_pb2.pimodivir_aa_subs]),
      flu_rimantadine_resistance_array = select_all([flu_antiviral_parser_h1_ha.rimantadine_aa_subs, flu_antiviral_parser_h3_ha.rimantadine_aa_subs, flu_antiviral_parser_n1_na.rimantadine_aa_subs, flu_antiviral_parser_n2_na.rimantadine_aa_subs, flu_antiviral_parser_pa.rimantadine_aa_subs, flu_antiviral_parser_pb1.rimantadine_aa_subs, flu_antiviral_parser_pb2.rimantadine_aa_subs, flu_antiviral_parser_h1n1_mp.rimantadine_aa_subs, flu_antiviral_parser_h3n2_mp.rimantadine_aa_subs]),
      flu_oseltamivir_resistance_array = select_all([flu_antiviral_parser_h1_ha.oseltamivir_aa_subs, flu_antiviral_parser_h3_ha.oseltamivir_aa_subs, flu_antiviral_parser_n1_na.oseltamivir_aa_subs, flu_antiviral_parser_n2_na.oseltamivir_aa_subs, flu_antiviral_parser_pa.oseltamivir_aa_subs, flu_antiviral_parser_pb1.oseltamivir_aa_subs, flu_antiviral_parser_pb2.oseltamivir_aa_subs]),
      flu_xofluza_resistance_array = select_all([flu_antiviral_parser_h1_ha.xofluza_aa_subs, flu_antiviral_parser_h3_ha.xofluza_aa_subs, flu_antiviral_parser_n1_na.xofluza_aa_subs, flu_antiviral_parser_n2_na.xofluza_aa_subs, flu_antiviral_parser_pa.xofluza_aa_subs, flu_antiviral_parser_pb1.xofluza_aa_subs, flu_antiviral_parser_pb2.xofluza_aa_subs]),
      flu_zanamivir_resistance_array = select_all([flu_antiviral_parser_h1_ha.zanamivir_aa_subs, flu_antiviral_parser_h3_ha.zanamivir_aa_subs, flu_antiviral_parser_n1_na.zanamivir_aa_subs, flu_antiviral_parser_n2_na.zanamivir_aa_subs, flu_antiviral_parser_pa.zanamivir_aa_subs, flu_antiviral_parser_pb1.zanamivir_aa_subs, flu_antiviral_parser_pb2.zanamivir_aa_subs])
  }
  output {
    # Flu antiviral resistance mutations
    String flu_A_315675_resistance = serialization.flu_A_315675_resistance
    String flu_amantadine_resistance = serialization.flu_amantadine_resistance
    String flu_compound_367_resistance = serialization.flu_compound_367_resistance
    String flu_favipiravir_resistance = serialization.flu_favipiravir_resistance
    String flu_fludase_resistance = serialization.flu_fludase_resistance
    String flu_L_742_001_resistance = serialization.flu_L_742_001_resistance
    String flu_laninamivir_resistance = serialization.flu_laninamivir_resistance
    String flu_peramivir_resistance = serialization.flu_peramivir_resistance
    String flu_pimodivir_resistance = serialization.flu_pimodivir_resistance
    String flu_rimantadine_resistance = serialization.flu_rimantadine_resistance
    String flu_oseltamivir_resistance = serialization.flu_oseltamivir_resistance
    String flu_xofluza_resistance = serialization.flu_xofluza_resistance
    String flu_zanamivir_resistance = serialization.flu_zanamivir_resistance
   
    # AA sunstitutions for various flu segments
    File? flu_h1_ha_aa_subs = aa_subs_h1_ha.aa_changes_tsv
    File? flu_h3_ha_aa_subs = aa_subs_h3_ha.aa_changes_tsv
    File? flu_n1_na_aa_subs = aa_subs_n1_na.aa_changes_tsv
    File? flu_n2_na_aa_subs = aa_subs_n2_na.aa_changes_tsv
    File? flu_pa_aa_subs = aa_subs_pa.aa_changes_tsv
    File? flu_pb1_aa_subs = aa_subs_pb1.aa_changes_tsv
    File? flu_pb2_aa_subs = aa_subs_pb2.aa_changes_tsv
    File? flu_h1n1_mp_aa_subs = aa_subs_h1n1_mp.aa_changes_tsv
    File? flu_h3n2_mp_aa_subs = aa_subs_h3n2_mp.aa_changes_tsv
  }
}