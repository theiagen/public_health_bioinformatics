version 1.0

task kleborate_klebsiella {
  # Inputs
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/kleborate:3.2.4_20251112"
    Int disk_size = 100
    Int cpu = 8
    Int memory = 16
    
    # Parameters
    String preset_organism = "kpsc" # Preset organism to run Kleborate klebsiella with {kpsc,kosc} (default: kpsc)
    Float min_percent_identity = 90.0 # Minimum alignment percent identity for main results (default: 90.0)
    Float min_percent_coverage = 80.0 #  Minimum alignment percent coverage for main results (default: 80.0)
    Float min_spurious_percent_identity = 80.0 # Minimum alignment percent identity for spurious results (default: 80.0)
    Float min_spurious_percent_coverage = 40.0 #  Minimum alignment percent coverage for spurious results (default: 40.0)
  }
  command <<<
    # Print and save date
    date | tee DATE
    
    # Print and save version
    kleborate --version | tee VERSION 

    # If preset_organism is kpsc, run kleborate for Klebsiella pneumoniae species complex
    if [ "~{preset_organism}" = "kpsc" ]; then
      
      kleborate \
      ~{'--klebsiella_pneumo_complex__mlst_min_identity ' + min_percent_identity} \
      ~{'--klebsiella_pneumo_complex__mlst_min_coverage ' + min_percent_coverage} \
      ~{'--klebsiella_pneumo_complex__amr_min_spurious_identity ' + min_spurious_percent_identity} \
      ~{'--klebsiella_pneumo_complex__amr_min_spurious_coverage ' + min_spurious_percent_coverage} \
      --outdir kleborate_results \
      --assemblies ~{assembly} \
      --preset ~{preset_organism} \
      --trim_headers

      # Check if kleborate_pneumo_complex_output.txt exists before moving and parsing
      if [ -f kleborate_results/klebsiella_pneumo_complex_output.txt ]; then

        mv kleborate_results/klebsiella_pneumo_complex_output.txt ~{samplename}_kleborate_kpsc_out.tsv
        mv kleborate_results/klebsiella_pneumo_complex_hAMRonization_output.txt ~{samplename}_hAMRonization_kleborate_out.tsv
      
      else
        #Create empty output file with message if output not found, parser will fail otherwise
        echo "No kleborate_pneumo_complex_output.txt found" >> ~{samplename}_kleborate_kpsc_out.tsv
        echo "None" >> ~{samplename}_kleborate_kpsc_out.tsv
        echo "No hAMRonization output for K. pneumoniae complex" >> ~{samplename}_hAMRonization_kleborate_out.tsv
      fi

      parse_kleborate_kleb.py ~{samplename}_kleborate_kpsc_out.tsv

    fi

    if [ "~{preset_organism}" = "kosc" ]; then
      
      kleborate \
      ~{'--klebsiella_oxytoca_complex__mlst_min_identity ' + min_percent_identity} \
      ~{'--klebsiella_oxytoca_complex__mlst_min_coverage ' + min_percent_coverage} \
      --outdir kleborate_results \
      --assemblies ~{assembly} \
      --preset ~{preset_organism} \
      --trim_headers

      # Check if kleborate_oxytoca_complex_output.txt exists before moving and parsing
      if [ -f kleborate_results/klebsiella_oxytoca_complex_output.txt ]; then

        mv kleborate_results/klebsiella_oxytoca_complex_output.txt ~{samplename}_kleborate_kosc_out.tsv

      else
        #Create empty output file with message if output not found, parser will fail otherwise
        echo "No kleborate_oxytoca_complex_output.txt found" >> ~{samplename}_kleborate_kosc_out.tsv
        echo "None" >> ~{samplename}_kleborate_kosc_out.tsv
      fi

      # Create placeholder hAMRonization file for kosc (not produced by this preset)
      echo "No hAMRonization output for K. oxytoca complex" >> ~{samplename}_hAMRonization_kleborate_out.tsv

      parse_kleborate_kosc.py ~{samplename}_kleborate_kosc_out.tsv

    fi

    
  >>>
  output {
    File kleborate_klebsiella_output_file = glob("~{samplename}_kleborate_*_out.tsv")[0]
    File kleborate_klebsiella_hAMRonization_output_file = "~{samplename}_hAMRonization_kleborate_out.tsv"
    String kleborate_klebsiella_version = read_string("VERSION")
    String kleborate_klebsiella_docker = docker
    String kleborate_klebsiella_mlst_sequence_type = read_string("MLST_SEQUENCE_TYPE")
    String kleborate_klebsiella_virulence_score = read_string("VIRULENCE_SCORE")
    String kleborate_klebsiella_resistance_score = read_string("RESISTANCE_SCORE")
    String kleborate_klebsiella_num_resistance_genes = read_string("NUM_RESISTANCE_GENES")
    String kleborate_klebsiella_bla_resistance_genes = read_string("BLA_RESISTANCE_GENES")
    String kleborate_klebsiella_esbl_resistance_genes = read_string("ESBL_RESISTANCE_GENES")
    String kleborate_klebsiella_key_resistance_genes = read_string("KEY_RESISTANCE_GENES")
    String kleborate_klebsiella_genomic_resistance_mutations = read_string("GENOMIC_RESISTANCE_MUTATIONS")
    String kleborate_klebsiella_klocus = read_string("K_LOCUS")
    String kleborate_klebsiella_ktype = read_string("K_TYPE")
    String kleborate_klebsiella_olocus = read_string("O_LOCUS")
    String kleborate_klebsiella_otype = read_string("O_TYPE")
    String kleborate_klebsiella_klocus_confidence = read_string("K_LOCUS_CONFIDENCE")
    String kleborate_klebsiella_olocus_confidence = read_string("O_LOCUS_CONFIDENCE")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}

task kleborate_ecoli {
  # Inputs
  input {
    File assembly
    String samplename
    Float min_percent_identity = 90.0 # Minimum alignment percent identity for main results (default: 90.0)
    Float min_percent_coverage = 80.0 #  Minimum alignment percent coverage for main results (default: 80.0)
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/kleborate:3.2.4_20251112"
    Int disk_size = 100
    Int cpu = 8
    Int memory = 16
  }
  command <<<
    # Print and save date
    date | tee DATE
    
    # Print and save version
    kleborate --version | tee VERSION 

    # Run Kleborate on the input assembly for Escherichia
    kleborate \
    ~{'--escherichia_mlst_achtman_min_identity ' + min_percent_identity} \
    ~{'--escherichia_mlst_achtman_min_coverage ' + min_percent_coverage} \
    --outdir kleborate_results \
    --assemblies ~{assembly} \
    --preset escherichia \
    --trim_headers

    # Check if output file exists before moving
    if [ -f kleborate_results/escherichia_output.txt ]; then
      mv kleborate_results/escherichia_output.txt ~{samplename}_escherichia_output.txt
    else
      # Create empty output file with message if output not found
      echo "No escherichia coli output found for escherichia module" > ~{samplename}_escherichia_output.txt
      echo "None" >> ~{samplename}_escherichia_output.txt
    fi
    # Parse the output for relevant fields
    parse_kleborate_ecoli.py ~{samplename}_escherichia_output.txt
  >>>
  output {
    File kleborate_ecoli_output_file = "~{samplename}_escherichia_output.txt"
    String kleborate_ecoli_version = read_string("VERSION")
    String kleborate_ecoli_lee_st = read_string("LEE_ST")
    String kleborate_ecoli_lee_lineage = read_string("LEE_LINEAGE")
    String kleborate_ecoli_pathotype = read_string("PATHOTYPE")
    String kleborate_ecoli_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}