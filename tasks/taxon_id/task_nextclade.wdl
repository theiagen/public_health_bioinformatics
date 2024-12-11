version 1.0

task nextclade {
  meta {
    description: "Nextclade classification of one sample. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
  }
  input {
    File genome_fasta
    File? root_sequence
    File? auspice_reference_tree_json
    File? qc_config_json
    File? gene_annotations_gff
    File? pcr_primers_csv
    File? virus_properties
    String docker = "us-docker.pkg.dev/general-theiagen/nextstrain/nextclade:2.14.0"
    String dataset_name
    String dataset_reference
    String dataset_tag
    Int disk_size = 50
    Int memory = 4
    Int cpu = 2
  }
  String basename = basename(genome_fasta, ".fasta")
  command <<<
    NEXTCLADE_VERSION="$(nextclade --version)"
    echo $NEXTCLADE_VERSION > NEXTCLADE_VERSION

    nextclade dataset get --name="~{dataset_name}" --reference="~{dataset_reference}" --tag="~{dataset_tag}" -o nextclade_dataset_dir --verbose
    set -e
    nextclade run \
        --input-dataset=nextclade_dataset_dir/ \
        ~{"--input-root-seq " + root_sequence} \
        ~{"--input-tree " + auspice_reference_tree_json} \
        ~{"--input-qc-config " + qc_config_json} \
        ~{"--input-gene-map " + gene_annotations_gff} \
        ~{"--input-pcr-primers " + pcr_primers_csv} \
        ~{"--input-virus-properties " + virus_properties}  \
        --output-json "~{basename}".nextclade.json \
        --output-tsv  "~{basename}".nextclade.tsv \
        --output-tree "~{basename}".nextclade.auspice.json \
        --output-all=. \
        "~{genome_fasta}"
  >>>
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3 
  }
  output {
    String nextclade_version = read_string("NEXTCLADE_VERSION")
    File nextclade_json = "~{basename}.nextclade.json"
    File auspice_json = "~{basename}.nextclade.auspice.json"
    File nextclade_tsv = "~{basename}.nextclade.tsv"
    String nextclade_docker = docker
    String nextclade_dataset_tag = "~{dataset_tag}"
  }
}

task nextclade_v3 {
  meta {
    description: "Nextclade classification of one sample. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
  }
  input {
    File genome_fasta
    File? auspice_reference_tree_json
    File? gene_annotations_gff
    File? nextclade_pathogen_json
    File? input_ref
    String docker = "us-docker.pkg.dev/general-theiagen/nextstrain/nextclade:3.3.1" 
    String dataset_name
    String verbosity = "warn" # other options are: "off" "error" "info" "debug" and "trace"
    String dataset_tag
    Int disk_size = 50
    Int memory = 4
    Int cpu = 2
  }
  String basename = basename(genome_fasta, ".fasta")
  command <<<
    # track version & print to log
    nextclade --version | tee NEXTCLADE_VERSION

    # --reference no longer used in v3. consolidated into --name and --tag
    nextclade dataset get \
      --name="~{dataset_name}" \
      --tag="~{dataset_tag}" \
      -o nextclade_dataset_dir \
      --verbosity ~{verbosity}

    # exit script/task upon error
    set -e

    # not necessary to include `--jobs <jobs>` in v3. Nextclade will use all available CPU threads by default. It's fast so I don't think we will need to change unless we see errors
    nextclade run \
      --input-dataset nextclade_dataset_dir/ \
      ~{"--input-ref " + input_ref} \
      ~{"--input-tree " + auspice_reference_tree_json} \
      ~{"--input-pathogen-json " + nextclade_pathogen_json} \
      ~{"--input-annotation " + gene_annotations_gff} \
      --output-json "~{basename}".nextclade.json \
      --output-tsv  "~{basename}".nextclade.tsv \
      --output-tree "~{basename}".nextclade.auspice.json \
      --output-all . \
      --verbosity ~{verbosity} \
      "~{genome_fasta}"
  >>>
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3 
  }
  output {
    String nextclade_version = read_string("NEXTCLADE_VERSION")
    File nextclade_json = "~{basename}.nextclade.json"
    File auspice_json = "~{basename}.nextclade.auspice.json"
    File nextclade_tsv = "~{basename}.nextclade.tsv"
    String nextclade_docker = docker
    String nextclade_dataset_tag = "~{dataset_tag}"
  }
}

task nextclade_output_parser {
  meta {
    description: "Python and bash codeblocks for parsing the output files from Nextclade."
  }
  input {
    File nextclade_tsv
    String docker = "us-docker.pkg.dev/general-theiagen/python/python:3.8.18-slim"
    Int disk_size = 50
    Int memory = 4
    Int cpu = 2
    String? organism
  }
  command <<<
    # Parse outputs using python3
    python3 <<CODE
    import csv

    with open("~{nextclade_tsv}", 'r') as tsv_file:
      tsv_reader = csv.reader(tsv_file, delimiter="\t")
      tsv_data = list(tsv_reader)

      if len(tsv_data) == 1:
        tsv_data.append(['NA']*len(tsv_data[0]))

      tsv_dict = dict(zip(tsv_data[0], tsv_data[1]))

      # function to write a field in the tsv_dict to a file for output
      def write_field_to_file(output_file_name, item_to_parse):
        with open(output_file_name, 'wt') as output_file:
          item = tsv_dict[item_to_parse]
          if item == '':
            item = 'NA'
          output_file.write(item)

      # combine 'clade_nextstrain' and 'clade_who' column if sars-cov-2
      if ("~{organism}" == "sars-cov-2"):
        with open("NEXTCLADE_CLADE", 'wt') as nextclade_clade:
          nc_clade = tsv_dict['clade_nextstrain']
          who_clade = tsv_dict['clade_who']
          if (nc_clade != who_clade) and (nc_clade != '') and (who_clade != ''):
            nc_clade = nc_clade + " (" + who_clade + ")"
          if nc_clade == '':
            nc_clade = 'NA'
          nextclade_clade.write(nc_clade)

      else:
        write_field_to_file("NEXTCLADE_CLADE", 'clade')

      write_field_to_file('NEXTCLADE_AASUBS', 'aaSubstitutions')
      write_field_to_file('NEXTCLADE_AADELS', 'aaDeletions')

      if 'lineage' in tsv_dict:
        write_field_to_file('NEXTCLADE_LINEAGE', 'lineage')
      elif 'Nextclade_pango' in tsv_dict:
        write_field_to_file('NEXTCLADE_LINEAGE', 'Nextclade_pango')

      write_field_to_file('NEXTCLADE_QC', 'qc.overallStatus')

      if ("~{organism}" == "flu"):
        # split the amino acid mutations by segment
        segments=["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

        def process_nc_aa_string(file_path):
          segments_dict = {segment: [] for segment in segments}

          with open(file_path, 'r') as file:
            nc_aa_string = file.read().strip()

          mutation_list = nc_aa_string.split(',')

          for single_mutation in mutation_list:
            segment, change = single_mutation.split(':')
            if segment in segments_dict and change != '':
              segments_dict[segment].append(single_mutation)

        if nc_aa_subs != '' and nc_aa_subs != 'NA':
          process_nc_aa_string("NEXTCLADE_AASUBS")

        if nc_aa_dels != '' and nc_aa_dels != 'NA':
          process_nc_aa_string("NEXTCLADE_AADELS")

        write_field_to_file
        for segment, values in segments_dict.items():
          with open(f"NEXTCLADE_AA_FLU_{segment}", 'w') as nextclade_aa_flu:
            nextclade_aa_flu.write(','.join(values))])
      else:
       # prevent WDL failures
       touch NEXTCLADE_AA_FLU_PB2 NEXTCLADE_AA_FLU_PB1 NEXTCLADE_AA_FLU_PA NEXTCLADE_AA_FLU_HA NEXTCLADE_AA_FLU_NP NEXTCLADE_AA_FLU_NA NEXTCLADE_AA_FLU_MP NEXTCLADE_AA_FLU_NS

    CODE
  >>>
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
  output {
    String nextclade_clade = read_string("NEXTCLADE_CLADE")
    String nextclade_aa_subs = read_string("NEXTCLADE_AASUBS")
    String nextclade_aa_dels = read_string("NEXTCLADE_AADELS")
    String nextclade_lineage = read_string("NEXTCLADE_LINEAGE")
    String nextclade_qc = read_string("NEXTCLADE_QC")
    # flu fields only
    String nextclade_aa_flu_pb2 = read_string("NEXTCLADE_AA_FLU_PB2")
    String nextclade_aa_flu_pb1 = read_string("NEXTCLADE_AA_FLU_PB1")
    String nextclade_aa_flu_pa = read_string("NEXTCLADE_AA_FLU_PA")
    String nextclade_aa_flu_ha = read_string("NEXTCLADE_AA_FLU_HA")
    String nextclade_aa_flu_np = read_string("NEXTCLADE_AA_FLU_NP")
    String nextclade_aa_flu_na = read_string("NEXTCLADE_AA_FLU_NA")
    String nextclade_aa_flu_mp = read_string("NEXTCLADE_AA_FLU_MP")
    String nextclade_aa_flu_ns = read_string("NEXTCLADE_AA_FLU_NS")
  }
}

task nextclade_add_ref {
  meta {
    description: "Nextclade task to add samples to either a user specified or a nextclade reference tree."
  }
  input {
    File genome_fasta
    File? reference_tree_json
    File? nextclade_pathogen_json
    File? gene_annotations_gff
    File? input_ref
    String docker = "us-docker.pkg.dev/general-theiagen/nextstrain/nextclade:3.3.1"
    String dataset_name
    String? dataset_tag
    String verbosity = "warn" # other options are: "off" "error" "info" "debug" and "trace"
    Int disk_size = 100
    Int memory = 4
    Int cpu = 2
  }
  String basename = basename(genome_fasta, ".fasta")
  command <<<
    # track version & print to log
    nextclade --version | tee NEXTCLADE_VERSION

    echo "DEBUG: downloading nextclade dataset..."
    nextclade dataset get \
      --name="~{dataset_name}" \
      ~{"--tag " + dataset_tag} \
      -o nextclade_dataset_dir \
      --verbosity ~{verbosity}

    # If no reference sequence is provided, use the reference tree from the dataset
    if [ -z "~{reference_tree_json}" ]; then
      echo "Default dataset reference tree JSON will be used"
      cp -v nextclade_dataset_dir/tree.json reference_tree.json
    else
      echo "User reference tree JSON will be used"
      cp -v ~{reference_tree_json} reference_tree.json
    fi

    tree_json="reference_tree.json"

    set -e
    echo "DEBUG: running nextclade..."
    nextclade run \
      --input-dataset nextclade_dataset_dir/ \
      --input-tree ${tree_json} \
      ~{"--input-pathogen-json " + nextclade_pathogen_json} \
      ~{"--input-annotation " + gene_annotations_gff} \
      ~{"--input-ref " + input_ref} \
      --output-json "~{basename}".nextclade.json \
      --output-tsv  "~{basename}".nextclade.tsv \
      --output-tree "~{basename}".nextclade.auspice.json \
      --output-all=. \
      "~{genome_fasta}"
  >>>
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
  output {
    String nextclade_version = read_string("NEXTCLADE_VERSION")
    File nextclade_json = "~{basename}.nextclade.json"
    File auspice_json = "~{basename}.nextclade.auspice.json"
    File nextclade_tsv = "~{basename}.nextclade.tsv"
    String nextclade_docker = docker
    File netclade_ref_tree = "reference_tree.json"
  }
}