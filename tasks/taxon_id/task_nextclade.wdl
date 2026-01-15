version 1.0

task nextclade_v3 {
  meta {
    description: "Nextclade classification of one sample. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
  }
  input {
    File genome_fasta
    File? custom_input_dataset
    File? auspice_reference_tree_json
    File? gene_annotations_gff
    File? nextclade_pathogen_json
    File? input_ref
    String docker = "us-docker.pkg.dev/general-theiagen/nextstrain/nextclade:3.16.0" 
    String? dataset_name
    String verbosity = "warn" # other options are: "off" "error" "info" "debug" and "trace"
    String? dataset_tag
    Int disk_size = 50
    Int memory = 4
    Int cpu = 2
  }
  String basename = basename(genome_fasta, ".fasta")
  command <<<
    # exit script/task upon error
    set -euo pipefail

    # track version & print to log
    nextclade --version | tee NEXTCLADE_VERSION

    # --reference no longer used in v3. consolidated into --name and --tag
    # if a custom input dataset is not provided, then use the dataset name and tag
    # ! -s (if file is not empty) is used to check if that dataset exists and is not empty
    if [ ! -s "~{custom_input_dataset}" ]; then
      nextclade dataset get \
        --name="~{dataset_name}" \
        --tag="~{dataset_tag}" \
        -o nextclade_dataset_dir \
        --verbosity ~{verbosity}
    fi

    # not necessary to include `--jobs <jobs>` in v3. Nextclade will use all available CPU threads by default. It's fast so I don't think we will need to change unless we see errors
    nextclade run \
      --input-dataset ~{default="nextclade_dataset_dir/" custom_input_dataset} \
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
    description: "Python codeblocks for parsing the output files from Nextclade."
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
          try:
            item = tsv_dict[item_to_parse]
            if item == '':
              item = 'NA'
          except:
            item = 'NA'
          output_file.write(item)

      # combine 'clade_nextstrain' and 'clade_who' column if sars-cov-2
      # this one is slightly more complicated so the function doesn't apply
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

      write_field_to_file('NEXTCLADE_LINEAGE', 'lineage')
      if 'Nextclade_pango' in tsv_dict:
        write_field_to_file('NEXTCLADE_LINEAGE', 'Nextclade_pango')

      write_field_to_file('NEXTCLADE_QC', 'qc.overallStatus')

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
  }
}

task nextclade_v3_set {
  meta {
    description: "Nextclade task to add samples to either a user specified or a nextclade reference tree."
  }
  input {
    Array[File] genome_fastas
    File? reference_tree_json
    File? pathogen_json
    File? gene_annotations_gff
    File? input_ref
    String docker = "us-docker.pkg.dev/general-theiagen/nextstrain/nextclade:3.14.5"
    String? dataset_name
    String? dataset_tag
    String verbosity = "warn" # other options are: "off" "error" "info" "debug" and "trace"
    Int disk_size = 100
    Int memory = 4
    Int cpu = 2
  }
  command <<<
    # track version & print to log
    nextclade --version | tee NEXTCLADE_VERSION

    # make the directory incase the dataset doesnt exist
    mkdir nextclade_dataset_dir/

    if [ ! -z "~{dataset_name}" ]; then
      echo "DEBUG: downloading nextclade dataset..."
      nextclade dataset get \
        --name="~{dataset_name}" \
        ~{"--tag " + dataset_tag} \
        -o nextclade_dataset_dir \
        --verbosity ~{verbosity}
    fi
    

    # If no reference sequence is provided, use the reference tree from the dataset
    if [ -z "~{reference_tree_json}" ]; then
      echo "Default dataset reference tree JSON will be used"

      # Identify the reference tree JSON by excluding the pathogen.json
      for file in $(find nextclade_dataset_dir/ -type f -name "*.json" ! -name "*pathogen.json"); do
        # if there is a tree.json (per Nextclade documentation), use that as the reference tree
        if [[ $file == *"tree.json" ]]; then
          reference_tree_json=$file
          break
        # otherwise, if there is a file that is not pathogen.json, use that as the reference tree
        else
          reference_tree_json=$file
        fi
      done

      if [[ ! -v reference_tree_json ]]; then
        echo "ERROR: No reference tree JSON found and no user reference tree provided"
        exit 1
      fi
      cp -v $reference_tree_json reference_tree.json
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
      ~{"--input-pathogen-json " + pathogen_json} \
      ~{"--input-annotation " + gene_annotations_gff} \
      ~{"--input-ref " + input_ref} \
      --output-json nextclade_output.json \
      --output-tsv  nextclade_output.tsv \
      --output-tree nextclade_output.auspice.json \
      --output-all=. \
      ~{sep=' ' genome_fastas}
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
    File nextclade_json = "nextclade_output.json"
    File auspice_json = "nextclade_output.auspice.json"
    File nextclade_tsv = "nextclade_output.tsv"
    String nextclade_docker = docker
    File nextclade_ref_tree = "reference_tree.json"
  }
}