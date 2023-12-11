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
      memory: "4 GB"
      cpu: 2
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
    }
}

task nextclade_output_parser {
    meta {
      description: "Python and bash codeblocks for parsing the output files from Nextclade."
    }
    input {
      File nextclade_tsv
      String docker = "python:slim"
      Int disk_size = 50
      String? organism
      Boolean? NA_segment
      String tamiflu_aa_substitutions = "NA:H275Y,NA:R292K"
    }
    command <<<
      # Set WDL input variable to input.tsv file
      cat "~{nextclade_tsv}" > input.tsv
      touch TAMIFLU_AASUBS
     
      # Parse outputs using python3
      python3 <<CODE
      import csv
      import codecs

      # list of aa substitutions linked with tamiflu resistance - with the possibility to
      # extend the list with other provived aa substitutions in the format "NA:V95A,NA:I97V"
      tamiflu_aa_subs = ["NA:V95A","NA:I97V","NA:E99A","NA:H101L","NA:G108E",
      "NA:Q116L","NA:V116A","NA:E119D","NA:E119G","NA:E119I","NA:E119V","NA:R136K",
      "NA:T146K","NA:T146P","NA:D151E","NA:N169S","NA:D179N","NA:D197N","NA:D198E",
      "NA:D198G","NA:D198N","NA:A200T","NA:I203M","NA:I203R","NA:I203V","NA:I221T",
      "NA:I222R","NA:I222V","NA:I223R","NA:I223V","NA:S227N","NA:S247N","NA:H255Y",
      "NA:E258Q","NA:H274N","NA:H274Y","NA:H275Y","NA:N275S","NA:H277Y","NA:R292K",
      "NA:N294S","NA:S334N","NA:R371K","NA:D432G","NA:H439P","NA:H439R"]

      tamiflu_aa_subs = tamiflu_aa_subs + "~{tamiflu_aa_substitutions}".split(',')

      def intersection(lst1, lst2):
        # returns intersection between nextclade identified aa substitutions and
        # tamiflu associated aa substitutions
        return list(set(lst1) & set(lst2))

      with codecs.open("./input.tsv",'r') as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        tsv_data = list(tsv_reader)

        if len(tsv_data) == 1:
          tsv_data.append(['NA']*len(tsv_data[0]))
        tsv_dict = dict(zip(tsv_data[0], tsv_data[1]))

        # combine 'clade_nextstrain' and 'clade_who' column if sars-cov-2, if false then parse 'clade' column
        if ("~{organism}" == "sars-cov-2"):
          with codecs.open("NEXTCLADE_CLADE", 'wt') as Nextclade_Clade:
            nc_clade = tsv_dict['clade_nextstrain']
            who_clade = tsv_dict['clade_who']
            if (nc_clade != who_clade) and (nc_clade != '') and (who_clade != ''):
              nc_clade = nc_clade + " (" + who_clade + ")"
            if nc_clade == '':
              nc_clade = 'NA'
            Nextclade_Clade.write(nc_clade)
        else:
          with codecs.open("NEXTCLADE_CLADE", 'wt') as Nextclade_Clade:
            nc_clade = tsv_dict['clade']
            if nc_clade == '':
              nc_clade = 'NA'
            Nextclade_Clade.write(nc_clade)

        with codecs.open("NEXTCLADE_AASUBS", 'wt') as Nextclade_AA_Subs:
          nc_aa_subs = tsv_dict['aaSubstitutions']
          if nc_aa_subs == '':
            nc_aa_subs = 'NA'
          else:
            # if organism is flu, return list of aa subs associated with tamiflu resistance
            if ("~{organism}" == "flu" and "~{NA_segment}" == "true"):
              tamiflu_subs = intersection(tamiflu_aa_subs, nc_aa_subs.split(','))
              with codecs.open("TAMIFLU_AASUBS", 'wt') as Tamiflu_AA_Subs:
                Tamiflu_AA_Subs.write(",".join(tamiflu_subs))
          Nextclade_AA_Subs.write(nc_aa_subs)

        with codecs.open("NEXTCLADE_AADELS", 'wt') as Nextclade_AA_Dels:
          nc_aa_dels = tsv_dict['aaDeletions']
          if nc_aa_dels == '':
            nc_aa_dels = 'NA'
          Nextclade_AA_Dels.write(nc_aa_dels)

        with codecs.open("NEXTCLADE_LINEAGE", 'wt') as Nextclade_Lineage:
          if 'lineage' in tsv_dict:
            nc_lineage = tsv_dict['lineage']
            if nc_lineage is None:
              nc_lineage = ""
          elif 'Nextclade_pango' in tsv_dict:
            nc_lineage = tsv_dict['Nextclade_pango']
            if nc_lineage is None:
              nc_lineage = ""
          else:
            nc_lineage = ""
          Nextclade_Lineage.write(nc_lineage)
        
        with codecs.open("NEXTCLADE_QC", 'wt') as Nextclade_QC:
          nc_qc = tsv_dict['qc.overallStatus']
          if nc_qc == '':
            nc_qc = 'NA'
          Nextclade_QC.write(nc_qc)
      CODE
    >>>
    runtime {
      docker: "~{docker}"
      memory: "4 GB"
      cpu: 2
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB"
      dx_instance_type: "mem1_ssd1_v2_x2"
      maxRetries: 3
    }
    output {
      String nextclade_clade = read_string("NEXTCLADE_CLADE")
      String nextclade_aa_subs = read_string("NEXTCLADE_AASUBS")
      String nextclade_tamiflu_aa_subs = read_string("TAMIFLU_AASUBS")
      String nextclade_aa_dels = read_string("NEXTCLADE_AADELS")
      String nextclade_lineage = read_string("NEXTCLADE_LINEAGE")
      String nextclade_qc = read_string("NEXTCLADE_QC")
    }
}

task nextclade_add_ref {
    meta {
      description: "Nextclade task to add samples to either a user specified or a nextclade reference tree."
    }
    input {
      File genome_fasta
      File? root_sequence
      File? reference_tree_json
      File? qc_config_json
      File? gene_annotations_gff
      File? pcr_primers_csv
      File? virus_properties
      String docker = "us-docker.pkg.dev/general-theiagen/nextstrain/nextclade:2.14.0"
      String dataset_name
      String? dataset_reference
      String? dataset_tag
      Int disk_size = 50
    }
    String basename = basename(genome_fasta, ".fasta")
    command <<<
        NEXTCLADE_VERSION="$(nextclade --version)"
        echo $NEXTCLADE_VERSION > NEXTCLADE_VERSION

        nextclade dataset get \
          --name="~{dataset_name}" \
          ~{"--reference " + dataset_reference} \
          ~{"--tag " + dataset_tag} \
          -o nextclade_dataset_dir \
          --verbose

        # If no referece sequence is provided, use the reference tree from the dataset
        if [ -z "~{reference_tree_json}" ]; then
          echo "Default dataset reference tree JSON will be used"
          cp nextclade_dataset_dir/tree.json reference_tree.json
        else
          echo "User reference tree JSON will be used"
          cp ~{reference_tree_json} reference_tree.json
        fi

        tree_json="reference_tree.json"

        set -e
        nextclade run \
            --input-dataset=nextclade_dataset_dir/ \
            ~{"--input-root-seq " + root_sequence} \
            --input-tree ${tree_json} \
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
      memory: "8 GB"
      cpu: 2
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
      File netclade_ref_tree = "reference_tree.json"
    }
}