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
    String nextclade_dataset_tag = "~{dataset_tag}"
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
    # Place holder for user defined aa substitutions
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

    # list of aa substitutions linked with tamiflu resistance
    tamiflu_aa_subs = ["HA:A28T","HA:D222G","HA:D225G","HA:G155E","HA:G158E","HA:K130E","HA:K133E","HA:K134E","HA:K140E","HA:K144E","HA:K153E","HA:K156E","HA:K234Q","HA:K238Q","HA:P194L","HA:R188K","HA:R192K","HA:R453M","HA:S138A","HA:T82K","HA:T92K","HA:V132A","HA:V135A","NA:A200T","NA:A201A","NA:A201T","NA:A245T","NA:A246T","NA:A272V","NA:A390E","NA:A395E","NA:D151D","NA:D151E","NA:D151G","NA:D151N","NA:D179G","NA:D179N","NA:D197E","NA:D197N","NA:D197Y","NA:D198D","NA:D198E","NA:D198G","NA:D198N","NA:D198Y","NA:D199E","NA:D199G","NA:D199N","NA:D213G","NA:D214G","NA:D344N","NA:D345N","NA:D432G","NA:E59G","NA:E99A","NA:E99D","NA:E99G","NA:E115V","NA:E117A","NA:E117D","NA:E117G","NA:E117V","NA:E118V","NA:E119A","NA:E119D","NA:E119E","NA:E119G","NA:E119I","NA:E119K","NA:E119V","NA:E222V","NA:E258Q","NA:E272Q","NA:E273Q","NA:E276D","NA:EG222","NA:G104E","NA:G108E","NA:G109E","NA:G140R","NA:G142R","NA:G145R","NA:G147E","NA:G147R","NA:G320E","NA:H101L","NA:H255Y","NA:H271Y","NA:H273Y","NA:H274H","NA:H274N","NA:H274Y","NA:H275H","NA:H275Y","NA:H276Y","NA:H277Y","NA:H439P","NA:H439R","NA:I97V","NA:I117M","NA:I117V","NA:I203L","NA:I203M","NA:I203R","NA:I203T","NA:I203V","NA:I219K","NA:I219L","NA:I219M","NA:I219R","NA:I219T","NA:I221L","NA:I221N","NA:I221T","NA:I222K","NA:I222L","NA:I222M","NA:I222N","NA:I222R","NA:I222T","NA:I222V","NA:I223K","NA:I223L","NA:I223M","NA:I223R","NA:I223T","NA:I223V","NA:I294V","NA:I314V","NA:I427T","NA:K130N","NA:K150N","NA:K273Q","NA:M372K","NA:M375K","NA:N44S","NA:N46S","NA:N142S","NA:N144K","NA:N146K","NA:N169S","NA:N199S","NA:N200S","NA:N220K","NA:N221K","NA:N275S","NA:N294S","NA:N295S","NA:N325K","NA:N329K","NA:N368K","NA:N369K","NA:N386K","NA:N390K","NA:P139S","NA:P141S","NA:Q116L","NA:Q136K","NA:Q136L","NA:Q313R","NA:R136K","NA:R150K","NA:R151W","NA:R152K","NA:R152W","NA:R193G","NA:R194G","NA:R221Q","NA:R222Q","NA:R224K","NA:R289K","NA:R290K","NA:R292K","NA:R293K","NA:R371K","NA:R374K","NA:S219T","NA:S227N","NA:S245N","NA:S246G","NA:S246N","NA:S246R","NA:S247G","NA:S247N","NA:S247P","NA:S247R","NA:S331R","NA:S334N","NA:S336N","NA:T146K","NA:T146P","NA:T148I","NA:T156I","NA:T157I","NA:V95A","NA:V96A","NA:V116A","NA:V215I","NA:V233M","NA:V234M","NA:V240I","NA:V241I","NA:Y142H","NA:Y144H","NA:Y155H"]
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