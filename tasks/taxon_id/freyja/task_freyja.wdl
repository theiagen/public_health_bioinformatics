version 1.0

task freyja_one_sample {
  input {
    File bamfile
    String samplename
    File reference_genome
    File? reference_gff
    String? freyja_pathogen
    File? freyja_barcodes
    File? freyja_lineage_metadata
    Boolean auto_adapt = false
    Float eps = 0.001 # set to mirror v2.0.1 default
    Float adapt = 0.0 # set to mirror v2.0.1 default
    Boolean update_db = false
    Boolean confirmed_only = false
    Boolean bootstrap = false
    Int number_bootstraps = 100 # set to mirror v1.5.3 default
    Int? depth_cutoff
    Int memory = 8
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/freyja:2.0.1"
    Int disk_size = 100
  }
  command <<<
  # capture version
  freyja --version | tee FREYJA_VERSION
  
  # update freyja reference files if specified
  if ~{update_db}; then 
      freyja update ~{"--pathogen " + freyja_pathogen} 2>&1 | tee freyja_update.log
      # check log files to ensure update did not fail
      if grep "FileNotFoundError.*lineagePaths.*" freyja_update.log; then 
        echo "Error in attempting to update Freyja files. Try increasing memory"
        >&2 echo "Killed"
        exit 1
      fi
      if grep "error" freyja_update.log; then
        grep "error" freyja_update.log | tail -1
        >&2 echo "Killed"
        exit 1
      fi

      freyja_usher_barcode_version="freyja update: $(date +"%Y-%m-%d")"
      freyja_metadata_version="freyja update: $(date +"%Y-%m-%d")"
  else
    # configure barcode    
    if [[ ! -z "~{freyja_barcodes}" ]]; then
      echo "User freyja usher barcodes identified; ~{freyja_barcodes} will be utilized for freyja demixing"
      freyja_usher_barcode_version=$(basename -- "~{freyja_barcodes}")
    else
      freyja_usher_barcode_version="unmodified from freyja container: ~{docker}"  
    fi
    # configure lineage metadata
    if [[ ! -z "~{freyja_lineage_metadata}" ]]; then
      echo "User lineage metadata; ~{freyja_lineage_metadata} will be utilized fre freyja demixing"
      freyja_metadata_version=$(basename -- "~{freyja_lineage_metadata}")
    else
      freyja_metadata_version="unmodified from freyja container: ~{docker}"
    fi
  fi
  
  # Capture reference file versions
  echo ${freyja_usher_barcode_version} | tee FREYJA_BARCODES
  echo ${freyja_metadata_version} | tee FREYJA_METADATA
  
  # Call variants and capture sequencing depth information
  echo "Running: freyja variants ~{bamfile} --variants ~{samplename}_freyja_variants.tsv --depths ~{samplename}_freyja_depths.tsv --ref ~{reference_genome}"
  freyja variants \
    ~{bamfile} \
    ~{"--annot " + reference_gff} \
    --variants ~{samplename}_freyja_variants.tsv \
    --depths ~{samplename}_freyja_depths.tsv \
    --ref ~{reference_genome}
  
  # Calculate Boostraps, if specified
  if ~{bootstrap}; then
    freyja boot \
    ~{"--pathogen " + freyja_pathogen} \
    ~{"--eps " + eps} \
    ~{"--meta " + freyja_lineage_metadata} \
    ~{"--barcodes " + freyja_barcodes} \
    ~{"--depthcutoff " + depth_cutoff} \
    ~{"--nb " + number_bootstraps } \
    ~{true='--confirmedonly' false='' confirmed_only} \
    ~{true='--autoadapt' false='' auto_adapt} \
    ~{samplename}_freyja_variants.tsv \
    ~{samplename}_freyja_depths.tsv \
    --output_base ~{samplename} \
    --boxplot pdf
  fi
  
  # Demix variants 
  echo "Running: freyja demix --eps ~{eps} ${freyja_barcode} ${freyja_metadata} ~{samplename}_freyja_variants.tsv ~{samplename}_freyja_depths.tsv --output ~{samplename}_freyja_demixed.tmp"
  freyja demix \
    ~{"--pathogen " + freyja_pathogen} \
    ~{'--eps ' + eps} \
    ~{'--meta ' + freyja_lineage_metadata} \
    ~{'--barcodes ' + freyja_barcodes} \
    ~{'--depthcutoff ' + depth_cutoff} \
    ~{true='--confirmedonly' false='' confirmed_only} \
    ~{if auto_adapt then '--autoadapt' else '--adapt ' + adapt} \
    ~{samplename}_freyja_variants.tsv \
    ~{samplename}_freyja_depths.tsv \
    --output ~{samplename}_freyja_demixed.tmp
  
  # Adjust output header
  echo -e "\t/~{samplename}" > ~{samplename}_freyja_demixed.tsv
  tail -n+2 ~{samplename}_freyja_demixed.tmp >> ~{samplename}_freyja_demixed.tsv

  if [ "~{freyja_pathogen}" == "SARS-CoV-2" ]; then
    if [ -f /opt/conda/envs/freyja-env/lib/python3.12/site-packages/freyja/data/usher_barcodes.feather ]; then
      mv /opt/conda/envs/freyja-env/lib/python3.12/site-packages/freyja/data/usher_barcodes.feather usher_barcodes.feather
    fi

    if [ -f /opt/conda/envs/freyja-env/lib/python3.12/site-packages/freyja/data/curated_lineages.json ]; then
      mv /opt/conda/envs/freyja-env/lib/python3.12/site-packages/freyja/data/curated_lineages.json curated_lineages.json
    fi
  fi

  #Output QC values to the Terra data table
  python <<CODE
  import csv
  import pandas as pd

  dataf = pd.read_csv("~{samplename}_freyja_demixed.tsv", sep="\t", header=None)
  dataf.columns = ["Attribute", "~{samplename}"]
  parsed_data = {
    "LIMS_ID": "~{samplename}"
  }
  
  #Want coverage output from freyja_demixed tsv file
  with open("~{samplename}_freyja_demixed.tsv",'r') as tsv_file:
    tsv_reader = csv.reader(tsv_file, delimiter="\t")
    for line in tsv_reader:
      if "coverage" in line[0]:
        with open("COVERAGE", 'wt') as coverage:
          coverage.write(line[1])
        parsed_data["coverage"] = dataf.loc[dataf['Attribute'] == "coverage", "~{samplename}"].values[0]
      if "lineages" in line[0]:
        with open("LINEAGES", 'wt') as lineages:
          lineages.write(line[1].replace(" ", ","))
        parsed_data["lineages"] = dataf.loc[dataf['Attribute'] == "lineages", "~{samplename}"].values[0]
      if "abundances" in line[0]:
        with open("ABUNDANCES", 'wt') as abundances:
          abundances.write(line[1].replace(" ", ","))
        parsed_data["abundances"] = dataf.loc[dataf['Attribute'] == "abundances", "~{samplename}"].values[0]
      if "resid" in line[0]:
        with open("RESID", 'wt') as resid:
          resid.write(line[1])
        parsed_data["resid"] = dataf.loc[dataf['Attribute'] == "resid", "~{samplename}"].values[0]
      if "summarized" in line[0]:
        with open("SUMMARIZED", 'wt') as summarized:
          summarized.write(line[1])
        parsed_data["summarized"] = dataf.loc[dataf['Attribute'] == "summarized", "~{samplename}"].values[0]
  
  # Initialize a list to store output rows
  output_data = []
  output_data.append(parsed_data)
  output_df = pd.DataFrame(output_data)
  output_df.to_csv("~{samplename}_freyja_demixed_parsed.tsv", sep='\t', index=False)
  
  CODE
  >>>
  output {
    File freyja_variants = "~{samplename}_freyja_variants.tsv"
    File freyja_depths = "~{samplename}_freyja_depths.tsv"
    File freyja_demixed = "~{samplename}_freyja_demixed.tsv"
    Float freyja_coverage = read_float("COVERAGE")
    File? freyja_update_log = "freyja_update.log"
    File? freyja_bootstrap_lineages = "~{samplename}_lineages.csv"
    File? freyja_bootstrap_lineages_pdf = "~{samplename}_lineages.pdf"
    File? freyja_bootstrap_summary = "~{samplename}_summarized.csv"
    File? freyja_bootstrap_summary_pdf = "~{samplename}_summarized.pdf"
    # capture barcode file if sars-cov-2
    File? freyja_sc2_barcode_file = "usher_barcodes.feather"
    File? freyja_sc2_lineage_metadata_file = "curated_lineages.json"
    String freyja_barcode_version = read_string("FREYJA_BARCODES")
    String freyja_metadata_version = read_string("FREYJA_METADATA")
    String freyja_version = read_string("FREYJA_VERSION")
    File freyja_demixed_parsed = "~{samplename}_freyja_demixed_parsed.tsv"
    String freyja_resid = read_string("RESID")
    String freyja_summarized = read_string("SUMMARIZED")
    String freyja_lineages = read_string("LINEAGES")
    String freyja_abundances = read_string("ABUNDANCES")
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: "~{docker}"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    maxRetries: 3
  }
}
