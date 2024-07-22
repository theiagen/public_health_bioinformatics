version 1.0

task freyja_one_sample {
  input {
    File primer_trimmed_bam
    String samplename
    File reference_genome
    File? freyja_usher_barcodes
    File? freyja_lineage_metadata
    Float? eps
    Float? adapt
    Boolean update_db = false
    Boolean confirmed_only = false
    Boolean bootstrap = false
    Int? number_bootstraps
    Int? depth_cutoff
    Int memory = 4
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/freyja:1.5.1-07_02_2024-01-27-2024-07-22"
    Int disk_size = 100
  }
  command <<<
  # capture version
  freyja --version | tee FREYJA_VERSION
  
  # update freyja reference files if specified
  if ~{update_db}; then 
      freyja update 2>&1 | tee freyja_update.log
      # check log files to ensure update did not fail
      if grep "FileNotFoundError.*lineagePaths.*" freyja_update.log
      then 
        echo "Error in attempting to update Freyja files. Try increasing memory"
        >&2 echo "Killed"
        exit 1
      fi
      if grep "error" freyja_update.log
      then
        grep "error" freyja_update.log | tail -1
        >&2 echo "Killed"
        exit 1
      fi
      # can't update barcodes in freyja 1.3.2; will update known issue is closed (https://github.com/andersen-lab/Freyja/issues/33)
      freyja_usher_barcode_version="freyja update: $(date +"%Y-%m-%d")"
      freyja_metadata_version="freyja update: $(date +"%Y-%m-%d")"
  else
    # configure barcode    
    if [[ ! -z "~{freyja_usher_barcodes}" ]]; then
      echo "User freyja usher barcodes identified; ~{freyja_usher_barcodes} will be utilized fre freyja demixing"
      freyja_usher_barcode_version=$(basename -- "~{freyja_usher_barcodes}")
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
  echo "Running: freyja variants ~{primer_trimmed_bam} --variants ~{samplename}_freyja_variants.tsv --depths ~{samplename}_freyja_depths.tsv --ref ~{reference_genome}"
  freyja variants \
    ~{primer_trimmed_bam} \
    --variants ~{samplename}_freyja_variants.tsv \
    --depths ~{samplename}_freyja_depths.tsv \
    --ref ~{reference_genome}
  
  # Calculate Boostraps, if specified
  if ~{bootstrap}; then
    freyja boot \
    ~{"--eps " + eps} \
    ~{"--meta " + freyja_lineage_metadata} \
    ~{"--barcodes " + freyja_usher_barcodes} \
    ~{"--depthcutoff " + depth_cutoff} \
    ~{"--nb " + number_bootstraps } \
    ~{true='--confirmedonly' false='' confirmed_only} \
    ~{samplename}_freyja_variants.tsv \
    ~{samplename}_freyja_depths.tsv \
    --output_base ~{samplename} \
    --boxplot pdf
  fi
  
  # Demix variants 
  echo "Running: freyja demix --eps ~{eps} ${freyja_barcode} ${freyja_metadata} ~{samplename}_freyja_variants.tsv ~{samplename}_freyja_depths.tsv --output ~{samplename}_freyja_demixed.tmp"
  freyja demix \
    ~{'--eps ' + eps} \
    ~{'--meta ' + freyja_lineage_metadata} \
    ~{'--barcodes ' + freyja_usher_barcodes} \
    ~{'--depthcutoff ' + depth_cutoff} \
    ~{true='--confirmedonly' false='' confirmed_only} \
    ~{'--adapt ' + adapt} \
    ~{samplename}_freyja_variants.tsv \
    ~{samplename}_freyja_depths.tsv \
    --output ~{samplename}_freyja_demixed.tmp
  
  # Adjust output header
  echo -e "\t/~{samplename}" > ~{samplename}_freyja_demixed.tsv
  tail -n+2 ~{samplename}_freyja_demixed.tmp >> ~{samplename}_freyja_demixed.tsv

  if [ -f /opt/conda/envs/freyja-env/lib/python3.10/site-packages/freyja/data/usher_barcodes.csv ]; then
    mv /opt/conda/envs/freyja-env/lib/python3.10/site-packages/freyja/data/usher_barcodes.csv usher_barcodes.csv
  fi
  
  #Output QC values to the Terra data table
  python <<CODE
  import csv
  #Want coverage output from freyja_demixed tsv file
  with open("~{samplename}_freyja_demixed.tsv",'r') as tsv_file:
    tsv_reader = csv.reader(tsv_file, delimiter="\t")
    for line in tsv_reader:
      if "coverage" in line[0]:
        with open("COVERAGE", 'wt') as coverage:
          coverage.write(line[1])

  CODE

  >>>
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: "~{docker}"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    maxRetries: 3
  }
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
    # capture barcode file - first is user supplied, second appears if the user did not supply a barcode file
    File freyja_barcode_file = select_first([freyja_usher_barcodes, "usher_barcodes.csv"])
    String freyja_barcode_version = read_string("FREYJA_BARCODES")
    String freyja_metadata_version = read_string("FREYJA_METADATA")
    String freyja_version = read_string("FREYJA_VERSION")
  }
}
