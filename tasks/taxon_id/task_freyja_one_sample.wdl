version 1.0

task freyja_one_sample {
  input {
    File primer_trimmed_bam
    String samplename
    File reference_genome
    File? freyja_usher_barcodes
    File? freyja_lineage_metadata
    Float? eps
    Boolean update_db = false
    Boolean confirmed_only = false
    Boolean bootstrap = false
    Int? number_bootstraps
    Int memory = 4
    String docker = "staphb/freyja:1.3.10"
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
    ~{"--nb " + number_bootstraps } \
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
    ~{true='--confirmedonly' false='' confirmed_only} \
    ~{samplename}_freyja_variants.tsv \
    ~{samplename}_freyja_depths.tsv \
    --output ~{samplename}_freyja_demixed.tmp
  # Adjust output header
  echo -e "\t/~{samplename}" > ~{samplename}_freyja_demixed.tsv
  tail -n+2 ~{samplename}_freyja_demixed.tmp >> ~{samplename}_freyja_demixed.tsv
  >>>
  runtime {
    memory: "~{memory} GB"
    cpu: 2
    docker: "~{docker}"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    maxRetries: 3
  }
  output {
    File freyja_variants = "~{samplename}_freyja_variants.tsv"
    File freyja_depths = "~{samplename}_freyja_depths.tsv"
    File freyja_demixed = "~{samplename}_freyja_demixed.tsv"
    File? freyja_update_log = "freyja_update.log"
    File? freyja_boostrap_lineages = "~{samplename}_lineages.csv"
    File? freyja_boostrap_lineages_pdf = "~{samplename}_lineages.pdf"
    File? freyja_boostrap_summary = "~{samplename}_summarized.csv"
    File? freyja_boostrap_summary_pdf = "~{samplename}_summarized.pdf"
    String freyja_barcode_version = read_string("FREYJA_BARCODES")
    String freyja_metadata_version = read_string("FREYJA_METADATA")
    String freyja_version = read_string("FREYJA_VERSION")
  }
}