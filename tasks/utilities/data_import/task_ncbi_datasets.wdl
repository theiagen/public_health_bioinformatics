version 1.0

task ncbi_datasets_download_genome_accession {
  input {
    String ncbi_accession
    Int cpu = 1
    Int memory = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:16.38.1" # not the latest version, but it's hard to keep up w/ the frequent releases
    Int disk_size = 50
    Boolean include_gbff = false
    Boolean include_gff3 = false
    Boolean use_ncbi_virus = false
  }
  meta {
    # added so that call caching is always turned off 
    volatile: true
  }
  command <<<
    date | tee DATE
    datasets --version | sed 's|datasets version: ||' | tee DATASETS_VERSION

    # if use_ncbi_virus is true, then use 'datasets download virus' sub-command
    if ~{use_ncbi_virus}; then
      # FYI the --include flag options is different for the virus sub-command. Does not include gbff or gff3 files
      datasets download virus genome accession \
        ~{ncbi_accession} \
        --filename ~{ncbi_accession}.zip \
        --include genome 

      # unfortunately have duplicate code blocks for renaming files because of the differences in 
      # output directory structure and output filenames between virus and non-virus downloads
      # :( 

      # unzip the archive and copy FASTA and JSON to PWD, rename in the process so output filenames are predictable 
      unzip ~{ncbi_accession}.zip
      cp -v ncbi_dataset/data/genomic.fna ./~{ncbi_accession}.fasta
      cp -v ncbi_dataset/data/data_report.jsonl ./~{ncbi_accession}.data_report.jsonl

    # otherwise, use the datasets download' sub-command
    else
      #### download FASTA file using ncbi_accession ####
      # '--assembly-version latest' ensures the most recent version is downloaded, not previous versions
      # NOTE: I have noticed that occasionally the same command may fail one moment with the error below and succeed the second time. usually.
      ### "Error: No assemblies found that match selection"
      # always include genome fna/fasta file
      # only include gbff or gff3 files if user sets booleans to true
      # FYI: currently removing this option "--assembly-version latest" since it's not yet available for downloading from ncbi virus
      # may be necessary to add back in later (for non-viruses) to ensure ONLY the latest version is downloaded and not all versions. Not clear to me what the default is.
      datasets download genome accession \
        ~{ncbi_accession} \
        --filename ~{ncbi_accession}.zip \
        --include genome \
        ~{true="--include gbff" false="" include_gbff} \
        ~{true="--include gff3" false="" include_gff3}

      # unzip the archive and copy FASTA and JSON to PWD, rename in the process so output filenames are predictable 
      unzip ~{ncbi_accession}.zip
      cp -v ncbi_dataset/data/~{ncbi_accession}*/~{ncbi_accession}*.fna ./~{ncbi_accession}.fasta
      cp -v ncbi_dataset/data/assembly_data_report.jsonl ./~{ncbi_accession}.data_report.jsonl

      # if GFF3 file exists, rename for output as a file
      if [ $(find . -maxdepth 4 -type f -iname "*.gff" | wc -l) -gt 0 ]; then
        echo ".gff file found, renaming output gff file ..."
        mv -v ncbi_dataset/data/~{ncbi_accession}*/*.gff ~{ncbi_accession}.gff
      fi

      # if GBFF file exists, rename for output as a file
      if [ $(find . -maxdepth 4 -type f -iname "*.gbff" | wc -l) -gt 0 ]; then
        echo ".gbff file found, renaming output gbff file ..."
        mv -v ncbi_dataset/data/~{ncbi_accession}*/*.gbff ~{ncbi_accession}.gbff
      fi
    fi
  >>>
  output {
    File ncbi_datasets_assembly_fasta = "~{ncbi_accession}.fasta"
    File? ncbi_datasets_gff3 = "~{ncbi_accession}.gff"
    File? ncbi_datasets_gbff = "~{ncbi_accession}.gbff"
    File ncbi_datasets_assembly_data_report_json = "~{ncbi_accession}.data_report.jsonl"
    String ncbi_datasets_version = read_string("DATASETS_VERSION")
    String ncbi_datasets_docker = docker
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}

task ncbi_datasets_viral_taxon_summary {
  input {
    String taxon_id # ideally this is SPECIES lvl taxon id
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:16.38.1" # not the latest version, but it's hard to keep up w/ the frequent releases
    Int cpu = 1
    Int memory = 4
    Int disk_size = 50
  }
  command <<<
    set -euo pipefail

    date | tee DATE
    datasets --version | sed 's|datasets version: ||' | tee DATASETS_VERSION

    datasets summary virus genome taxon ~{taxon_id} \
      --limit 100 \
      --complete-only \
      --refseq \
      --as-json-lines | \

    dataformat tsv virus-genome \
      --fields accession,completeness,is-annotated,length,sourcedb,virus-name,virus-tax-id \
       > ~{taxon_id}.ncbi_viral_genome_summary.tsv

    # take the average of of all genome lengths across each row
    awk -F'\t' '{
      if (NR > 1) { sum += $4; count++ }
    }
    END {
      if (count) { print int(sum / count) } else { print "No matching taxon id found" }
    }' ~{taxon_id}.ncbi_viral_genome_summary.tsv > ~{taxon_id}.AVG_GENOME_LENGTH

  >>>
  output {
    File taxon_summary_tsv = "~{taxon_id}.ncbi_viral_genome_summary.tsv"
    Int avg_genome_length = read_string("~{taxon_id}.AVG_GENOME_LENGTH")
    String ncbi_datasets_version = read_string("DATASETS_VERSION")
    String ncbi_datasets_docker = docker
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 1
  }
}

task ncbi_datasets_identify {
  input {
    String taxon # can be taxon id (int) or organism name (string)
    String? rank # limit input taxon to the user specified rank
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:16.38.1"
    Int cpu = 1
    Int memory = 4
    Int disk_size = 50
  }
  command <<<
    set -euo pipefail

    date | tee DATE
    datasets --version | sed 's|datasets version: ||' | tee DATASETS_VERSION

    datasets summary taxonomy taxon ~{taxon} \
      ~{"--rank " + rank} \
      --as-json-lines | \

    dataformat tsv taxonomy \
      --template tax-summary > ncbi_taxon_summary.tsv

    # extract the taxid ($2) tax name ($3) and rank ($5) from the output tsv file
    # skip the header line. if the input taxon rank is invalid (too specific), there will be more than one line listed here.
    awk -F'\t' 'NR == 2 {print $2}' ncbi_taxon_summary.tsv > TAXON_ID
    awk -F'\t' 'NR == 2 {print $3}' ncbi_taxon_summary.tsv > TAXON_NAME
    awk -F'\t' 'NR == 2 {print $5}' ncbi_taxon_summary.tsv > TAXON_RANK

    # check if ncbi_taxon_summary.tsv is longer than 2 lines (header + 1 data line)
    # if so, then the taxon rank (either calculated or provided) is too specific and we need to exit with an error
    if [ $(wc -l < ncbi_taxon_summary.tsv) -gt 2 ]; then
      echo echo "Error: input taxon rank '~{rank}' is not valid (too specific) for taxon: '~{taxon}'." 1>&2
      exit 1
    fi
  >>>
  output {
    String taxon_id = read_string("TAXON_ID")
    String taxon_name = read_string("TAXON_NAME")
    String taxon_rank = read_string("TAXON_RANK")
    String ncbi_datasets_version = read_string("DATASETS_VERSION")
    String ncbi_datasets_docker = docker
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}