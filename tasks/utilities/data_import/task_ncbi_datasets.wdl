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
      if [ ! -s ncbi_dataset/data/genomic.fna ]; then
        echo "ERROR: no assemblies found for accession: ~{ncbi_accession}"
        echo "FAIL" > NCBIDATASETS_STATUS
      else
        echo "PASS" > NCBIDATASETS_STATUS
        cp -v ncbi_dataset/data/genomic.fna ./~{ncbi_accession}.fasta
        cp -v ncbi_dataset/data/data_report.jsonl ./~{ncbi_accession}.data_report.jsonl

        # acquire the taxon id for the accession
        datasets summary virus genome accession \
          ~{ncbi_accession} --as-json-lines | \
        dataformat tsv virus-genome --fields virus-name,virus-tax-id | \
        tail -n+2 > accession_taxonomy.tsv

        cut -f 1 accession_taxonomy.tsv > TAXON_NAME
        cut -f 2 accession_taxonomy.tsv > TAXON_ID
      fi
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
      if [ ! -s ncbi_dataset/data/~{ncbi_accession}*/~{ncbi_accession}*.fna ]; then
        echo "ERROR: no assemblies found for accession: ~{ncbi_accession}"
        echo "FAIL" > NCBIDATASETS_STATUS
      else
        echo "PASS" > NCBIDATASETS_STATUS

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

        # acquire the taxon id for the accession
        datasets summary genome accession \
          ~{ncbi_accession} --as-json-lines | \
        dataformat tsv genome --fields organism-name,organism-tax-id | \
        tail -n+2 > accession_taxonomy.tsv

        cut -f 1 accession_taxonomy.tsv > TAXON_NAME
        cut -f 2 accession_taxonomy.tsv > TAXON_ID
      fi
    fi
  >>>
  output {
    File? ncbi_datasets_assembly_fasta = "~{ncbi_accession}.fasta"
    String ncbi_datasets_status = read_string("NCBIDATASETS_STATUS")
    File? ncbi_datasets_gff3 = "~{ncbi_accession}.gff"
    File? ncbi_datasets_gbff = "~{ncbi_accession}.gbff"
    File? ncbi_datasets_assembly_data_report_json = "~{ncbi_accession}.data_report.jsonl"
    String? taxon_name = read_string("TAXON_NAME")
    String? taxon_id = read_string("TAXON_ID")
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

task ncbi_datasets_download_genome_taxon {
  input {
    String taxon
    Int cpu = 1
    Int memory = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:16.38.1" # not the latest version, but it's hard to keep up w/ the frequent releases
    Int disk_size = 50
    Boolean refseq = false
    Boolean complete = false
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
      datasets summary \
        virus \
        genome taxon ~{'"' + taxon + '"'} \
        --limit 1 \
        ~{true="--complete-only" false="" complete} \
        ~{true="--refseq" false="" refseq} \
        --as-json-lines | \
      dataformat tsv virus-genome \
        --fields accession
         > ncbi_genome_summary.tsv

      # report if the file is empty
      if [ ! -s "ncbi_genome_summary.tsv" ]; then
        echo "ERROR: no assemblies found for taxon: ~{taxon}"
        echo "FAIL: no assemblies found" > NCBIDATASETS_STATUS
        echo "N/A" > NCBI_ACCESSION
      else
        ncbi_accession=$(tail -1 ncbi_genome_summary.tsv | cut -f 1)
        echo "PASS" > NCBIDATASETS_STATUS
        echo "DEBUG: ncbi_accession: ${ncbi_accession}"
        echo $ncbi_accession > NCBI_ACCESSION

        # FYI the --include flag options is different for the virus sub-command. Does not include gbff or gff3 files
        datasets download virus genome accession \
          ${ncbi_accession} \
          --filename ncbi_taxon.zip \
          --include genome

        # unzip the archive and copy FASTA and JSON to PWD, rename in the process so output filenames are predictable
        unzip ncbi_taxon.zip
        cp -v ncbi_dataset/data/genomic.fna ./ncbi_taxon.fasta
        cp -v ncbi_dataset/data/data_report.jsonl ./ncbi_taxon.data_report.jsonl
      fi
    # otherwise, use the datasets download' sub-command
    else
      datasets summary \
        genome taxon ~{'"' + taxon + '"'} \
        --limit 1 \
        ~{true="--assembly-level complete" false="" complete} \
        ~{true="--assembly-source RefSeq" false="" refseq} \
        --as-json-lines | \
      dataformat tsv genome \
        --fields accession \
         > ncbi_genome_summary.tsv

      # report if the file is empty
      if [ ! -s "ncbi_genome_summary.tsv" ]; then
        echo "ERROR: no assemblies found for taxon: ~{taxon}"
        echo "FAIL: no assemblies found" > NCBIDATASETS_STATUS
        echo "N/A" > NCBI_ACCESSION
      else
        ncbi_accession=$(tail -1 ncbi_genome_summary.tsv | cut -f 1)
        echo "DEBUG: ncbi_accession: ${ncbi_accession}"
        echo $ncbi_accession > NCBI_ACCESSION
        echo "PASS" > NCBIDATASETS_STATUS
        #### download FASTA file using ncbi_accession ####
        datasets download genome accession \
          ${ncbi_accession} \
          --filename ncbi_taxon.zip \
          --include genome \
          ~{true="--include gbff" false="" include_gbff} \
          ~{true="--include gff3" false="" include_gff3}

        # unzip the archive and copy FASTA and JSON to PWD, rename in the process so output filenames are predictable 
        unzip ncbi_taxon.zip
        cp -v ncbi_dataset/data/${ncbi_accession}*/${ncbi_accession}*.fna ./ncbi_taxon.fasta
        cp -v ncbi_dataset/data/assembly_data_report.jsonl ./ncbi_taxon.data_report.jsonl

        # if GFF3 file exists, rename for output as a file
        if [ $(find . -maxdepth 4 -type f -iname "*.gff" | wc -l) -gt 0 ]; then
          echo ".gff file found, renaming output gff file ..."
          mv -v ncbi_dataset/data/${ncbi_accession}*/*.gff ncbi_taxon.gff
        fi

        # if GBFF file exists, rename for output as a file
        if [ $(find . -maxdepth 4 -type f -iname "*.gbff" | wc -l) -gt 0 ]; then
          echo ".gbff file found, renaming output gbff file ..."
          mv -v ncbi_dataset/data/${ncbi_accession}*/*.gbff ncbi_taxon.gbff
        fi
      fi
    fi
  >>>
  output {
    File? ncbi_datasets_assembly_fasta = "ncbi_taxon.fasta"
    String ncbi_datasets_accession = read_string("NCBI_ACCESSION")
    String ncbi_datasets_status = read_string("NCBIDATASETS_STATUS")
    File? ncbi_datasets_gff3 = "ncbi_taxon.gff"
    File? ncbi_datasets_gbff = "ncbi_taxon.gbff"
    File? ncbi_datasets_assembly_data_report_json = "ncbi_taxon.data_report.jsonl"
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

task ncbi_datasets_genome_summary {
  input {
    String taxon # can be taxon id (int) or organism name (string)
    Boolean use_ncbi_virus = false
    Int? summary_limit 
    Boolean complete = true
    Boolean refseq = true
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:16.38.1" # not the latest version, but it's hard to keep up w/ the frequent releases
    Int cpu = 1
    Int memory = 4
    Int disk_size = 50
  }
  command <<<
    set -euo pipefail

    date | tee DATE
    datasets --version | sed 's|datasets version: ||' | tee DATASETS_VERSION

    if ~{use_ncbi_virus}; then
      datasets summary \
        virus \
        genome taxon ~{'"' + taxon + '"'} \
        ~{"--limit " + summary_limit} \
        ~{true="--complete-only" false="" complete} \
        ~{true="--refseq" false="" refseq} \
        --as-json-lines | \
      dataformat tsv virus-genome \
        --fields accession,completeness,is-annotated,length,sourcedb,virus-name,virus-tax-id \
         > ncbi_genome_summary.tsv
    else
      datasets summary \
        genome taxon ~{'"' + taxon + '"'} \
        ~{"--limit " + summary_limit} \
        ~{true="--assembly-level complete" false="" complete} \
        ~{true="--assembly-source RefSeq" false="" refseq} \
        --as-json-lines | \
      dataformat tsv genome \
        --fields accession,completeness,is-annotated,length,sourcedb,oragnism-name,organism-tax-id \
         > ncbi_genome_summary.tsv
    fi 


    # take the average of of all genome lengths across each row
    awk -F'\t' '{
      if (NR > 1) { sum += $4; count++ }
    }
    END {
      if (count) { print int(sum / count) } else { print "No matching taxon id found" }
    }' ncbi_genome_summary.tsv > AVG_GENOME_LENGTH

  >>>
  output {
    File taxon_summary_tsv = "ncbi_genome_summary.tsv"
    Int avg_genome_length = read_string("AVG_GENOME_LENGTH")
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