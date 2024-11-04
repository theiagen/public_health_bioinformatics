version 1.0

task get_srr_from_biosamples {
  input {
    String biosample_accession 
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0"
    Int disk_size = 10
    Int cpu = 2
    Int memory = 8
  }

  meta {
    volatile: true
  }

  command <<< 
    mkdir -p metadata_output
    date -u | tee DATE

    # Debug output to show the biosample being processed
    echo "Fetching metadata for BioSample: ${biosample_accession}"

    # Use fastq-dl to fetch metadata only
    fastq-dl --accession ~{biosample_accession} --provider SRA --outdir metadata_output --only-download-metadata --verbose

    if [[ -f metadata_output/fastq-run-info.tsv ]]; then
      echo "Metadata written for ${biosample_accession}:"
      echo "TSV content:"
      cat metadata_output/fastq-run-info.tsv

      # Extract the SRR accession (assuming it's in the first column)
      SRR_accessions=$(awk -F'\t' 'NR>1 {print $1}' metadata_output/fastq-run-info.tsv)
      echo "Extracted SRR accessions: ${SRR_accessions}"

      # Output the SRR accessions as a single string
      echo "${SRR_accessions}" > metadata_output/srr_accession.txt
    else
      echo "No metadata found for ${biosample_accession}"
      exit 1
    fi
  >>>

  output {
    String srr_accession = read_string("metadata_output/srr_accession.txt")
  }

  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disk: disk_size + " GB"
    preemptible: 1
  }
}
