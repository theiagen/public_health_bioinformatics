version 1.0

task fetch_srr_accession {
  input {
    String sample_accession 
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

    # Debug output to show the sample being processed
    echo "Fetching metadata for sample accession: ${sample_accession}"

    # Use fastq-dl to fetch metadata only
    fastq-dl --accession ~{sample_accession} --outdir metadata_output --only-download-metadata --verbose


    if [[ -f metadata_output/fastq-run-info.tsv ]]; then
      echo "Metadata written for ${sample_accession}:"
      echo "TSV content:"
      cat metadata_output/fastq-run-info.tsv

      # Extract the SRR accession (It is typically in the first column)
      SRR_accessions=$(awk -F'\t' 'NR>1 {print $1}' metadata_output/fastq-run-info.tsv)
      if [[ -z "${SRR_accessions}" ]]; then
        echo "No SRR accession found for ${sample_accession}" > metadata_output/srr_accession.txt
      else
        echo "Extracted SRR accessions: ${SRR_accessions}"
        echo "${SRR_accessions}" > metadata_output/srr_accession.txt
      fi
    else
      echo "No metadata found for ${sample_accession}"
      echo "No SRR accession found" > metadata_output/srr_accession.txt
    fi
  >>>
  output {
    String srr_accession = read_string("metadata_output/srr_accession.txt")
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}
