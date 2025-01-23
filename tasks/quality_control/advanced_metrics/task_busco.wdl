version 1.0

task busco {
  meta {
    description: "Run BUSCO on assemblies"
  }
  input {
    File assembly
    String samplename
    String gcp_bucket = "gs://theiagen-large-public-files-rp/terra/databases/busco"
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/busco:gcp_dev3"
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
    Boolean eukaryote = false
  }
  command <<<
    set -e
    # get version
    busco --version | tee "VERSION"


    project_id=$(gcloud config get-value project)
    echo "Project ID: ${project_id}"

    # Test bucket access (will fail if permissions are incorrect)
    echo "Testing bucket access..."
    gsutil -u ${project_id} ls ~{gcp_bucket} || echo "Warning: Unable to list bucket contents"

    # run busco
    # -i input assembly
    # -m geno for genome input
    # -o output file tag
    # --auto-lineage-euk looks at only eukaryotic organisms
    # --auto-lineage-prok looks at only prokaryotic organisms; default
    busco \
      -i ~{assembly} \
      -c ~{cpu} \
      -m geno \
      -o ~{samplename} \
      --gcp_bucket ~{gcp_bucket} \
      --gcp_project ${project_id} \
      ~{true='--auto-lineage-euk' false='--auto-lineage-prok' eukaryote}

    # check for existence of output file; otherwise display a string that says the output was not created
    if [ -f ~{samplename}/short_summary.specific.*.~{samplename}.txt ]; then

      # grab the database version and format it according to BUSCO recommendations
      # pull line out of final specific summary file
      # cut out the database name and date it was created
      # sed is to remove extra comma and to add parentheses around the date and remove all tabs
      # finally write to a file called DATABASE
      cat ~{samplename}/short_summary.specific.*.~{samplename}.txt | grep "dataset is:" | cut -d' ' -f 6,9 | sed 's/,//; s/ / (/; s/$/)/; s|[\t]||g' | tee DATABASE
      
      # extract the results string; strip off all tab and space characters; write to a file called BUSCO_RESULTS
      cat ~{samplename}/short_summary.specific.*.~{samplename}.txt | grep "C:" | sed 's|[\t]||g; s| ||g' | tee BUSCO_RESULTS

      # rename final output file to predictable name
      cp -v ~{samplename}/short_summary.specific.*.~{samplename}.txt ~{samplename}_busco-summary.txt
    else
      echo "BUSCO FAILED" | tee BUSCO_RESULTS
      echo "NA" > DATABASE
    fi
  >>>
  output {
    String busco_version = read_string("VERSION")
    String busco_docker = docker
    String busco_database = read_string("DATABASE")
    String busco_results = read_string("BUSCO_RESULTS")
    File? busco_report = "~{samplename}_busco-summary.txt"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}