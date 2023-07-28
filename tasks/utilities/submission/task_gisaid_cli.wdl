version 1.0

task gisaid_upload {
  input {
    # output from Mercury workflows
    File concatenated_fastas
    File concatenated_metadata
    
    # provided by GISAID
    String client_id

    # authentication token will either be generated or must be provided
    File? gisaid_credentials
    File? authentication_file
    
    String frameshift_notification = "catch_novel" 
    # options: 
    # "catch_all" = Notify me about ALL DETECTED FRAMESHIFTS in this submission for reconfirmation of affected sequences
    # "catch_novel" = Notify me only about NOT PREVIOUSLY REPORTED FRAMESHIFTS in this submission for reconfirmation of affected sequences
    # "catch_none" = I confirm ANY FRAMESHIFTS in this submission and request their release without confirmation by a curator
    
    Int disk_size = 100
  }
  command <<<
    # capture GISAID CLI version
    cli3 version > GISAID_VERSION

    # generate an authentication token if username and password are provided
    if [[ -e "~{gisaid_credentials}" ]]; then
      echo "here"
      username=$(cut -f1 "~{gisaid_credentials}")
      password=$(cut -f2 "~{gisaid_credentials}")
      
      echo "trying to authenticate..."
      cli3 authenticate --username "${username}" \
        --token gisaid.authtoken \
        --password "${password}" \
        --client_id "~{client_id}" \
        --force > submission_log.txt

      cat submission_log.txt

    else # otherwise, authentication token must be provided
      cat ~{authentication_file} > gisaid.authtoken
    fi

    ls

    # upload to GISAID
    cli3 upload \
      --token gisaid.authtoken \
      --metadata ~{concatenated_metadata} \
      --fasta ~{concatenated_fastas} \
      --frameshift "~{frameshift_notification}" \
      --failed "failed.txt" >> submission_log.txt

    # TO-DO: parse out the EPI_ISL accessions for easy addition to Terra table
  >>>
  output {
    String gisaid_cli_version = read_string("GISAID_VERSION")
    File gisaid_logs = "submission_log.txt"
    File? failed_uploads = "failed.txt"
  }
  runtime {
    cpu: 1
    memory: "2 GB"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    docker: "us-docker.pkg.dev/general-theiagen/broadinstitute/gisaid-cli:3.0"
  }
}