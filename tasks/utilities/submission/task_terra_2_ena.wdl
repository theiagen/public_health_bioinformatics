task terra_reads_to_ena {
  input {
    Array[String] sample_names
    File input_table
    # Column name configuration for read columns
    String read1_column_name = "read1"
    String read2_column_name = "read2"
    String bam_column_name = "bam_file"
    String cram_column_name = "cram_file"
    File? column_mappings_file 
    # ENA submission parameters
    String study_accession
    String webin_username
    String webin_password
    String? center_name
    # Read metadata parameters - can be provided as overrides if necessary
    String? library_strategy
    String? library_source
    String? library_selection
    String? platform
    String? instrument
    # Submission parameters
    Boolean test_submission = false
    Int parallel_cores = 1
    Boolean allow_missing_metadata = false
    Boolean force_submission = false
    # Resource parameters
    Int memory_gb = 16
    Int cpu = 8
    Int disk_size_gb = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra_to_ena:0.4"
  }
  
  command <<<
    set -euo pipefail
    
    mkdir -p data
    mkdir -p output
    mkdir -p logs
    
    # Copy the input table to the working directory, this will come from the sample registry
    cp "~{input_table}" input_table.tsv
    INPUT_TABLE="input_table.tsv"
    
    # Convert sample names array to comma-separated list if provided
    if [[ -n "~{sep=',' sample_names}" ]]; then
        SAMPLE_LIST="--samples ~{sep=',' sample_names}"
    else
        SAMPLE_LIST=""
    fi
    
    echo "Step 1: Preparing metadata spreadsheet..."
    
    # Set up column mappings parameter if provided, make sure the file is copied to the working directory
    COLUMN_MAPPINGS_PARAM=""
    if [[ -n "~{column_mappings_file}" ]]; then
      cp "~{column_mappings_file}" column_mappings.tsv
      COLUMN_MAPPINGS_PARAM="--column-mappings column_mappings.tsv"
    fi
    
    # Run the metadata preparation script
    echo "Preparing ENA metadata..."
    prepare_ena_data.py \
      --input "$INPUT_TABLE" \
      --output "output/ena_reads_spreadsheet.txt" \
      --excluded "output/excluded_samples.tsv" \
      --file-paths "output/file_paths.json" \
      --sample-id-column "~{table_name}_id" \
      $SAMPLE_LIST \
      $COLUMN_MAPPINGS_PARAM \
      ~{true="--allow-missing" false="" allow_missing_metadata} \
      --read1-column "~{read1_column_name}" \
      --read2-column "~{read2_column_name}" \
      --bam-column "~{bam_column_name}" \
      --cram-column "~{cram_column_name}" \
      --study-accession "~{study_accession}" \
      ~{"--library-strategy" + library_strategy} \
      ~{"--library-source" + library_source} \
      ~{"--library-selection" + library_selection} \
      ~{"--platform" + platform} \
      ~{"--instrument" + instrument}
    
    echo "Step 2: Localizing data files..."
    
    # Run the file localization script
    echo "Localizing files..."
    localize_files --file-paths "output/file_paths.json"
    
    echo "Step 3: Validating submission with ENA..."
    
    # Setting these up for now and making them more explicit
    TEST_FLAG=""
    if [ ~{test_submission} == true ]; then
      TEST_FLAG="-t"
    fi
    
    CENTER_FLAG=""
    if [ -n "~{center_name}" ]; then
      CENTER_FLAG="-c \"~{center_name}\""
    fi
    
    # Run validation first
    echo "Starting ENA validation..."
    bulk_webincli.py \
      -u ~{webin_username} \
      -p ~{webin_password} \
      -g reads \
      -w /scripts/webin-cli.jar \
      -s output/ena_reads_spreadsheet.txt \
      -d . \
      -m validate \
      $TEST_FLAG \
      -pc ~{parallel_cores} \
      $CENTER_FLAG > logs/validation.log 2>&1
    
    # Check if validation failed
    VALIDATION_FAILED=0
    
    # Check for validation reports indicating failures
    VALIDATION_REPORTS=$(find . -path "*/validate/*.report" 2>/dev/null)
    if [ -n "$VALIDATION_REPORTS" ]; then
      echo "Found validation report files. Checking for errors..."
      for report in $VALIDATION_REPORTS; do
        if grep -q "ERROR" "$report"; then
          echo "Validation errors found in $report"
          VALIDATION_FAILED=1
          mkdir -p output/validation_reports
          cp "$report" "output/validation_reports/$(basename $report)"
        fi
      done
    fi
    
    # Check for webin-cli.report files
    WEBIN_REPORTS=$(find . -path "*/submissions/webin-cli.report" 2>/dev/null)
    if [ -n "$WEBIN_REPORTS" ]; then
      echo "Found webin-cli.report files. Checking for errors..."
      for report in $WEBIN_REPORTS; do
        if grep -q "ERROR" "$report"; then
          echo "Validation errors found in $report"
          VALIDATION_FAILED=1
        fi
        mkdir -p output/validation_reports
        cp "$report" "output/validation_reports/$(basename $(dirname $(dirname $report)))_webin-cli.report"
      done
    fi
    
    # Also check log output for errors
    if grep -q "VALIDATION FAILED" logs/validation.log || grep -q "ERROR" logs/validation.log; then
      echo "Validation failed according to log output."
      VALIDATION_FAILED=1
    fi
    
    if [ $VALIDATION_FAILED -eq 0 ]; then
      echo "Validation succeeded."
    else
      echo "Validation failed. Check the validation reports for details."
    fi
    
    # Create validation report
    echo "=== ENA Validation Report ===" > output/validation_report.txt
    echo "Timestamp: $(date)" >> output/validation_report.txt
    echo "Test mode: ~{test_submission}" >> output/validation_report.txt
    
    if [ $VALIDATION_FAILED -eq 1 ]; then
      echo "Status: FAILED" >> output/validation_report.txt
      echo "See logs/validation.log and failed_validation.txt for details" >> output/validation_report.txt
      cat failed_validation.txt >> output/validation_report.txt 2>/dev/null || echo "No detailed validation errors found" >> output/validation_report.txt
    else
      echo "Status: PASSED" >> output/validation_report.txt
    fi
    
    echo "Step 4: Submitting to ENA..."
    
    # Run the actual submission
    run_bulk_webincli \
      -u ~{webin_username} \
      -p ~{webin_password} \
      -g reads \
      -w /scripts/webin-cli.jar \
      -s output/ena_reads_spreadsheet.txt \
      -d . \
      -m submit \
      -pc ~{parallel_cores} \
      $TEST_FLAG \
      $CENTER_FLAG > logs/submission.log 2>&1
    
    # Create a submission report
    echo "=== ENA Submission Report ===" > output/submission_report.txt
    echo "Timestamp: $(date)" >> output/submission_report.txt
    echo "Test submission: ~{test_submission}" >> output/submission_report.txt
    echo "Parallel cores: ~{parallel_cores}" >> output/submission_report.txt
    echo "" >> output/submission_report.txt
    
    # Find receipt XML files based on Webin CLI output structure
    mkdir -p output/receipt_xmls
    find . -path "*/submit/receipt.xml" > receipt_xmls.txt || echo "No receipt XMLs found"
    
    
    # Collect all logs and reports
    cp failed_validation.txt output/failed_validation.txt 2>/dev/null || touch output/failed_validation.txt
    cp logs/validation.log output/validation.log
    cp logs/submission.log output/submission.log
  >>>
  
  output {
    # Metadata files
    File metadata_spreadsheet = "output/ena_reads_spreadsheet.txt"
    File excluded_samples = "output/excluded_samples.tsv"
    File file_paths_mapping = "output/file_paths.json"
    Array[File] manifests = glob("output/manifests/*")
  }
  
  runtime {
    docker: docker
    memory: "${memory_gb} GB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} SSD"
    preemptible: 0
  }
}