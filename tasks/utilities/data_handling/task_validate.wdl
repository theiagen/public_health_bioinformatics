version 1.0

task theiavalidate {
  input {
    File datatable1_tsv
    File datatable2_tsv
    String columns_to_compare
    String output_prefix

    File? validation_criteria_tsv
    File? column_translation_tsv
    String? na_values
    Boolean debug_output = false
    
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/theiavalidate:1.1.2"
    Int disk_size = 10
    Int memory = 4
    Int cpu = 2
  }
  String datatable1_name = basename(datatable1_tsv)
  String datatable2_name = basename(datatable2_tsv)
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    # grab theiavalidate version
    theiavalidate.py -v > VERSION

    # run theiavalidate
    theiavalidate.py \
      ~{datatable1_tsv} \
      ~{datatable2_tsv} \
      --columns_to_compare ~{columns_to_compare} \
      --output_prefix ~{output_prefix} \
      ~{"--na_values " + na_values} \
      ~{"--validation_criteria " + validation_criteria_tsv} \
      ~{"--column_translation " + column_translation_tsv} \
      ~{true="--debug" false="--verbose" debug_output}
  >>>
  runtime {
    docker: docker 
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    String theiavalidate_version = read_string("VERSION")
    File summary_pdf_report = "~{output_prefix}_summary.pdf"
    File summary_html_report = "~{output_prefix}_summary.html"
    File filtered_input_table1 = "filtered_~{datatable1_name}"
    File filtered_input_table2 = "filtered_~{datatable2_name}"
    File exact_differences = "~{output_prefix}_exact_differences.tsv"
    File? validation_criteria_differences = "~{output_prefix}_validation_criteria_differences.tsv"
    Array[File]? diffs = glob("*diffs*/*")
  }
}
