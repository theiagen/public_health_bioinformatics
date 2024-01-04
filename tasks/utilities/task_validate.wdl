version 1.0

task export_two_tsvs {
  input {
    String terra_project1
    String? terra_project2
    String terra_workspace1
    String? terra_workspace2
    String datatable1
    String datatable2
    Int disk_size = 10
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    python3 /scripts/export_large_tsv/export_large_tsv.py --project ~{terra_project1} --workspace ~{terra_workspace1} --entity_type ~{datatable1} --tsv_filename "~{datatable1}.tsv"

    # check if second project is provided; if not, use first
    if [[ -z "~{terra_project2}" ]]; then
      PROJECT2="~{terra_project1}"
    else
      PROJECT2="~{terra_project2}"
    fi

    # check if second workspace is provided; if not, use first
    if [[ -z "~{terra_workspace2}" ]]; then
      WORKSPACE2="~{terra_workspace1}"
    else
      WORKSPACE2="~{terra_workspace2}"
    fi

    python3 /scripts/export_large_tsv/export_large_tsv.py --project ${PROJECT2} --workspace ${WORKSPACE2} --entity_type ~{datatable2} --tsv_filename "~{datatable2}.tsv"

    if [[ $(wc -l ~{datatable1} | cut -f1 -d' ') -eq $(wc -l ~{datatable2} | cut -f1 -d' ') ]]; then
      echo true | tee CONTINUE
    else 
      echo false | tee CONTINUE
    fi
  >>>
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "1 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File datatable1_tsv = "~{datatable1}.tsv"
    File datatable2_tsv = "~{datatable2}.tsv"
    Boolean same_table_length = read_boolean("CONTINUE")
  }
}

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

    Int disk_size = 10
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
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/theiavalidate:0.0.1"
    memory: "4 GB"
    cpu: 2
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
  }
}