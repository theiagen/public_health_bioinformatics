version 1.0

task freyja_microreact {
  meta {
    description: "Converts a long-format Freyja TSV into a Microreact-ready file; emits an empty file when every sample is below the coverage threshold."
  }
  input {
    File freyja_parsed_format_tsv
    String freyja_plot_name
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/freyja-microreact:1.0.1"
    Int disk_size = 50
    Int memory = 4
    Int cpu = 2
  }
  command <<<
    # if every sample fell below the coverage threshold, the upstream long-format
    # task writes the literal string "all samples are below coverage" — short-circuit
    # to an empty microreact file in that case so the workflow keeps moving
    if grep -q "all samples are below coverage" ~{freyja_parsed_format_tsv}; then
      echo "All samples are below coverage threshold, creating empty microreact file"
      touch ~{freyja_plot_name}.microreact
    else
      freyja_microreact.py ~{freyja_parsed_format_tsv} --output ~{freyja_plot_name}.microreact
    fi
  >>>
  output {
    File freyja_microreact_output = "~{freyja_plot_name}.microreact"
    String freyja_microreact_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
  }
}
