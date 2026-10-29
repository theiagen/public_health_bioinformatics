version 1.0

task minimap2 {
  meta {
    description: "Align a query genome to a reference with minimap2"
  }
  input {
    File query1
    File? query2
    File reference
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/minimap2:2.22" # newer versions seem to be bugged (infinite loop)
    String mode = "asm20"
    Boolean output_sam = false
    Boolean long_read_flags = false
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    # Preset options - https://lh3.github.io/minimap2/minimap2.html
    # Version capture
    minimap2 --version | tee VERSION

    if [ -z "~{query2}" ] ; then
      INPUT_QUERY="~{query1}"
    else
      INPUT_QUERY="~{query1} ~{query2}"
    fi

    # Run minimap2 - output can be sam or paf file depending on ~{output_sam}
    minimap2 \
      ~{true="-a" false="" output_sam} \
      ~{true="-L --cs --MD" false="" long_read_flags} \
      -x "~{mode}" \
      -t "~{cpu}" \
      "~{reference}" \
      ${INPUT_QUERY} > "~{samplename}"_minimap2.out

    # count mapped records: SAM records without the unmapped flag (4), or any PAF line
    if ~{output_sam}; then
      mapped_count=$(awk -F'\t' '!/^@/ && int($2 / 4) % 2 == 0 {n++} END {print n + 0}' "~{samplename}"_minimap2.out)
    else
      mapped_count=$(awk 'NF {n++} END {print n + 0}' "~{samplename}"_minimap2.out)
    fi
    echo "DEBUG: ${mapped_count} mapped records"
    if [[ "${mapped_count}" -eq 0 ]]; then
      echo "ERROR: no reads mapped to the provided reference sequence" | tee MINIMAP2_WARNING
    else
      touch MINIMAP2_WARNING
    fi

  >>>
  output {
    File minimap2_out = "~{samplename}_minimap2.out"
    String minimap2_version = read_string("VERSION")
    String minimap2_docker = "~{docker}"
    String minimap2_warning = read_string("MINIMAP2_WARNING")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
