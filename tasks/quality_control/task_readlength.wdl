version 1.0

task readlength {
  input {
    File read1
    File read2
    Int memory = 8
    String docker = "quay.io/staphb/bbtools:38.76"
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE

    readlength.sh in=~{read1} > STDOUT_FORWARD
    readlength.sh in=~{read2} > STDOUT_REVERSE

    avg_forward=$(cat STDOUT_FORWARD | grep "#Avg:" | cut -f 2)
    avg_reverse=$(cat STDOUT_REVERSE | grep "#Avg:" | cut -f 2)

    result=$(awk "BEGIN { printf \"%.2f\", ($avg_forward +  $avg_reverse ) / 2 }")
    echo $result | tee AVERAGE_READ_LENGTH
  >>>
  output {
    Float average_read_length = read_string("AVERAGE_READ_LENGTH")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}
