version 1.0

task ncbi_scrub_pe {
  input {
    File read1
    File read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 4
  }
  command <<<
    set -euo pipefail

    # date and version control
    date | tee DATE

    status="PASS"

    # corrupt/truncated gzip
    for file in "~{read1}" "~{read2}"; do
      if ! zcat -f "$file" >/dev/null 2>&1; then
        status="FAIL: corrupt or truncated gzip: $(basename "$file")"
      fi
    done

    if [[ "$status" == "PASS" ]]; then
      r1_lines=$(zcat -f "~{read1}" | wc -l)
      r2_lines=$(zcat -f "~{read2}" | wc -l)
      echo "DEBUG: R1 reads: $((r1_lines / 4)); R2 reads: $((r2_lines / 4))"

      if (( r1_lines == 0 || r2_lines == 0 )); then
        status="FAIL: empty input (R1 lines: $r1_lines, R2 lines: $r2_lines)"
      elif (( r1_lines % 4 != 0 || r2_lines % 4 != 0 )); then
        status="FAIL: truncated FASTQ (R1 lines: $r1_lines, R2 lines: $r2_lines)"
      elif (( r1_lines != r2_lines )); then
        status="FAIL: R1/R2 read counts differ (R1: $((r1_lines / 4)), R2: $((r2_lines / 4)))"
      fi
    fi

    echo "$status" | tee STATUS

    if [[ "$status" != "PASS" ]]; then
      echo "ERROR: $status" >&2
      echo 0 > SPOTS_REMOVED
      # make empty files
      echo -n "" | gzip > ~{samplename}_R1_dehosted.fastq.gz
      echo -n "" | gzip > ~{samplename}_R2_dehosted.fastq.gz
      exit 0
    fi

    # if compressed, unzip read files as scrub tool does not take in .gz fastq files, and interleave them
    # paste command takes 4 lines at a time and merges them into a single line with tabs
    # tr substitutes the tab separators from paste into new lines, effectively interleaving the reads and keeping the FASTQ format
    # Important: To ensure that the reads are interleaved correctly, the reads must be in the same order in both files
    echo "DEBUG: Interleaving reads with paste..."
    paste <(zcat -f "~{read1}" | paste - - - -) <(zcat -f "~{read2}" | paste - - - -) | tr '\t' '\n' > interleaved.fastq
    # dehost reads
    # -x Remove spots instead of default 'N' replacement.
    # -s ; Input is (collated) interleaved paired-end(read) file AND you wish both reads masked or removed.
    echo "DEBUG: Running HRRT..."
    echo "DEBUG: /opt/scrubber/scripts/scrub.sh -p ~{cpu} -x -s -i interleaved.fastq"
    /opt/scrubber/scripts/scrub.sh -p ~{cpu} -x -s -i interleaved.fastq > STDOUT 2> STDERR

    # |& tail -n1 | awk -F" " '{print $1}' > SPOTS_REMOVED ; capture the number of spots removed by fetching the first item on the last line of stderr
    cat STDERR |& tail -n1 | awk -F" " '{print $1}' > SPOTS_REMOVED

    # split interleaved reads and compress files
    # paste takes 8 at a time and merges them into a single line with tabs
    # grouping each pair of reads into one line
    echo "DEBUG: Splitting interleaved dehosted reads..."
    paste - - - - - - - - < interleaved.fastq.clean | cut -f 1-4 | tr '\t' '\n' | gzip > "~{samplename}_R1_dehosted.fastq.gz"
    paste - - - - - - - - < interleaved.fastq.clean | cut -f 5-8 | tr '\t' '\n' | gzip > "~{samplename}_R2_dehosted.fastq.gz"
  >>>
  output {
    File read1_dehosted = "~{samplename}_R1_dehosted.fastq.gz"
    File read2_dehosted = "~{samplename}_R2_dehosted.fastq.gz"
    Int human_spots_removed = read_int("SPOTS_REMOVED")
    String ncbi_scrub_docker = docker
    String ncbi_scrub_status = read_string("STATUS")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task ncbi_scrub_se {
  input {
    File read1
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 4
  }
  command <<<
    # date and version control
    date | tee DATE

    status="PASS"

    # corrupt/truncated gzip
    if ! zcat -f "~{read1}" > read1.fastq 2>/dev/null; then
      status="FAIL: corrupt or truncated gzip: $(basename "~{read1}")"
    fi

    if [[ "$status" == "PASS" ]]; then
      r1_lines=$(wc -l < read1.fastq)
      echo "DEBUG: R1 reads: $((r1_lines / 4))"

      if (( r1_lines == 0 )); then
        status="FAIL: empty input (R1 lines: $r1_lines)"
      elif (( r1_lines % 4 != 0 )); then
        status="FAIL: truncated FASTQ (R1 lines: $r1_lines)"
      fi
    fi

    echo "$status" | tee STATUS

    if [[ "$status" != "PASS" ]]; then
      echo "ERROR: $status" >&2
      echo 0 > SPOTS_REMOVED
      # make empty files
      echo -n "" | gzip > ~{samplename}_R1_dehosted.fastq.gz
      exit 0
    fi

    # dehost reads
    # -x Remove spots instead of default 'N' replacement.
    # |& tail -n1 | awk -F" " '{print $1}' > SPOTS_REMOVED ; capture the number of spots removed by fetching the first item on the last line of stderr
    echo "DEGUB: Running HRRT..."
    echo "DEBUG: /opt/scrubber/scripts/scrub.sh -p ~{cpu} -x -i read1.fastq"
    /opt/scrubber/scripts/scrub.sh -p ~{cpu} -x -i read1.fastq > STDOUT 2> STDERR

    #tail -n1 STDERR | awk -F" " '{print $1}' > SPOTS_REMOVED ; capture the number of spots removed by fetching the first item on the last line of stderr
    tail -n1 STDERR | awk -F" " '{print $1}' > SPOTS_REMOVED

    # gzip dehosted reads
    echo "DEGUB: Renaming dehosted reads..."
    gzip read1.fastq.clean -c > ~{samplename}_R1_dehosted.fastq.gz
  >>>
  output {
    File read1_dehosted = "~{samplename}_R1_dehosted.fastq.gz"
    Int human_spots_removed = read_int("SPOTS_REMOVED")
    String ncbi_scrub_docker = docker
    String ncbi_scrub_status = read_string("STATUS")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 0
  }
}
