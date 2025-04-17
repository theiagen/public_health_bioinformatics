version 1.0

task polypolish {
  input {
    String samplename
    File assembly_fasta
    File read1_sam # these files need to be aligned to the draft assembly with `-a` flag
    File read2_sam # see also the task_bwa.wdl#bwa_all task

    Int illumina_polishing_rounds = 1 # Default: 1 round of polishing
    String? pair_orientation # default: auto
    Float? low_percentile_threshold # default: 0.1
    Float? high_percentile_threshold # default: 99.9

    Float? fraction_invalid # default: 0.2
    Float? fraction_valid # default 0.5
    Int? maximum_errors # default: 10
    Int? minimum_depth # default: 5
    Boolean careful = false # ignore any reads with multiple alignments

    Int cpu = 1 # polypolish is single-threaded
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/polypolish:0.6.0"
    Int memory = 8
  }
  command <<<
    set -euo pipefail

    polypolish --version | tee VERSION

    # Initial input for polishing
    polished_assembly="~{assembly_fasta}"

    for i in $(seq 1 ~{illumina_polishing_rounds}); do
      echo "Starting Polypolish round $i..."

      # Filter SAM files
      polypolish filter \
        --in1 ~{read1_sam} \
        --in2 ~{read2_sam} \
        --out1 ${polished_assembly}_filtered1.sam \
        --out2 ${polished_assembly}_filtered2.sam \
        ~{"--orientation " + pair_orientation} \
        ~{"--low " + low_percentile_threshold} \
        ~{"--high " + high_percentile_threshold}

      # Perform polishing
      polypolish polish ${polished_assembly} \
        ${polished_assembly}_filtered1.sam \
        ${polished_assembly}_filtered2.sam \
        ~{"--fraction_invalid " + fraction_invalid} \
        ~{"--fraction_valid " + fraction_valid} \
        ~{"--max_errors " + maximum_errors} \
        ~{"--min_depth " + minimum_depth} \
        ~{true="--careful" false="" careful} \
        > "${polished_assembly}_round${i}.fasta"

      polished_assembly="${polished_assembly}_round${i}.fasta"
    done

    # Final polished output
    mv "${polished_assembly}" "~{samplename}_final_polished.fasta"
  >>>
  output {
    String polypolish_version = read_string("VERSION")
    File polished_assembly = "~{samplename}_final_polished.fasta"
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}