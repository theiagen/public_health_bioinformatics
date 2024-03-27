version 1.0

task stxtyper {
  input {
    File assembly
    String samplename
    String docker = "kapsakcj/stxtyper:8328c4d"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    # capture date
    date | tee DATE

    # capture version info
    stxtyper --version | tee VERSION.txt

    echo "DEBUG: running StxTyper now..."
    # run StxTyper on assembly; may need to add/remove options in the future if they change
    stxtyper \
      --nucleotide ~{assembly} \
      --name ~{samplename} \
      --output ~{samplename}_stxtyper.tsv

    # parse output TSV
    echo "DEBUG: Parsing StxTyper output TSV..."
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 2 >stxtyper_target_contig.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 3 >stxtyper_stx_type.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 4 >stxtyper_stx_operon_status.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 5 >stxtyper_combined_identity.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 6 >stxtyper_target_start.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 7 >stxtyper_target_stop.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 8 >stxtyper_target_strand.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 9 >stxtyper_a_reference.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 10 >stxtyper_a_identity.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 11 >stxtyper_a_coverage.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 12 >stxtyper_b_reference.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 13 >stxtyper_b_identity.txt
    tail -n 1 ~{samplename}_stxtyper.tsv | cut -f 14 >stxtyper_b_coverage.txt

    echo "DEBUG: Finished parsing StxTyper output TSV."
  >>>
  output {
    File stxtyper_report = "~{samplename}_stxtyper.tsv"
    String stxtyper_docker = docker
    String stxtyper_version = read_string("VERSION.txt")
    String stxtyper_target_contig = read_string("stxtyper_target_contig.txt")
    String stxtyper_stx_type = read_string("stxtyper_stx_type.txt")
    String stxtyper_stx_operon_status = read_string("stxtyper_stx_operon_status.txt")
    String stxtyper_combined_identity = read_string("stxtyper_combined_identity.txt")
    String stxtyper_target_start = read_string("stxtyper_target_start.txt")
    String stxtyper_target_stop = read_string("stxtyper_target_stop.txt")
    String stxtyper_a_reference = read_string("stxtyper_a_reference.txt")
    String stxtyper_a_identity = read_string("stxtyper_a_identity.txt")
    String stxtyper_a_coverage = read_string("stxtyper_a_coverage.txt")
    String stxtyper_b_reference = read_string("stxtyper_b_reference.txt")
    String stxtyper_b_identity = read_string("stxtyper_b_identity.txt")
    String stxtyper_b_coverage = read_string("stxtyper_b_coverage.txt")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
