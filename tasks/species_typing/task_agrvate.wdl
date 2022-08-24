version 1.0

task agrvate {
  meta {
    description: "Rapid identification of Staphylococcus aureus agr locus type and agr operon variants."
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/agrvate:1.0.2--hdfd78af_0"
    Int? cpu = 1

    # Parameters
    # --typing_only    agr typing only. Skips agr operon extraction and frameshift detection
    Boolean typing_only = false
  }
  command <<<
    echo $(agrvate -v 2>&1) | sed 's/agrvate v//;' | tee VERSION
    agrvate \
        ~{true="--typing_only" false="" typing_only} \
        -i $fasta_name
    cp results/~{samplename}-summary.tab ~{samplename}.tsv
    tar -czvf ~{samplename}.tar.gz results/
  >>>
  output {
    File agrvate_summary = "~{samplename}.tsv"
    File agrvate_results = "~{samplename}.tar.gz"
    String agrvate_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
