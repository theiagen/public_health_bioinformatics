version 1.0

task pasty {
  input {
    File  assembly
    String  samplename
    Int min_pident = 95
    Int min_coverage = 95
    String docker = "staphb/pasty:1.0.2"
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE
    pasty --version > VERSION && sed -i -e 's/pasty\, version //' VERSION  
    pasty \
    --assembly ~{assembly} \
    --min_pident ~{min_pident} \
    --min_coverage ~{min_coverage} \
    --prefix ~{samplename} \
    --outdir .  
    awk 'FNR==2' "~{samplename}.tsv" | cut -d$'\t' -f2 > SEROGROUP
    awk 'FNR==2' "~{samplename}.tsv" | cut -d$'\t' -f3 > COVERAGE
    awk 'FNR==2' "~{samplename}.tsv" | cut -d$'\t' -f4 > FRAGMENTS
    awk 'FNR==2' "~{samplename}.tsv" | cut -d$'\t' -f5 > COMMENT
  >>>
  output {
    String pasty_serogroup = read_string("SEROGROUP")
    Float pasty_serogroup_coverage = read_float("COVERAGE")
    Int pasty_serogroup_fragments = read_int("FRAGMENTS")
    File pasty_summary_tsv = "~{samplename}.tsv"
    File pasty_blast_hits = "~{samplename}.blastn.tsv"
    File pasty_all_serogroups = "~{samplename}.details.tsv"
    String pasty_version = read_string("VERSION")
    String pasty_pipeline_date = read_string("DATE")
    String pasty_docker = docker
    String pasty_comment = read_string("COMMENT")
  }
  runtime {
    docker: "~{docker}"
    memory: "4 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}