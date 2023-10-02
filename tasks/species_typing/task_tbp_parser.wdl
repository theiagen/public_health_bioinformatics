version 1.0

task tbp_parser {
  input {
    File tbprofiler_json
    File tbprofiler_bam
    File tbprofiler_bai
    String samplename

    String? sequencing_method
    String? operator
    Int min_depth = 10
    Int coverage_threshold = 100
    Boolean tbp_parser_debug = false

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:1.0.1"
    Int disk_size = 100
    Int memory = 4
    Int cpu = 1
  }
  command <<<
    # get version
    python3 /tbp-parser/tbp_parser/tbp_parser.py --version | tee VERSION

    # run tbp-parser
    python3 /tbp-parser/tbp_parser/tbp_parser.py ~{tbprofiler_json} ~{tbprofiler_bam} \
      ~{"--sequencing_method " + sequencing_method} \
      ~{"--operator " + operator} \
      ~{"--min_depth " + min_depth} \
      ~{"--coverage_threshold " + coverage_threshold} \
      --output_prefix ~{samplename} \
      ~{true="--debug" false="--verbose" tbp_parser_debug}

    # get genome percent coverage for the entire reference genome length over min_depth
    genome=$(samtools depth -J ~{tbprofiler_bam} | awk -F "\t" '{if ($3 >= ~{min_depth}) print;}' | wc -l )
    python3 -c "print ( ($genome / 4411532 ) * 100 )" | tee GENOME_PC
  >>>
  output {
    File tbp_parser_looker_report_csv = "~{samplename}.looker_report.csv"
    File tbp_parser_laboratorian_report_csv = "~{samplename}.laboratorian_report.csv"
    File tbp_parser_lims_report_csv = "~{samplename}.lims_report.csv"
    File tbp_parser_coverage_report = "~{samplename}.percent_gene_coverage.csv"
    Float tbp_parser_genome_percent_coverage = read_float("GENOME_PC")
    String tbp_parser_version = read_string("VERSION")
    String tbp_parser_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 1
  }
}