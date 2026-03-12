version 1.0

task tbp_parser {
  input {
    String samplename
    # required/positional arguments
    File tbprofiler_json
    File tbprofiler_bam
    File tbprofiler_bai
    # file arguments
    File? config
    File? coverage_bed
    File? err_coverage_bed
    File? lims_report_format_yml
    File? gene_database_yml
    # QC arguments
    Int min_depth = 10 # default 10; NOTE must explicity set this default for awk to properly interpret threshold
    Float? min_percent_coverage  # default 1.0
    Int? min_read_support # default 10
    Float? min_frequency # default 0.1
    Float? min_percent_loci_covered # default 0.7
    # tNGS-specific arguments
    Boolean tngs_data = false
    Boolean use_err_as_brr = false
    String? tngs_read_support_boundaries # default "10,10"
    String? tngs_frequency_boundaries # default "0.1,0.1"
    # text arguments
    String? sequencing_method
    String? operator
    Map[String, String]? find_and_replace # ex) '{"rifampicin": "rifampin", "mmpR5": "Rv0678", "p.0?": ""}'
    # logging arguments
    Boolean tbp_parser_debug = true
    # WDL runtime arguments
    Int cpu = 1
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:3.0.0-dev"
    Int memory = 4
  }
  command <<<
    # get version
    python3 /tbp-parser/tbp_parser/tbp_parser_main.py --version | tee VERSION

    # run tbp-parser
    python3 /tbp-parser/tbp_parser/tbp_parser_main.py ~{tbprofiler_json} ~{tbprofiler_bam} \
      ~{"--config " + config} \
      ~{"--coverage_bed " + coverage_bed} \
      ~{"--err_coverage_bed " + err_coverage_bed} \
      ~{"--lims_report_format_yml " + lims_report_format_yml} \
      ~{"--gene_database_yml " + gene_database_yml} \
      ~{"--min_depth " + min_depth} \
      ~{"--min_percent_coverage " + min_percent_coverage} \
      ~{"--min_read_support " + min_read_support} \
      ~{"--min_frequency " + min_frequency} \
      ~{"--min_percent_loci_covered " + min_percent_loci_covered} \
      ~{true="--tngs" false="" tngs_data} \
      ~{true="--use_err_as_brr" false="" use_err_as_brr} \
      ~{"--tngs_read_support_boundaries " + tngs_read_support_boundaries} \
      ~{"--tngs_frequency_boundaries " + tngs_frequency_boundaries} \
      ~{"--sequencing_method " + sequencing_method} \
      ~{"--operator " + operator} \
      ~{if defined(find_and_replace) then "--find_and_replace '~{find_and_replace}'" else ""} \
      ~{true="--debug" false="--verbose" tbp_parser_debug} \
      --output_prefix ~{samplename}

    python <<CODE
    import pysam
    import sys

    tbprofiler_bam = "~{tbprofiler_bam}"
    min_depth = int("~{min_depth}")
    is_tngs = "~{tngs_data}" == "true"

    positions = 0
    covered = 0
    depth_sum = 0

    def parse_mpileup(output, min_depth):
        """Parse mpileup output lines, return (positions, covered, depth_sum)"""
        positions = covered = depth_sum = 0
        for line in output.strip().split("\n"):
            if not line:
                continue
            cols = line.split("\t")
            depth = int(cols[3])
            positions += 1
            if depth >= min_depth:
                covered += 1
                depth_sum += depth
        return positions, covered, depth_sum

    if is_tngs:
        # extract chromosome name from BAM
        bam = pysam.AlignmentFile(tbprofiler_bam, "rb")
        chromosome = bam.references[0]
        bam.close()

        # iterate BED file regions (1-based coordinates) and run mpileup for each region, accumulating coverage metrics
        with open("~{coverage_bed}") as bed:
            for line in bed:
                line = line.strip()
                if not line:
                    continue
                cols = line.split("\t")
                start, stop = cols[1], cols[2]
                output = pysam.mpileup(
                    "-a",
                    "--count-orphans",
                    "--min-BQ", "0",
                    "-r", f"{chromosome}:{start}-{stop}",
                    tbprofiler_bam,
                )

                p, c, d = parse_mpileup(output, min_depth)
                positions += p
                covered += c
                depth_sum += d

    else:
        # WGS: single mpileup over entire BAM
        output = pysam.mpileup(
            "-a",
            "--count-orphans",
            "--min-BQ", "0",
            tbprofiler_bam,
        )
        positions, covered, depth_sum = parse_mpileup(output, min_depth)

        if positions != 4411532:
            print(f"WARNING: Expected 4411532 positions (geneome length) for M. tuberculosis H37Rv reference genome, but found {positions} positions in mpileup output")

    print(f"DEBUG: Number of positions: {positions}")
    print(f"DEBUG: Number of covered positions: {covered}")
    print(f"DEBUG: Total covered depth: {depth_sum}")

    if positions == 0:
        print("WARNING: No positions found in mpileup output. Setting genome percent coverage and average depth to 0")
        genome_pc = 0.0
        avg_depth = 0.0
    else:
        genome_pc = (covered / positions) * 100
        avg_depth = depth_sum / positions

    with open("GENOME_PC", "w") as f:
        print(genome_pc, file=f)
    with open("AVG_DEPTH", "w") as f:
        print(avg_depth, file=f)

    print(f"DEBUG: GENOME_PC = {genome_pc}%")
    print(f"DEBUG: AVG_DEPTH = {avg_depth}")
    CODE
  >>>
  output {
    File tbp_parser_looker_report_csv = "~{samplename}.looker_report.csv"
    File tbp_parser_laboratorian_report_csv = "~{samplename}.laboratorian_report.csv"
    File tbp_parser_lims_report_csv = "~{samplename}.lims_report.csv"
    File tbp_parser_lims_report_transposed_csv = "~{samplename}.lims_report.transposed.csv"
    File tbp_parser_locus_coverage_report = "~{samplename}.locus_coverage_report.csv"
    File? tbp_parser_target_coverage_report = "~{samplename}.target_coverage_report.csv"
    Float tbp_parser_genome_percent_coverage = read_float("GENOME_PC")
    Float tbp_parser_average_genome_depth = read_float("AVG_DEPTH")
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
