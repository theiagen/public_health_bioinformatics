version 1.0

task tbp_parser {
  input {
    String samplename
    # required arguments
    File tbprofiler_json
    File tbprofiler_bam
    File tbprofiler_bai
    File coverage_bed
    # file arguments
    File? tbprofiler_db_bed # the BED file describing the TBProfiler database the variants were called against
    File? config
    File? err_coverage_bed
    File? lims_report_format_yml
    File? gene_database_yml
    # QC arguments
    Int min_depth = 10 # default 10; NOTE must explicity set this default in order to determine GENOME_PC and AVERAGE_DEPTH
    Float? min_percent_coverage  # default 1.0
    Int? min_read_support # default 10
    Float? min_frequency # default 0.1
    Float? min_percent_loci_covered # default 0.7
    Boolean skip_input_validation = false # skip checking that every gene/drug pair in the BED and LIMS files is present in the gene database
    # tNGS-specific arguments
    Boolean tngs_data = false
    Boolean use_err_for_qc = false
    Boolean resolve_overlapping_regions = false
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
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:3.1.0-dev"
    Int memory = 16
  }
  command <<<
    # get version
    tbp-parser --version | tee VERSION

    # tbp-parser v3.1.0 requires `--gene_database_yml`
    # Build one from the TBProfiler database BED if one was not provided.
    gene_database_yml="~{gene_database_yml}"
    if [[ ! -s "${gene_database_yml}" ]]; then
      if [[ ! -s "~{tbprofiler_db_bed}" ]]; then
        echo "ERROR: no '--gene_database_yml' provided and no 'tbprofiler_db_bed' to build one from; supply one of the two"
        exit 1
      fi

      echo "No '--gene_database_yml' provided; building one from the TBProfiler database BED"
      tbp-parser build_gene_db \
        --db_bed "~{tbprofiler_db_bed}" \
        --output "gene_db.yml" \
        ~{true="--debug" false="" tbp_parser_debug}

      gene_database_yml="gene_db.yml"
    fi

    # tbp-parser v3.1.0 requires `--lims_report_format_yml`
    # Build one from the gene database resolved above if one was not provided.
    lims_report_format_yml="~{lims_report_format_yml}"
    if [[ ! -s "${lims_report_format_yml}" ]]; then
      echo "No '--lims_report_format_yml' provided; deriving one from ${gene_database_yml}"
      tbp-parser build_lims_fmt \
        --gene_database_yml "${gene_database_yml}" \
        --output lims_report_fmt.yml \
        ~{true="--debug" false="" tbp_parser_debug}

      lims_report_format_yml="lims_report_fmt.yml"
    fi

    # run tbp-parser
    tbp-parser parse \
      --input_json ~{tbprofiler_json} \
      --input_bam ~{tbprofiler_bam} \
      --coverage_bed ~{coverage_bed} \
      --gene_database_yml "${gene_database_yml}" \
      --lims_report_format_yml "${lims_report_format_yml}" \
      ~{"--config " + config} \
      ~{"--err_coverage_bed " + err_coverage_bed} \
      ~{"--min_depth " + min_depth} \
      ~{"--min_percent_coverage " + min_percent_coverage} \
      ~{"--min_read_support " + min_read_support} \
      ~{"--min_frequency " + min_frequency} \
      ~{"--min_percent_loci_covered " + min_percent_loci_covered} \
      ~{true="--skip_input_validation" false="" skip_input_validation} \
      ~{true="--tngs" false="" tngs_data} \
      ~{true="--use_err_for_qc" false="" use_err_for_qc} \
      ~{true="--resolve_overlapping_regions" false="" resolve_overlapping_regions} \
      ~{"--tngs_read_support_boundaries " + tngs_read_support_boundaries} \
      ~{"--tngs_frequency_boundaries " + tngs_frequency_boundaries} \
      ~{"--sequencing_method " + sequencing_method} \
      ~{"--operator " + operator} \
      ~{if defined(find_and_replace) then "--find_and_replace '~{find_and_replace}'" else ""} \
      ~{true="--debug" false="" tbp_parser_debug} \
      --output_prefix ~{samplename}

    python <<CODE
    import pysam

    tbprofiler_bam = "~{tbprofiler_bam}"
    coverage_bed = "~{coverage_bed}"
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
        with open(coverage_bed) as bed:
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
    File tbp_parser_locus_coverage_report_csv = "~{samplename}.locus_coverage_report.csv"
    File? tbp_parser_target_coverage_report_csv = "~{samplename}.target_coverage_report.csv"
    File? tbp_parser_generated_gene_database_yml = "gene_db.yml"
    File? tbp_parser_generated_lims_report_format_yml = "lims_report_fmt.yml"
    File tbp_parser_log = "~{samplename}.log"
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
    maxRetries: 2
    preemptible: 1
  }
}
