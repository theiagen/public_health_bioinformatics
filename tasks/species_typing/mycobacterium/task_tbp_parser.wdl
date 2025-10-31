version 1.0

task tbp_parser {
  input {
    File tbprofiler_json
    File tbprofiler_bam
    File tbprofiler_bai
    String samplename

    File? config
    String? sequencing_method
    String? operator

    Int min_depth = 10 # default 10
    Float? min_frequency # default 0.1
    Int? min_read_support # default 10
    
    Float? min_percent_coverage # default 100 (--min_percent_coverage)
    File? coverage_regions_bed
  
    Boolean add_cycloserine_lims = false
    Boolean tbp_parser_debug = true
    Boolean tngs_data = false    

    Float? rrs_frequency # default 0.1
    Int? rrs_read_support # default 10
    Float? rrl_frequency # default 0.1
    Int? rrl_read_support # default 10
    Float? rpob449_frequency # default 0.1
    Float? etha237_frequency # default 0.1
    File? expert_rule_regions_bed
    
    Int cpu = 1
    Int disk_size = 100   
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:2.6.0"
    Int memory = 4
  }
  command <<<
    # explicitly set covereage regions bed file defaults for wgs/tngs data
    if [[ -z "~{coverage_regions_bed}" ]]; then
      if [[ "~{tngs_data}" == "true" ]]; then
        coverage_regions_bed="/tbp-parser/data/tngs-primer-regions.bed"
      else
        coverage_regions_bed="/tbp-parser/data/tbdb-modified-regions.bed"
      fi
    else
      coverage_regions_bed="~{coverage_regions_bed}"
    fi

    # get version
    python3 /tbp-parser/tbp_parser/tbp_parser.py --version | tee VERSION

    # run tbp-parser
    python3 /tbp-parser/tbp_parser/tbp_parser.py ~{tbprofiler_json} ~{tbprofiler_bam} \
      ~{"--config " + config} \
      ~{"--sequencing_method " + sequencing_method} \
      ~{"--operator " + operator} \
      ~{"--min_depth " + min_depth} \
      ~{"--min_frequency " + min_frequency} \
      ~{"--min_read_support " + min_read_support} \
      ~{"--min_percent_coverage " + min_percent_coverage} \
      --coverage_regions $coverage_regions_bed \
      ~{"--tngs_expert_regions " + expert_rule_regions_bed} \
      ~{"--rrs_frequency " + rrs_frequency} \
      ~{"--rrs_read_support " + rrs_read_support} \
      ~{"--rrl_frequency " + rrl_frequency} \
      ~{"--rrl_read_support " + rrl_read_support} \
      ~{"--rpob449_frequency " + rpob449_frequency} \
      ~{"--etha237_frequency " + etha237_frequency} \
      --output_prefix ~{samplename} \
      ~{true="--debug" false="--verbose" tbp_parser_debug} \
      ~{true="--tngs" false="" tngs_data} \
      ~{true="--add_cs_lims" false="" add_cycloserine_lims}

    # set default genome percent coverage and average depth to 0 to prevent failures
    echo 0.0 > GENOME_PC
    echo 0.0 > AVG_DEPTH

    if [[ "~{tngs_data}" == "true" ]]; then
      # extract chromosome name from bam file
      chromosome=$(samtools idxstats ~{tbprofiler_bam} | cut -f 1 | head -1)

      # initialize counters
      filtered_primer_positions=0
      total_primer_positions=0
      total_primer_depth=0

      # iterate through file and calculate coverage for each row in the bedfile (1-based)
      while read -r line; do
        # pull out the important fields from the bedfile
        start=$(echo "$line" | cut -f 2)
        stop=$(echo "$line" | cut -f 3)
        gene=$(echo "$line" | cut -f 4)

        # count filtered positions, total positions, and total depth for the region
        read \
          filtered_count \
          total_count \
          total_depth \
          < <(
            samtools depth -a -J -r "$chromosome:$start-$stop" "~{tbprofiler_bam}" | \
            awk -v min_depth="~{min_depth}" '
              BEGIN { filtered_count=0; total_count=0; total_depth=0 }
              {
                total_count++;
                if ($3 >= min_depth) {
                  filtered_count++; total_depth += $3
                }
              }
              END { print filtered_count, total_count, total_depth }
            '
          )

        # accumulate counts across all primer regions
        filtered_primer_positions=$((filtered_primer_positions + filtered_count))
        total_primer_positions=$((total_primer_positions + total_count))
        total_primer_depth=$((total_primer_depth + total_depth))
      done < "$coverage_regions_bed"

      # prevents division by zero if no coverage across any primer regions
      if [[ $filtered_primer_positions -eq 0 ]]; then
        echo "No coverage across any tNGS regions found."
      else
        python3 -c "print ( ($filtered_primer_positions / $total_primer_positions ) * 100 )" | tee GENOME_PC
        # get average depth for all primer regions
        python3 -c "print ( $total_primer_depth / $filtered_primer_positions )" | tee AVG_DEPTH
      fi

    else
      # get genome percent coverage for the entire reference genome length over min_depth
      genome=$(samtools depth -a -J ~{tbprofiler_bam} | awk -F "\t" -v min_depth=~{min_depth} '{if ($3 >= min_depth) print;}' | wc -l )
      python3 -c "print ( ($genome / 4411532 ) * 100 )" | tee GENOME_PC
      # get genome average depth
      samtools depth -a -J ~{tbprofiler_bam} | awk -F "\t" '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }' | tee AVG_DEPTH
    fi

    # add sample id to the beginning of the coverage report
    awk '{s=(NR==1)?"Sample_accession_number,":"~{samplename},"; $0=s$0}1' ~{samplename}.percent_gene_coverage.csv > tmp.csv && mv -f tmp.csv ~{samplename}.percent_gene_coverage.csv
  >>>
  output {
    File tbp_parser_looker_report_csv = "~{samplename}.looker_report.csv"
    File tbp_parser_laboratorian_report_csv = "~{samplename}.laboratorian_report.csv"
    File tbp_parser_lims_report_csv = "~{samplename}.lims_report.csv"
    File tbp_parser_coverage_report = "~{samplename}.percent_gene_coverage.csv"
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
