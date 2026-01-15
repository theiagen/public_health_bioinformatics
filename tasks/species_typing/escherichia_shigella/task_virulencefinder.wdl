version 1.0

task virulencefinder {
  input {
    File? assembly
    ## only allowing assembly input as CGE web app functionality for read inputs cannot be recreated with current staphb docker
    #File? read1
    #File? read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/virulencefinder:2.0.4"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
    Float? min_percent_coverage
    Float? min_percent_identity
    String database = "virulence_ecoli"
    # determine what version of the software to run
    #Boolean paired_end
    #Boolean assembly_only
    #Boolean ont_data
  }
  command <<<
    # capture date and version
    date | tee DATE

    # create temporary intermediate file directory
    mkdir tmp

    # run virulence finder
    #  -i input file(s)
    #  -o output directory
    #  -x extended output (needed to make the results_tab file)
    #  -tmp temporary directory for intermediate results
    #  -l sets threshold for min coverage
    #  -t sets threshold for min blast identity
    #  -d sets a specific database

    # if {assembly_only} || ! {paired_end} || {ont_data}; then
      # run assembly version if SE, ONT, or assembly
    virulencefinder.py \
      -i ~{assembly} \
      -o . \
      -tmp tmp \
      ~{'-l ' + min_percent_coverage} \
      ~{'-t ' + min_percent_identity} \
      ~{'-d' + database} \
      -x 
    # else 
    #   # run the reads version if PE data
    #   virulencefinder.py \
    #     -i {read1} {read2} \
    #     -o . \
    #     -tmp tmp \
    #     ~{'-l ' + min_percent_coverage} \
    #     ~{'-t ' + min_percent_identity} \
    #     ~{'-d' + database} \
    #     -x
    # fi
   
    # rename output file
    mv results_tab.tsv ~{samplename}_results_tab.tsv

    # parse 2nd column (Virulence factor) of _tab.tsv file into a comma-separated string
    tail -n+2 ~{samplename}_results_tab.tsv | awk '{print $2}' | uniq | paste -s -d, - | tee VIRULENCE_FACTORS

    # if virulencefinder fails/no hits are found, the results_tab.tsv file will still exist but only be the header
    # check to see if VIRULENCE_FACTORS is just whitespace
    # if so, say that no virulence factors were found instead
    if ! grep -q '[^[:space:]]' VIRULENCE_FACTORS ; then 
      echo "No virulence factors found" | tee VIRULENCE_FACTORS
    fi
  >>>
  output {
    File virulencefinder_report_tsv = "~{samplename}_results_tab.tsv"
    String virulencefinder_docker = docker
    String virulencefinder_hits = read_string("VIRULENCE_FACTORS")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible:  0
  }
}
