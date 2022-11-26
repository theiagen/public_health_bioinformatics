version 1.0

task busco {
  meta {
    description: "Run BUSCO on assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "ezlabgva/busco:v5.3.2_cv1"
    Boolean eukaryote = false
  }
  command <<<
    # get version
    busco --version | tee "VERSION"
 
    # run busco
    # -i input assembly
    # -m geno for genome input
    # -o output file tag
    # --auto-lineage-euk looks at only eukaryotic organisms
    # --auto-lineage-prok looks at only prokaryotic organisms; default
    busco \
      -i ~{assembly} \
      -m geno \
      -o ~{samplename} \
      ~{true='--auto-lineage-euk' false='--auto-lineage-prok' eukaryote}

    # check for existence of output file; otherwise display a string that says the output was not created
    if [ -f ~{samplename}/short_summary.specific.*.~{samplename}.txt ]; then

      # grab the database version and format it according to BUSCO recommendations
      cat ~{samplename}/short_summary.specific.*.~{samplename}.txt | grep "dataset is:" | cut -d' ' -f 6,9 | sed 's/,//' | sed 's/ / (/' | sed 's/$/)/' | tee DATABASE
      
      # extract the results string
      cat ~{samplename}/short_summary.specific.*.~{samplename}.txt | grep "C:" | tee BUSCO_RESULTS

      cp ~{samplename}/short_summary.specific.*.~{samplename}.txt ~{samplename}_busco-summary.txt
    else
      echo "BUSCO FAILED" | tee BUSCO_RESULTS
      echo "NA" > DATABASE
    fi
  >>>
  output {
    String busco_version = read_string("VERSION")
    String busco_database = read_string("DATABASE")
    String busco_results = read_string("BUSCO_RESULTS")
    File? busco_report = "~{samplename}_busco-summary.txt"
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}