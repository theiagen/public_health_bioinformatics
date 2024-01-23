version 1.0

task augur_refine {
  input {
    File aligned_fasta
    File draft_augur_tree
    File metadata
    String build_name

    Int? gen_per_year # number of generations per year (default: 50)
    Float? clock_rate # fixed clock rate
    Float? clock_std_dev # standard deviation of the fixed clock_rate estimate
    Boolean keep_root = true # do not reroot the tree; use it as-is (overrides anything specified by --root)
    String? root # rooting mechanism ("best", "least-squares", "min_dev", "oldest", etc)
    Boolean? covariance # account for covariation when estimating rates and/or rerooting (default: true)
    Boolean keep_polytomies = false # don't attempt to resolve polytomies
    Int? precision # precision used to determine the number of grid points; options: 0 (rough) to 3 (ultra fine); default 'auto'
    Boolean date_confidence = true # calculate confidence intervals for node dates
    String date_inference = "marginal" # assign internal nodes to their marginally most likley dates (joint, marginal)
    String? branch_length_inference # branch length mode of treetime to use (auto, joint, marginal, input; default: auto)
    String? coalescent # coalescent time scale in units of inverse clock rate (float), optimize as scalar ("opt") or skyline (skyline)
    Int? clock_filter_iqd = 4 # remove tips that deviate more than n_iqd interquartile ranges from the root-to-tip vs time regression
    String divergence_units = "mutations" # units in which sequence divergences is exported ("mutations" or "mutations-per-site")

    Int disk_size = 100
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur refine \
      --tree "~{draft_augur_tree}" \
      --alignment "~{aligned_fasta}" \
      --metadata "~{metadata}" \
      --output-tree "~{build_name}_timetree.nwk" \
      --output-node-data "~{build_name}_branch_lengths.json" \
      --timetree \
      ~{"--clock-rate " + clock_rate} \
      ~{"--clock-std-dev " + clock_std_dev} \
      ~{"--coalescent " + coalescent} \
      ~{"--clock-filter-iqd " + clock_filter_iqd} \
      ~{"--gen-per-year " + gen_per_year} \
      ~{"--root " + root} \
      ~{"--precision " + precision} \
      ~{"--date-inference " + date_inference} \
      ~{"--branch-length-inference " + branch_length_inference} \
      ~{"--divergence-units " + divergence_units} \
      ~{true="--covariance" false="--no-covariance" covariance} \
      ~{true="--keep-root" false="" keep_root} \
      ~{true="--keep-polytomies" false="" keep_polytomies} \
      ~{true="--date-confidence" false="" date_confidence} 
  >>>
  output {
    File refined_tree   = "~{build_name}_timetree.nwk"
    File branch_lengths = "~{build_name}_branch_lengths.json"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/biocontainers/augur:22.0.2--pyhdfd78af_0"
    memory: "50 GB"
    cpu : 2
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem3_ssd1_v2_x8"
    preemptible: 0
    maxRetries: 3
  }
}
