version 1.0

task augur_align {
  input {
    File assembly_fasta
    File reference_fasta
    Boolean fill_gaps = false
    Int cpus = 64
    Int mem_size = 32
    Int disk_size = 750
  }
  command <<<
    # capture version information
    augur version > VERSION

    # run augur align
    augur align \
      --sequences ~{assembly_fasta} \
      --nthreads ~{cpus} \
      --reference-sequence ~{reference_fasta} \
      ~{true="--fill-gaps" false="" fill_gaps}
  >>>
  output {
    File aligned_fasta = "alignment.fasta"
    String augur_version = read_string("VERSION")
  }
  runtime {
    docker: "quay.io/staphb/augur:16.0.3"
    memory: mem_size + " GB"
    cpu :   cpus
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    preemptible: 0
    dx_instance_type: "mem3_ssd1_v2_x36"
    maxRetries: 3
  }
}

task augur_tree {
  input {
    File aligned_fasta
    String build_name
    String method = "iqtree" # possible choices: fasttree, raxml, iqtree
    String substitution_model = "GTR" # only available for iqtree
    File? exclude_sites # file name of one-based sites to exclude for raw tree building
    String? tree_builder_args # additional tree builder arguments
    Boolean override_default_args = false # override default tree builder args instead of augmenting them

    Int cpus = 64
    Int disk_size = 750
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur tree \
      --alignment "~{aligned_fasta}" \
      --output "~{build_name}_~{method}.nwk" \
      --method "~{method}" \
      --substitution-model ~{substitution_model} \
      ~{"--exclude-sites " + exclude_sites} \
      ~{"--tree-builder-args " + tree_builder_args} \
      ~{true="--override-default-args" false="" override_default_args} \
      --nthreads auto
  >>>
  output {
    File aligned_tree  = "~{build_name}_~{method}.nwk"
  }
  runtime {
    docker: "quay.io/staphb/augur:16.0.3"
    memory: "32 GB"
    cpu: cpus
    disks: "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x36"
    preemptible: 0
    maxRetries: 3
  }
}

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
    Int clock_filter_iqd = 4 # remove tips that deviate more than n_iqd interquartile ranges from the root-to-tip vs time regression
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
    docker: "quay.io/staphb/augur:16.0.3"
    memory: "50 GB"
    cpu : 2
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem3_ssd1_v2_x8"
    preemptible: 0
    maxRetries: 3
  }
}

task augur_frequencies {
  input {
    File refined_tree
    File metadata
    String method = "kde" # options: diffusion, kde
    String build_name

    Float? min_date # date to begin frequency calculations
    Float? max_date # date to end frequency calculations
    Int? pivot_interval # number of units between pivots (default: 3)
    String? pivot_interval_units # options: months (default), weeks
    Float? narrow_bandwidth # the bandwidth for the narrow KDE (default: 0.08333333333333333)
    Float? wide_bandwidth # the bandwidth for the wide KDE (default: 0.25)
    Float? proportion_wide # the proportion of the wide bandwidth to use in the KDE mixture model (default: 0.2)
    Float? minimal_frequency # minimal all-time frequencies for a trajectory to be estimates (default: 0.05)
    Float? stiffness # parameter penalizing curbature of the frequency trajectory (default: 10.0)
    Float? inertia # determines how frequencies continue in absence of data: if inertia = 0 (default), go flat, if inertia = 1, continue trend
    Boolean censored = false # calculate censored frequencies at each pivot
    Boolean include_internal_nodes = false # calculate frequencies for internal nodes as well as tips

    Int mem_size = 30
    Int disk_size = 100
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur frequencies \
      --method "~{method}" \
      --tree "~{refined_tree}" \
      --metadata "~{metadata}" \
      ~{'--min-date ' + min_date} \
      ~{'--max-date ' + max_date} \
      ~{'--pivot-interval ' + pivot_interval} \
      ~{'--pivot-interval-units ' + pivot_interval_units} \
      ~{'--narrow-bandwidth ' + narrow_bandwidth} \
      ~{'--wide-bandwidth ' + wide_bandwidth} \
      ~{'--proportion-wide ' + proportion_wide} \
      ~{'--narrow-bandwidth ' + narrow_bandwidth} \
      ~{'--minimal-frequency ' + minimal_frequency} \
      ~{'--stiffness ' + stiffness} \
      ~{'--inertia ' + inertia} \
      ~{true='--censored' false='' censored} \
      ~{true='--include-internal-nodes' false='' include_internal_nodes} \
      --output "~{build_name}_tip-frequencies.json"
  >>>
  output {
    File tip_frequencies_json = "~{build_name}_tip-frequencies.json"
  }
  runtime {
    docker: "quay.io/staphb/augur:16.0.3"
    memory: mem_size + " GB"
    cpu: 4
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem3_ssd2_x4"
    preemptible: 0
    maxRetries: 3
  }
}

task augur_ancestral {
  input {
    File refined_tree
    File aligned_fasta
    String build_name

    String inference = "joint" # options: joint, marginal
    Boolean keep_ambiguous = false # do not infer nucleotides at ambiguous (N) sites
    Boolean infer_ambiguous = false # infer nucleotides at ambiguous sites and replace with most likely
    Boolean keep_overhangs = false # do not infer nucleotides for gaps on either side of the alignment

    Int disk_size = 50
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur ancestral \
      --tree "~{refined_tree}" \
      --alignment "~{aligned_fasta}" \
      --output-node-data "~{build_name}_nt_muts.json" \
      --output-sequences "~{build_name}_ancestral_sequences.fasta" \
      --inference ~{default="joint" inference} \
      ~{true="--keep-ambiguous" false="" keep_ambiguous} \
      ~{true="--infer-ambiguous" false="" infer_ambiguous} \
      ~{true="--keep-overhangs" false="" keep_overhangs} 
  >>>
  output {
    File ancestral_nt_muts_json = "~{build_name}_nt_muts.json"
    File ancestral_sequences = "~{build_name}_ancestral_sequences.fasta"
  }
  runtime {
    docker: "quay.io/staphb/augur:16.0.3"
    memory: "50 GB"
    cpu: 4
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" 
    dx_instance_type: "mem3_ssd1_v2_x8"
    preemptible: 0
    maxRetries: 3
  }
}

task augur_translate {
  input {
    File refined_tree
    File ancestral_nt_muts_json
    File reference_genbank
    String build_name

    File? genes # a file containing list of genes to translate (from nucleotides to amino acids)

    Int disk_size = 50
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur translate \
      --tree "~{refined_tree}" \
      --ancestral-sequences "~{ancestral_nt_muts_json}" \
      --reference-sequence "~{reference_genbank}" \
      ~{"--genes " + genes} \
      --output-node-data "~{build_name}_aa_muts.json"
  >>>
  output {
    File translated_aa_muts_json = "~{build_name}_aa_muts.json"
  }
  runtime {
    docker: "quay.io/staphb/augur:16.0.3"
    memory: "2 GB"
    cpu : 1
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 0
    maxRetries: 3
  }
}

task augur_clades {
  input {
    File refined_tree
    File ancestral_nt_muts_json
    File translated_aa_muts_json
    File reference_fasta
    String build_name

    File clades_tsv # tsv file containing clade definitions by amino acid
    Int disk_size = 50
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur clades \
      --tree "~{refined_tree}" \
      --mutations "~{ancestral_nt_muts_json}" "~{translated_aa_muts_json}" \
      --reference "~{reference_fasta}" \
      --clades "~{clades_tsv}" \
      --output-node-data "~{build_name}_clades.json"
  >>>
  output {
    File clade_assignments_json = "~{build_name}_clades.json"
  }
  runtime {
    docker: "quay.io/staphb/augur:16.0.3"
    memory: "2 GB"
    cpu :   1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" 
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 0
    maxRetries: 3
  }
}

task augur_export {
  input {
    File refined_tree
    File metadata
    Array[File] node_data_jsons
    String build_name

    File? auspice_config # auspice configuration file
    String? title # title to be displayed by Auspice
    File? description_md # markdown file with description of build and/or acknowledgements
    File? colors_tsv # custom color definitions, one per line
    File? lat_longs_tsv # latitudes and longitudes for geography traits
    Boolean include_root_sequence = false # export an additional json containing the root sequence used to identify mutations
  
    Int disk_size = 100
  }
  command <<<
    augur export v2 \
      --tree ~{refined_tree} \
      --metadata ~{metadata} \
      --node-data ~{sep=' ' node_data_jsons} \
      --output ~{build_name}_auspice.json \
      ~{"--auspice-config " + auspice_config} \
      ~{"--title " + title} \
      ~{"--description " + description_md} \
      ~{"--colors " + colors_tsv} \
      ~{"--lat-longs " + lat_longs_tsv} \
      ~{true="--include-root-sequence " false=""  include_root_sequence}
  >>>
  output {
    File auspice_json = "~{build_name}_auspice.json"
    File? root_sequence_json = "~{build_name}_auspice_root-sequence.json"
  }
  runtime {
    docker: "quay.io/staphb/augur:16.0.3"
    memory: "64 GB"
    cpu :   4
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem3_ssd1_v2_x4"
    preemptible: 0
    maxRetries: 3
  }
}

