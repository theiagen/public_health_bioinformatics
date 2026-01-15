version 1.0

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

    Int memory = 30
    Int disk_size = 100
    Int cpu = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/augur:31.5.0"
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
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem3_ssd2_x4"
    preemptible: 0
    maxRetries: 3
  }
}