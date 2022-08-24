version 1.0

task emmtyper {
  meta {
    description: "emm-typing of Streptococcus pyogenes assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/emmtyper:0.2.0--py_0"
    Int? cpu = 2

    # Parameters
    # --workflow [blast|pcr]      Choose workflow  [default: blast]
    # --cluster-distance INTEGER  Distance between cluster of matches to consider as different clusters.  [default: 500]
    # --percent-identity INTEGER      [BLAST] Minimal percent identity of sequence.  [default: 95]
    # --culling-limit INTEGER         [BLAST] Total hits to return in a position. [default: 5]
    # --mismatch INTEGER              [BLAST] Threshold for number of mismatch to allow in BLAST hit.  [default: 4]
    # --align-diff INTEGER            [BLAST] Threshold for difference between alignment length and subject length in BLAST hit.  [default: 5]
    # --gap INTEGER                   [BLAST] Threshold gap to allow in BLAST hit. [default: 2]
    # --min-perfect INTEGER           [isPcr] Minimum size of perfect match at 3' primer end.  [default: 15]
    # --min-good INTEGER              [isPcr] Minimum size where there must be 2 matches for each mismatch.  [default: 15]
    # --max-size INTEGER              [isPcr] Maximum size of PCR product. [default: 2000]

    String wf = "blast"
    Int cluster_distance = 500
    Int percid = 95
    Int culling_limit = 5
    Int mismatch = 4
    Int align_diff = 5
    Int gap = 2
    Int min_perfect = 15
    Int min_good = 15
    Int max_size = 2000
  }
  command <<<
    echo $(emmtyper --version 2>&1) | sed 's/^.*emmtyper v//' | tee VERSION
    emmtyper \
      ~{'--workflow' + wf} \
      ~{'--cluster-distance' + cluster_distance} \
      ~{'--percent-identity' + percid} \
      ~{'--culling-limit' + culling_limit} \
      ~{'--mismatch' + mismatch} \
      ~{'--align-diff' + align_diff} \
      ~{'--gap' + gap} \
      ~{'--min-perfect' + min_perfect} \
      ~{'--min-good' + min_good} \
      ~{'--max-size' + max_size} \
      ~{assembly} \
      > ~{samplename}.tsv
  >>>
  output {
    File emmtyper_results = "~{samplename}.tsv"
    String emmtyper_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
