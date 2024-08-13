version 1.0

# task adapted with some modifications from Neranjan Perera's GAS_identification task_emmtyper.wdl
# https://github.com/neranjan007/GAS_identification/blob/454ef3b0cc8a90950b48342cde87136962f9adb1/tasks/task_emmtyper.wdl

task emmtyper {
  meta {
    description: "emm-typing of Streptococcus pyogenes assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/emmtyper:0.2.0--py_0"
    Int cpu = 2
    Int memory = 8
    Int disk_size = 50

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
      ~{'--workflow ' + wf} \
      ~{'--cluster-distance ' + cluster_distance} \
      ~{'--percent-identity ' + percid} \
      ~{'--culling-limit ' + culling_limit} \
      ~{'--mismatch ' + mismatch} \
      ~{'--align-diff ' + align_diff} \
      ~{'--gap ' + gap} \
      ~{'--min-perfect ' + min_perfect} \
      ~{'--min-good ' + min_good} \
      ~{'--max-size ' + max_size} \
      --output-format verbose \
      ~{assembly} \
      > ~{samplename}_emmtyper.tsv

      # emm type is in column 4 for verbose output format
      awk -F "\t" '{print $4}' ~{samplename}_emmtyper.tsv > EMM_TYPE
  >>>
  output {
    String emmtyper_emm_type = read_string("EMM_TYPE")
    File emmtyper_results_tsv = "~{samplename}_emmtyper.tsv"
    String emmtyper_version = read_string("VERSION")
    String emmtyper_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
