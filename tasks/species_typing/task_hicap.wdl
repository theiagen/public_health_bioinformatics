version 1.0

task hicap {
  meta {
    description: "Identify cap locus serotype and structure in your Haemophilus influenzae assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/hicap:1.0.3--py_0"
    Int? cpu = 4

    # Parameters
    # --gene_coverage GENE_COVERAGE                 Minimum percentage coverage to consider a single gene complete. [default: 0.80]
    # --gene_identity GENE_IDENTITY                 Minimum percentage identity to consider a single gene complete. [default: 0.70]
    # --broken_gene_length BROKEN_GENE_LENGTH       Minimum length to consider a broken gene. [default: 60]
    # --broken_gene_identity BROKEN_GENE_IDENTITY   Minimum percentage identity to consider a broken gene. [default: 0.80]
    Float gene_coverage = 0.8
    Float gene_identity = 0.7
    Int broken_gene_length = 60
    Float broken_gene_identity = 0.8
    Boolean full_sequence = false
    Boolean debug = false
  }
  command <<<
    echo $( hicap --version 2>&1 ) | sed 's/^.*hicap //' | tee VERSION
    hicap \
      --query_fp ~{assembly} \
      ~{'--gene_coverage' + gene_coverage} \
      ~{'--gene_identity' + gene_identity} \
      ~{'--broken_gene_length' + broken_gene_length} \
      ~{'--broken_gene_identity' + broken_gene_identity} \
      ~{true="--full_sequence" false="" full_sequence} \
      ~{true="--debug" false="" debug} \
      --threads ~{cpu} \
      -o ./

    if [ ! -f ${samplename}.tsv ]; then
      # No hits, make a file to say so for downstream merging
      echo "isolate<TAB>predicted_serotype<TAB>attributes<TAB>genes_identified<TAB>locus_location<TAB>region_I_genes<TAB>region_II_genes<TAB>region_III_genes<TAB>IS1016_hits" | sed 's/<TAB>/\t/g' > ${samplename}.tsv
      echo "~{samplename}<TAB>cap_not_found<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-<TAB>-" | sed 's/<TAB>/\t/g' >> ~{samplename}.tsv
    else
      sed -i 's/#isolate/isolate/' ~{samplename}.tsv
    fi
  >>>
  output {
    File hicap_results = "~{samplename}.tsv"
    File hicap_genbank = "~{samplename}.gbk"
    File hicap_image = "~{samplename}.svg"
    String hicap_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
