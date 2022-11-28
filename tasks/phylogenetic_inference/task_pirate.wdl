version 1.0

task pirate {
  input {
    Array[File] prokka_gff
    String cluster_name
    Boolean? align # align all genes and produce core/pangenome alignments
    String? steps = "50,60,70,80,90,95,98" # % identity thresholds to use for pangenome construction [default: 50,60,70,80,90,95,98]
    String? features = "CDS" # features to use for pangenome construction [default: CDS]
    Boolean? nucl = false # CDS are not translated to AA sequence [default: off]
    String? panopt # additional arguments to pass to pangenome_contruction
    Int memory = 32
    Int cpu = 4
    String docker_image = "quay.io/biocontainers/pirate:1.0.5--hdfd78af_0"
  }
  command <<<
  
  # date and version control
  date | tee DATE
  PIRATE -v | tee VERSION

  # pirate requires the directory containing the gff files as input
  mkdir INPUT_DIR
  ln -s ~{sep=' ' prokka_gff} INPUT_DIR

  # run pirate on input gff
  PIRATE \
  --input INPUT_DIR \
  --output PIRATE \
  ~{'--steps ' + steps} \
  ~{'--features ' + features} \
  ~{true="--nucl" false="" nucl} \
  ~{true="--align" false="" align} \
  ~{'--pan-opt ' + panopt} \
  ~{'--threads ' + cpu} 
  
  # rename outputs with cluster name 
  mv PIRATE/PIRATE.pangenome_summary.txt PIRATE/~{cluster_name}_pangenome_summary.txt
  mv PIRATE/PIRATE.log PIRATE/~{cluster_name}.log
  mv PIRATE/PIRATE.gene_families.ordered.tsv PIRATE/~{cluster_name}_gene_families.ordered.tsv
  mv PIRATE/PIRATE.unique_alleles.tsv PIRATE/~{cluster_name}_unique_alleles.tsv
  mv PIRATE/binary_presence_absence.fasta PIRATE/~{cluster_name}_binary_presence_absence.fasta
  mv PIRATE/binary_presence_absence.nwk PIRATE/~{cluster_name}_binary_presence_absence.nwk
  mv PIRATE/pangenome.gfa PIRATE/~{cluster_name}_pangenome.gfa
  mv PIRATE/pangenome_alignment.fasta PIRATE/~{cluster_name}_pangenome_alignment.fasta
  mv PIRATE/pangenome_alignment.gff PIRATE/~{cluster_name}_pangenome_alignment.gff
  mv PIRATE/core_alignment.fasta PIRATE/~{cluster_name}_core_alignment.fasta
  mv PIRATE/core_alignment.gff PIRATE/~{cluster_name}_core_alignment.gff

  >>>
  output {
    File pirate_pangenome_summary = "PIRATE/~{cluster_name}_pangenome_summary.txt"
    File pirate_gene_families_ordered = "PIRATE/~{cluster_name}_gene_families.ordered.tsv"
    File pirate_unique_alleles = "PIRATE/~{cluster_name}_unique_alleles.tsv"
    File pirate_binary_fasta = "PIRATE/~{cluster_name}_binary_presence_absence.fasta"
    File pirate_binary_tree = "PIRATE/~{cluster_name}_binary_presence_absence.nwk"
    File pirate_pangenome_gfa = "PIRATE/~{cluster_name}_pangenome.gfa" 
    File pirate_pangenome_alignment_fasta = "PIRATE/~{cluster_name}_pangenome_alignment.fasta" 
    File pirate_pangenome_alignment_gff = "PIRATE/~{cluster_name}_pangenome_alignment.gff" 
    File pirate_core_alignment_fasta = "PIRATE/~{cluster_name}_core_alignment.fasta" 
    File pirate_core_alignment_gff = "PIRATE/~{cluster_name}_core_alignment.gff" 
    String pirate_docker_image = docker_image
  } 
  runtime {
    docker: "~{docker_image}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
