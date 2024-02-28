version 1.0

task ngmaster {
  meta {
    description: "Multi-antigen sequence typing for Neisseria gonorrhoeae"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ngmaster:1.0.0"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    ngmaster --version 2>&1 | sed 's/^.*ngmaster //' | tee VERSION

    # run ngmaster on input assembly
    # unfortunately ngmaster 1.0.0 fails when either mincov or minid flags are supplied (this is with different install strategies too - bioconda & manually)
    # so we're forced to stick with default minid of 90 and mincov of 10. https://github.com/MDU-PHL/ngmaster/issues/39
    # ngmaster --comments also does not work
    ngmaster \
      ~{assembly} \
      > ~{samplename}.ngmaster.tsv

    # parse output TSV
    # first one is tricky since MLSTs are in the 3rd column, separated by a /
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $3}' | cut -d '/' -f 1 | tee NGMAST_SEQUENCE_TYPE
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $3}' | cut -d '/' -f 2 | tee NGSTAR_SEQUENCE_TYPE
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $4}' | tee NGMAST_PORB
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $5}' | tee NGMAST_TBPB
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $6}' | tee NGSTAR_PENA
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $7}' | tee NGSTAR_MTRR
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $8}' | tee NGSTAR_PORB
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $9}' | tee NGSTAR_PONA
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $10}' | tee NGSTAR_GYRA
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $11}' | tee NGSTAR_PARC
    tail -n 1 ~{samplename}.ngmaster.tsv | awk '{print $12}' | tee NGSTAR_23S
  >>>
  output {
    File ngmaster_tsv = "~{samplename}.ngmaster.tsv"
    String ngmaster_version = read_string("VERSION")
    # NG-MAST scheme's MLST and alleles (only 2 loci)
    String ngmaster_ngmast_sequence_type = read_string("NGMAST_SEQUENCE_TYPE")
    String ngmaster_ngmast_porB_allele = read_string("NGMAST_PORB")
    String ngmaster_ngmast_tbpB_allele = read_string("NGMAST_TBPB")
    # NG-STAR scheme's MLST and alleles (7 loci)
    String ngmaster_ngstar_sequence_type = read_string("NGSTAR_SEQUENCE_TYPE")
    String ngmaster_ngstar_penA_allele = read_string("NGSTAR_PENA")
    String ngmaster_ngstar_mtrR_allele = read_string("NGSTAR_MTRR")
    String ngmaster_ngstar_porB_allele = read_string("NGSTAR_PORB")
    String ngmaster_ngstar_ponA_allele = read_string("NGSTAR_PONA")
    String ngmaster_ngstar_gyrA_allele = read_string("NGSTAR_GYRA")
    String ngmaster_ngstar_parC_allele = read_string("NGSTAR_PARC")
    String ngmaster_ngstar_23S_allele = read_string("NGSTAR_23S")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}
