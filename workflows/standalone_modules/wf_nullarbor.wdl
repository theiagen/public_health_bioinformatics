version 1.0

import "../../tasks/utilities/task_nullarbor.wdl" as nullarbor

workflow nullarbor_workflow {
  meta {
    description: "Nullarbor workflow"
  }
  input {
    String? run_name
    File ref_genome
    Array[File] read1
    Array[File] read2
    Array[String] samplename
    String? tree_builder
    String? assembler
    String? taxoner
    String? docker
    Int? memory
    Int? cpu
    String? kraken1_db
    String? kraken2_db
    String? centrifuge_db
  }
  call nullarbor.nullarbor_tsv as nullarbor_task {
    input:
      run_name = run_name,
      ref_genome = ref_genome,
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      tree_builder = tree_builder,
      assembler = assembler,
      taxoner = taxoner,
      docker = docker,
      memory = memory,
      cpu = cpu,
      kraken1_db = kraken1_db,
      kraken2_db = kraken2_db,
      centrifuge_db = centrifuge_db
  }
  output {
    # Version Capture
    String nullarbor_version =  nullarbor_task.nullarbor_version
    String nullarbor_docker = nullarbor_task.nullarbor_docker
    String nullarbor_analysis_date = nullarbor_task.analysis_date
    File nullarbor_report = nullarbor_task.nullarbor_report
    File nullarbor_outputdir = nullarbor_task.nullarbor_output_dir
  }
}