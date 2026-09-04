version 1.0

task dorado {
  input {
    File read1
    File unpolished_fasta
    String samplename

    String? dorado_model = "dna_r10.4.1_e8.2_400bps_sup@v5.0.0_polish_rl"
    Boolean ignore_read_groups = true
    Boolean auto_detect_model = true
    Boolean use_bacteria = true

    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/dorado:2.1.2"
  }
  command <<<
    set -euo pipefail

    dorado --version 2&> VERSION

    echo "DEBUG: moving the unpolished_fasta to the local dir and building the index locally"
    cp ~{unpolished_fasta} unpolished.fasta
    samtools faidx unpolished.fasta

    echo "DEBUG: aligning the reads to the unpolished fasta" > ~{samplename}.log
    dorado aligner \
      unpolished.fasta \
      ~{read1} \
      -t ~{cpu} \
      --no-sort \
      > ~{samplename}.bam

    echo "DEBUG: sorting and indexing the alignment" >> ~{samplename}.dorado_log
    samtools sort \
      -@ ~{cpu} \
      ~{samplename}.bam > ~{samplename}_sorted.bam
    samtools index ~{samplename}_sorted.bam


    if ~{auto_detect_model}; then
      echo "DEBUG: model will be auto-detected by Dorado" >> ~{samplename}.log
      model=""
    else
      echo "DEBUG: model auto-detection was not enabled; using the user-specified or default model (~{dorado_model}" >> ~{samplename}.log
      model="--model-override /models/~{dorado_model}"

      # confirm model is in /models dir
      if [[ ! -d /models/~{dorado_model} ]]; then
        echo "WARNING: model (~{dorado_model}) is not included in the model directory; downloading the model" >> ~{samplename}.log

        dorado download \
          --model ~{dorado_model} \
          --models-directory /models
      fi
    fi

    dorado polish \
      ~{samplename}_sorted.bam \
      unpolished.fasta \
      -v -t ~{cpu} \
      --models-directory /models \
      ${model} \
      ~{true="--ignore-read-groups" false="" ignore_read_groups} \
      ~{true="--bacteria" false="" use_bacteria} \
      > ~{samplename}_polished.fasta \
      2> >(tee ~{samplename}.log >&2)

  >>>
  output {
    File polished_fasta = "~{samplename}_polished.fasta"
    File dorado_log = "~{samplename}.log"
    String dorado_version = read_string("VERSION")
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: memory + " GB"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
