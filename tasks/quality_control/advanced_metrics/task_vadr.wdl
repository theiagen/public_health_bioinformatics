version 1.0

task vadr {
  meta {
    description: "Runs NCBI's Viral Annotation DefineR for annotation and QC. See https://github.com/ncbi/vadr/wiki/Coronavirus-annotation"
  }
  input {
    File genome_fasta
    String vadr_opts = "--noseqnamemax --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --out_allfasta"
    File vadr_model_file = "gs://theiagen-public-resources-rp/reference_data/databases/vadr_models/vadr-models-sarscov2-1.3-2.tar.gz"
    Int assembly_length_unambiguous
    Int skip_length = 10000
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.4"
    Int min_length = 50
    Int max_length = 30000
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
  }
  String out_base = basename(basename(basename(genome_fasta, ".fasta"), ".fa"), ".fna")
  command <<<
    set -euo pipefail

    if [ ~{assembly_length_unambiguous} -gt ~{skip_length} ]; then

      # extract the model file
      mkdir -p model_dir
      tar -C model_dir -xzf ~{vadr_model_file}

      # sometimes the model files are in a subdirectory and we need to find/move them.
      # the .minfo file is created by the v-build.pl command and is always in a valid model directory
      model_file_paths=$(find model_dir -type f -name "*.minfo")

      echo "DEBUG: Location(s) of '*.minfo' model files: "
      echo -e "${model_file_paths} \n"

      if [ -z "$model_file_paths" ]; then
        echo "ERROR: No model files found in the extracted model directory."
        exit 1
      fi

      # sometimes there can be multiple '*.minfo' files further nested in the model directory.
      # get the outermost (least nested) directory containing '*.minfo' model files.
      # then count the number of forward slashes and sort them to find the least nested path.
      top_model_file_path=$(echo "${model_file_paths}" | awk -F "/" '{print NF-1, $0}' | sort -n | head -n1 | cut -d' ' -f2)

      echo "DEBUG: Using least nested model file path: "
      echo -e "${top_model_file_path} \n"

      # get the directory containing the top-level model files
      top_model_dir=$(dirname "${top_model_file_path}")

      # vadr will expect the model files to be in this outermost directory
      mv "${top_model_dir}"/* model_dir/

      # remove any empty directories if they exist
      find model_dir -type d -empty -delete

      # remove terminal ambiguous nucleotides
      /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
        "~{genome_fasta}" \
        --minlen ~{min_length} \
        --maxlen ~{max_length} \
        > "~{out_base}_trimmed.fasta"

      # run VADR
      # --split and --cpu must be used in conjuction
      v-annotate.pl \
        --split --cpu ~{cpu} \
        --mdir model_dir \
        ~{vadr_opts} \
        "~{out_base}_trimmed.fasta" \
        "~{out_base}" \
        -v

      # package everything for output
      tar -C "~{out_base}" -czvf "~{out_base}.vadr.tar.gz" .

      # package up FASTA files into zip file for output. Note: this will work whether the --out_allfasta flag is included or not (there are just more when the option is used)
      mkdir -v vadr_fasta_files
      cp -v ~{out_base}/*.fa vadr_fasta_files
      zip ~{out_base}_vadr-fasta-files.zip vadr_fasta_files/*.fa 

      # prep alerts into a tsv file for parsing
      cut -f 5 "~{out_base}/~{out_base}.vadr.alt.list" | tail -n +2 > "~{out_base}.vadr.alerts.tsv"
      cat "~{out_base}.vadr.alerts.tsv" | wc -l > NUM_ALERTS

      # rename sequence classification summary file to end in txt
      mv -v "~{out_base}/~{out_base}.vadr.sqc" "~{out_base}/~{out_base}.vadr.sqc.txt"

    else
      echo "VADR skipped due to poor assembly; assembly length (unambiguous) = ~{assembly_length_unambiguous}" > NUM_ALERTS
    fi

  >>>
  output {
    File? feature_tbl_pass = "~{out_base}/~{out_base}.vadr.pass.tbl"
    File? feature_tbl_fail = "~{out_base}/~{out_base}.vadr.fail.tbl"
    File? classification_summary_file = "~{out_base}/~{out_base}.vadr.sqc.txt"
    String num_alerts = read_string("NUM_ALERTS")
    File? alerts_list = "~{out_base}/~{out_base}.vadr.alt.list"
    File? outputs_tgz = "~{out_base}.vadr.tar.gz"
    File? vadr_fastas_zip_archive = "~{out_base}_vadr-fasta-files.zip"
    String vadr_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    cpu: cpu
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}