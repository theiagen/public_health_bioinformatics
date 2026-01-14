version 1.0

task irma {
  input {
    File read1
    File? read2
    String seq_method
    String samplename
    Boolean keep_ref_deletions = true # set DEL_TYPE config in irma_config.sh to "DEL" if false, "NNN" if true
    Int minimum_consensus_support = 50 # IRMA default is 1, but matching MIRA standards for ONT = 50 and ILMN = 30 via defaults at theiacov workflow level WDLs: https://cdcgov.github.io/MIRA/articles/sequence-qc.html
    Int minimum_read_length = 75 # matching default for TheiaCoV_Illumina_PE; NOTE: IRMA's default is 125 bp
    Int minimum_average_consensus_allele_quality = 10 # IRMA default is 0, we are matching MIRA standards for both ONT and ILMN: https://cdcgov.github.io/MIRA/articles/sequence-qc.html
    Float minimum_ambiguous_threshold = 0.20
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/irma:1.3.1-dev"
    Int memory = 16
    Int cpu = 4
    Int disk_size = 100
  }
  command <<<

    # capture irma vesion
    IRMA | head -n1 | awk -F' ' '{ print "IRMA " $5 }' | tee VERSION

    # set -euo pipefail to avoid silent failure; must happen AFTER running "IRMA" since it throws exit code 1
    set -euo pipefail

    ### IRMA configuration ###
  
    # set how to handle deletions
    if ~{keep_ref_deletions}; then 
      export DEL_TYPE="NNN"
    else # default in WDL and IRMA
      # IRMA docs state: If sites are completely missing during read gathering use the reference seed (REF), delete by ambiguation (NNN), or just remove (DEL)
      export DEL_TYPE="DEL"
    fi

    python3 <<CODE
    import os 
    # CPU config
    num_cpus_actual = os.cpu_count()
    del_type = os.environ.get("DEL_TYPE")
    with open("irma_config.sh", "w") as conf_file:
      lines = [f"SINGLE_LOCAL_PROC={num_cpus_actual}\n",
              f"DOUBLE_LOCAL_PROC={num_cpus_actual // 2}\n",
              "MIN_CONS_SUPPORT=~{minimum_consensus_support}\n",
              "MIN_CONS_QUALITY=~{minimum_average_consensus_allele_quality}\n",
              "MIN_AMBIG=~{minimum_ambiguous_threshold}\n",
              "TMP=$(pwd)\n",
              f"DEL_TYPE='{del_type}'\n",
              "ALIGN_PROG='BLAT'\n",
              "MIN_LEN=~{minimum_read_length}"]
      conf_file.writelines(lines)
      conf_file.close()
    CODE

    # run IRMA
    # set IRMA module depending on sequencing technology
    if [[ ~{seq_method} == "OXFORD_NANOPORE" ]]; then
      IRMA "FLU-minion" "~{read1}" ~{samplename} --external-config irma_config.sh
    else
      # else, assume Illumina paired-end data as input
      IRMA "FLU" "~{read1}" "~{read2}" ~{samplename} --external-config irma_config.sh
    fi

    if compgen -G "~{samplename}/*fasta"; then
      # capture some IRMA log & config files; rename to use .tsv suffix instead of .txt
      mv -v ~{samplename}/tables/READ_COUNTS.txt ~{samplename}/tables/READ_COUNTS.tsv
      mv -v ~{samplename}/logs/run_info.txt ~{samplename}/logs/run_info.tsv

      # look at list of files that match the above pattern, grab the first one, and extract the type from the filename. We expect: ~{samplename}/B_HA.fasta
      irma_type="$(find ~{samplename}/*.fasta -type f | head -n1 | xargs -n1 basename | cut -d_ -f1)"
      echo "Type_$irma_type" > IRMA_TYPE

      # Helper script to handle the creation of QC summary, renaming and creation of secondary FASTA files, and acquisition of subtype
      python3 /scripts/irma_helper.py -d ~{samplename} -t "${irma_type}" -s ~{samplename} -v

      # rename IRMA outputs to include samplename. Example: "B_HA.fasta" -> "sample0001_HA.fasta"
      echo "DEBUG: Renaming IRMA output VCFs, FASTAs, and BAMs to include samplename...."
      for irma_out in ~{samplename}/*{.vcf,.fasta,.bam}; do
        new_name="~{samplename}_"$(basename "${irma_out}" | cut -d "_" -f2- )
        mv -v "${irma_out}" "${new_name}"
      done

      if ls "~{samplename}"_HA?*.bam 1> /dev/null 2>&1; then
        echo "DEBUG: Renaming HA BAM files...."
        for file in "~{samplename}"_HA?*.bam; do
          mv -v "$file" "${file%_HA*.bam}_HA.bam"
        done
      fi
      if ls "~{samplename}"_NA?*.bam 1> /dev/null 2>&1; then
        echo "DEBUG: Renaming NA BAM files...."
        for file in "~{samplename}"_NA?*.bam; do
          mv -v "$file" "${file%_NA*.bam}_NA.bam"
        done
      fi

    else
      echo "No IRMA assembly generated for flu type prediction" | tee IRMA_TYPE
      echo "No subtype predicted by IRMA" > IRMA_SUBTYPE.txt
      echo "No subtype notes" > IRMA_SUBTYPE_NOTES.txt
      echo "Exiting IRMA task early since no IRMA assembly was generated."
      exit 0
    fi
  >>>
  output {
    # all of these FASTAs are derived from the amended_consensus/*.fa files produced by IRMA
    File? irma_assembly_fasta = "~{samplename}/amended_consensus/~{samplename}.irma.consensus.fasta"
    File? irma_assembly_fasta_concatenated = "~{samplename}/amended_consensus/~{samplename}.irma.consensus.concatenated.fasta"
    File? seg_ha_assembly = "~{samplename}/amended_consensus/~{samplename}_HA.fasta"
    File? seg_na_assembly = "~{samplename}/amended_consensus/~{samplename}_NA.fasta"
    File? seg_pa_assembly = "~{samplename}/amended_consensus/~{samplename}_PA.fasta"
    File? seg_pb1_assembly = "~{samplename}/amended_consensus/~{samplename}_PB1.fasta"
    File? seg_pb2_assembly = "~{samplename}/amended_consensus/~{samplename}_PB2.fasta"
    File? seg_mp_assembly = "~{samplename}/amended_consensus/~{samplename}_MP.fasta"
    File? seg_np_assembly = "~{samplename}/amended_consensus/~{samplename}_NP.fasta"
    File? seg_ns_assembly = "~{samplename}/amended_consensus/~{samplename}_NS.fasta"
    
    # adding these "padded" assemblies as outputs to be passed to VADR and MAFFT (antiviral substitutions tasks)
    # we may remove these outputs in the future if IRMA code is updated to not output periods in the consensus sequences
    File? irma_assembly_fasta_padded = "padded_assemblies/~{samplename}.irma.consensus.pad.fasta"
    File? irma_assembly_fasta_concatenated_padded = "padded_assemblies/~{samplename}.irma.consensus.concatenated.pad.fasta"
    File? seg_ha_assembly_padded = "padded_assemblies/~{samplename}_HA.pad.fasta"
    File? seg_na_assembly_padded = "padded_assemblies/~{samplename}_NA.pad.fasta"
    File? seg_pa_assembly_padded = "padded_assemblies/~{samplename}_PA.pad.fasta"
    File? seg_pb1_assembly_padded = "padded_assemblies/~{samplename}_PB1.pad.fasta"
    File? seg_pb2_assembly_padded = "padded_assemblies/~{samplename}_PB2.pad.fasta"
    File? seg_mp_assembly_padded = "padded_assemblies/~{samplename}_MP.pad.fasta"
    File? seg_np_assembly_padded = "padded_assemblies/~{samplename}_NP.pad.fasta"
    File? seg_ns_assembly_padded = "padded_assemblies/~{samplename}_NS.pad.fasta"

    # Output MIRA-like QC summary and all variants called by IRMA for each segment.
    File? irma_qc_summary_tsv = "~{samplename}/~{samplename}_irma_qc_summary.tsv"
    File? irma_all_snvs_tsv = "~{samplename}/tables/~{samplename}_irma_all_variants.tsv"
    File? irma_all_insertions_tsv = "~{samplename}/tables/~{samplename}_irma_all_insertions.tsv"
    File? irma_all_deletions_tsv = "~{samplename}/tables/~{samplename}_irma_all_deletions.tsv"

    String irma_type = read_string("IRMA_TYPE")
    String irma_subtype = read_string("IRMA_SUBTYPE.txt")
    String irma_subtype_notes = read_string("IRMA_SUBTYPE_NOTES.txt")
    Array[File] irma_plurality_consensus_assemblies = glob("~{samplename}*.fasta")
    Array[File] irma_vcfs = glob("~{samplename}*.vcf")
    Array[File] irma_bams = glob("~{samplename}*.bam")
    String irma_docker = docker
    String irma_version = read_string("VERSION")
    # tracking this for outputting min depth threshold used for calling consensus nucleotides in IRMA
    Int irma_minimum_consensus_support = minimum_consensus_support
    File? irma_read_counts_tsv = "~{samplename}/tables/READ_COUNTS.tsv"
    File? irma_run_info_tsv = "~{samplename}/logs/run_info.tsv"
    File? irma_nr_read_counts = "~{samplename}/logs/NR_COUNTS_log.txt"
    # for now just adding bams for these segments for mean coverage calculation
    File? seg_ha_bam = "~{samplename}_HA.bam"
    File? seg_na_bam = "~{samplename}_NA.bam"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible:  0
  }
}