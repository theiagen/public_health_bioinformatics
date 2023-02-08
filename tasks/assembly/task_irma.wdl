version 1.0

task irma {
  input {
    File read1
    File read2
    String samplename
    Boolean keep_ref_deletions = true
    String irma_module = "FLU"
    String read_basename = basename(read1)
    String docker = "quay.io/staphb/irma:1.0.3"
    Int memory = 8
    Int cpu = 4
    Int disk_size = 100
  }
  command <<<
    date | tee DATE
    #capture reads as bash variables
    read1=~{read1}
    read2=~{read2}
    # set cat command based on compression
    if [[ "~{read1}" == *".gz" ]] ; then
      cat_reads="zcat"
    else
      cat_reads="cat"
    fi
    # capture irma vesion
    IRMA | head -n1 | awk -F' ' '{ print "IRMA " $5 }' | tee VERSION
    # set config if needed
    if ~{keep_ref_deletions}; then 
      touch irma_config.sh
      echo 'DEL_TYPE="NNN"' >> irma_config.sh
      echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
    fi
    # format reads, if needed
    read_header=$(${cat_reads} ~{read1} | head -n1)
    if ! [[ "${read_header}" =~ @(.+?)[_[:space:]][123]:.+ ]]; then
      echo "Read headers may lead to IRMA failure; reformatting to meet IRMA input requirements"
      sra_id=$(echo "~{read_basename}" | awk -F "_" '{ print $1 }')
      eval "${cat_reads} ~{read1}" | awk '{print (NR%4 == 1) ? "@'${sra_id}'-" ++i " 1:1" : $0}' | gzip -c > "${sra_id}-irmafix_R1.fastq.gz"
      eval "${cat_reads} ~{read2}" | awk '{print (NR%4 == 1) ? "@'${sra_id}'-" ++i " 2:2" : $0}' | gzip -c > "${sra_id}-irmafix_R2.fastq.gz"
      #modify read variables
      read1="${sra_id}-irmafix_R1.fastq.gz"
      read2="${sra_id}-irmafix_R2.fastq.gz"
    else
      echo "Read headers match IRMA formatting requirements"
    fi
    # run IRMA 
    IRMA "~{irma_module}" "${read1}" "${read2}" ~{samplename}
    # capture IRMA type
    if compgen -G "~{samplename}/*fasta"; then
      echo "Type_"$(basename "$(echo "$(find ~{samplename}/*.fasta | head -n1)")" | cut -d_ -f1) > IRMA_TYPE
      # cat consensus assemblies
      cat ~{samplename}/*.fasta > ~{samplename}.irma.consensus.fasta
    else
      echo "No IRMA assembly generated for flu type prediction" >> IRMA_TYPE
    fi
    # rename IRMA outputs
    for irma_out in ~{samplename}/*{.vcf,.fasta,.bam}; do
      new_name="~{samplename}_"$(basename "${irma_out}" | cut -d "_" -f2- )
      echo "New name: ${new_name}; irma_out: ${irma_out}"
      mv "${irma_out}" "${new_name}"
    done
    # capture type A subtype
    if compgen -G "~{samplename}_HA*.fasta"; then # check if HA segment exists
      if [[ "$(ls ~{samplename}_HA*.fasta)" == *"HA_H"* ]]; then # if so, grab H-type if one is identified in assembly header
        subtype="$(basename ~{samplename}_HA*.fasta | awk -F _ '{print $NF}' | cut -d. -f1)" # grab H-type from last value in under-score-delimited filename
      fi
      # format HA segment to target output name
      mv "~{samplename}"_HA*.fasta "~{samplename}"_HA.fasta
    fi
    if compgen -G "~{samplename}_NA*.fasta" && [[ "$(ls ~{samplename}_NA*.fasta)" == *"NA_N"* ]]; then # check if NA segment exists with an N-type identified in header
       subtype+="$(basename ~{samplename}_NA*.fasta | awk -F _ '{print $NF}' | cut -d. -f1)" # grab N-type from last value in under-score-delimited filename
    fi
    if ! [ -z "${subtype}" ]; then 
      echo "${subtype}" > IRMA_SUBTYPE
    else
      echo "No subtype predicted by IRMA" > IRMA_SUBTYPE
    fi
  >>>
  output {
    File? irma_assembly_fasta = "~{samplename}.irma.consensus.fasta"
    File? seg4_ha_assembly = "~{samplename}_HA.fasta"
    String irma_type = read_string("IRMA_TYPE")
    String irma_subtype = read_string("IRMA_SUBTYPE")
    Array[File] irma_assemblies = glob("~{samplename}*.fasta")
    Array[File] irma_vcfs = glob("~{samplename}*.vcf")
    Array[File] irma_bams = glob("~{samplename}*.bam")
    String irma_version = read_string("VERSION")
    String irma_pipeline_date = read_string("DATE")
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