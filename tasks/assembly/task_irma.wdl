version 1.0

task irma {
  input {
    File read1
    File? read2
    String seq_method
    String samplename
    Boolean keep_ref_deletions = true
    String read_basename = basename(read1)
    String docker = "us-docker.pkg.dev/general-theiagen/cdcgov/irma:v1.1.5"
    Int memory = 16
    Int cpu = 4
    Int disk_size = 100
  }
  command <<<
    date | tee DATE

    num_cpus_actual=$(nproc)
    echo "DEBUG: Number of CPUs available: ${num_cpus_actual}. Setting this in irma_config.sh file..."
    echo "SINGLE_LOCAL_PROC=${num_cpus_actual}" > irma_config.sh
    # set this variable to half the value of num_cpus_actual, as per the IRMA documentation: https://wonder.cdc.gov/amd/flu/irma/configuration.html
    echo "DOUBLE_LOCAL_PROC=$((${num_cpus_actual}/2))" >> irma_config.sh

    # this is done so that IRMA used PWD as the TMP directory instead of /tmp/root that it tries by default; cromwell doesn't allocate much disk space here (64MB or some small amount)
    echo "DEBUG: creating an optional IRMA configuration file to set TMP directory to $(pwd)"
    echo "TMP=$(pwd)" >> irma_config.sh

    #capture reads as bash variables
    read1=~{read1}
    if [[ "~{read2}" ]]; then 
      read2=~{read2}
    fi

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
      # attempting this to try to use Ns in final output FASTAs instead of periods; output assembly should have .pad.fa suffix
      # TODO test again and look at .pad.fa files
      #echo "ALIGN_AMENDED=1" >> irma_config.sh
      #echo "ASSEM_REF=1" >> irma_config.sh
      #echo "PADDED_CONSENSUS=1" >> irma_config.sh
    fi

    # format reads, if needed
    read_header=$(${cat_reads} ~{read1} | head -n1)
    if ! [[ "${read_header}" =~ @(.+?)[_[:space:]][123]:.+ ]]; then
      echo "Read headers may lead to IRMA failure; reformatting to meet IRMA input requirements"
      sra_id=$(echo "~{read_basename}" | awk -F "_" '{ print $1 }')
      eval "${cat_reads} ~{read1}" | awk '{print (NR%4 == 1) ? "@'${sra_id}'-" ++i " 1:1" : $0}' | gzip -c > "${sra_id}-irmafix_R1.fastq.gz"
      read1="${sra_id}-irmafix_R1.fastq.gz"
      if [[ "~{read2}" ]]; then 
        eval "${cat_reads} ~{read2}" | awk '{print (NR%4 == 1) ? "@'${sra_id}'-" ++i " 2:2" : $0}' | gzip -c > "${sra_id}-irmafix_R2.fastq.gz"
        read2="${sra_id}-irmafix_R2.fastq.gz"
      fi     
    else
      echo "Read headers match IRMA formatting requirements"
    fi

    echo "DEBUG: Custom irma_config.sh file contents:"
    cat irma_config.sh

    # run IRMA
    # set IRMA module depending on sequencing technology
    if [[ ~{seq_method} == "OXFORD_NANOPORE" ]]; then
      IRMA "FLU-minion" "${read1}" ~{samplename} --external-config irma_config.sh
    else
      # else, assume Illumina paired-end data as input
      IRMA "FLU" "${read1}" "${read2}" ~{samplename} --external-config irma_config.sh
    fi

    # capture IRMA type
    if compgen -G "~{samplename}/*fasta"; then
      # look at list of files that match the above pattern, grab the first one, and extract the type from the filename. We expect: ~{samplename}/B_HA.fasta
      echo "Type_"$(basename "$(echo "$(find ~{samplename}/*.fasta | head -n1)")" | cut -d_ -f1) > IRMA_TYPE
      # set irma_type bash variable which is used later
      irma_type=$(cat IRMA_TYPE)
      # concatenate consensus assemblies into single file with all genome segments
      echo "DEBUG: creating IRMA FASTA file containing all segments...."
      cat ~{samplename}/*.fasta > ~{samplename}.irma.consensus.fasta
      echo "DEBUG: editing IRMA FASTA file to include sample name in FASTA headers...."
      sed -i "s/>/>~{samplename}_/g" ~{samplename}.irma.consensus.fasta

      echo "DEBUG: creating copy of consensus FASTA with periods replaced by Ns...."
      # use sed to create copy of FASTA file where periods are replaced by Ns, except in the FASTA header lines that begin with '>'
      sed '/^>/! s/\./N/g' ~{samplename}.irma.consensus.fasta > ~{samplename}.irma.consensus.pad.fasta
    else
      echo "No IRMA assembly generated for flu type prediction" | tee IRMA_TYPE
    fi

    # rename IRMA outputs to include samplename. Example: "B_HA.fasta" -> "sample0001_HA.fasta"
    echo "DEBUG: Renaming IRMA output VCFs, FASTAs, and BAMs to include samplename...."
    for irma_out in ~{samplename}/*{.vcf,.fasta,.bam}; do
      new_name="~{samplename}_"$(basename "${irma_out}" | cut -d "_" -f2- )
      mv -v "${irma_out}" "${new_name}"
    done
    
    # capture type A subtype
    if compgen -G "~{samplename}_HA*.fasta"; then # check if HA segment exists
      # NOTE: this if block does not get triggered for Flu B samples, because they do not include a subtype in FASTA filename for HA or NA segment
      if [[ "$(ls ~{samplename}_HA*.fasta)" == *"HA_H"* ]]; then # if so, grab H-type if one is identified in assembly header
        subtype="$(basename ~{samplename}_HA*.fasta | awk -F _ '{print $NF}' | cut -d. -f1)" # grab H-type from last value in under-score-delimited filename
        # rename HA FASTA file to not include subtype in filename)
        echo "DEBUG: renaming HA FASTA file to not include subtype in filename...."
        mv -v ~{samplename}_HA*.fasta ~{samplename}_HA.fasta
      fi
      echo "DEBUG: Running sed to change HA segment FASTA header now..."
      # format HA segment to target output name and rename header to include the samplename. Example FASTA header change: ">B_HA" -> ">sample0001_B_HA"
      sed -i "1s/>/>~{samplename}_/" "~{samplename}_HA.fasta"
    fi

    # if there is a file that matches the pattern AND the file contains an N-type in the header, grab the N-type
    # NOTE: this does not get triggered for Flu B samples, because they do not include a subtype in FASTA filename for HA or NA segment
    if compgen -G "~{samplename}_NA*.fasta" && [[ "$(ls ~{samplename}_NA*.fasta)" == *"NA_N"* ]]; then # check if NA segment exists with an N-type identified in header
      subtype+="$(basename ~{samplename}_NA*.fasta | awk -F _ '{print $NF}' | cut -d. -f1)" # grab N-type from last value in under-score-delimited filename 
      # rename NA FASTA file to not include subtype in filename)
      echo "DEBUG: renaming NA FASTA file to not include subtype in filename...."
      mv -v ~{samplename}_NA*.fasta ~{samplename}_NA.fasta
    fi

    echo "DEBUG: Running sed to change NA FASTA header now..."
    # format NA segment to target output name and rename header to include the samplename
    sed -i "1s/>/>~{samplename}_/" "~{samplename}_NA.fasta"

    # if bash variable "subtype" is not empty, write it to a file; 
    # otherwise, write a message indicating no subtype was predicted
    if [ -n "${subtype}" ]; then 
      echo "${subtype}" | tee IRMA_SUBTYPE
    else
      echo "No subtype predicted by IRMA" | tee IRMA_SUBTYPE
    fi

    # if "subtype" is "Type_B" then write a note indicating that IRMA does not differentiate between Victoria and Yamagata lineages
    if [[ "${irma_type}" == "Type_B" ]]; then
      echo "IRMA does not differentiate Victoria and Yamagata Flu B lineages. See abricate_flu_subtype output column" | tee IRMA_SUBTYPE_NOTES
    else
      # create empty file
      touch IRMA_SUBTYPE_NOTES
    fi

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
    
    echo "DEBUG: Creating padded FASTA files for each individual segment...."
    # this loop looks in the PWD for files ending in .fasta, and creates a copy of the file with periods replaced by Ns and dashes are deleted (since they represent gaps)
    for FASTA in ./*.fasta; do
      # if the renamed file ends in .fasta; then create copy of the file with periods replaced by Ns and dashes removed in all lines except FASTA headers that begin with '>'
        echo "DEBUG: creating padded FASTA file for ${FASTA}...."
        sed '/^>/! s/\./N/g' "${FASTA}" > "${FASTA%.fasta}.temp.fasta"
        sed '/^>/! s/-//g' "${FASTA%.fasta}.temp.fasta" > "${FASTA%.fasta}.pad.fasta"
    done
  >>>
  output {
    File? irma_assembly_fasta = "~{samplename}.irma.consensus.fasta"
    File? seg_ha_assembly = "~{samplename}_HA.fasta"
    File? seg_na_assembly = "~{samplename}_NA.fasta"
    File? seg_pa_assembly = "~{samplename}_PA.fasta"
    File? seg_pb1_assembly = "~{samplename}_PB1.fasta"
    File? seg_pb2_assembly = "~{samplename}_PB2.fasta"
    File? seg_mp_assembly = "~{samplename}_MP.fasta"
    File? seg_np_assembly = "~{samplename}_NP.fasta"
    File? seg_ns_assembly = "~{samplename}_NS.fasta"
    # adding these "padded" assemblies as outputs to be passed to VADR and MAFFT (antiviral substitutions tasks)
    # we may remove these outputs in the future if IRMA code is updated to not output periods in the consensus sequences
    File? irma_assembly_fasta_padded = "~{samplename}.irma.consensus.pad.fasta"
    File? seg_ha_assembly_padded = "~{samplename}_HA.pad.fasta"
    File? seg_na_assembly_padded = "~{samplename}_NA.pad.fasta"
    File? seg_pa_assembly_padded = "~{samplename}_PA.pad.fasta"
    File? seg_pb1_assembly_padded = "~{samplename}_PB1.pad.fasta"
    File? seg_pb2_assembly_padded = "~{samplename}_PB2.pad.fasta"
    File? seg_mp_assembly_padded = "~{samplename}_MP.pad.fasta"
    File? seg_np_assembly_padded = "~{samplename}_NP.pad.fasta"
    File? seg_ns_assembly_padded = "~{samplename}_NS.pad.fasta"

    String irma_type = read_string("IRMA_TYPE")
    String irma_subtype = read_string("IRMA_SUBTYPE")
    String irma_subtype_notes = read_string("IRMA_SUBTYPE_NOTES")
    Array[File] irma_assemblies = glob("~{samplename}*.fasta")
    Array[File] irma_vcfs = glob("~{samplename}*.vcf")
    Array[File] irma_bams = glob("~{samplename}*.bam")
    String irma_docker = docker
    String irma_version = read_string("VERSION")
    String irma_pipeline_date = read_string("DATE")
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