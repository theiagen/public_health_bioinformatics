version 1.0

task ncbi_prep_one_sample {
  input {
    #required files
    File assembly_fasta
    File read1_dehosted
    File read2_dehosted
    #required metadata
    String assembly_method
    String bioproject_accession
    String collecting_lab
    String collection_date
    String country
    String design_description
    String dehosting_method
    String filetype
    String host_disease
    String host
    String host_sci_name
    String instrument_model
    String isolation_source
    String library_id
    String library_layout
    String library_selection
    String library_source
    String library_strategy
    String organism
    String seq_platform
    String state
    String submission_id
    #optional metadata
    String? amplicon_primer_scheme
    String? amplicon_size
    String? biosample_accession = "{populate_with_bioSample_accession}"
    String? gisaid_accession
    String? gisaid_organism="hCoV-2019"
    String? patient_age
    String? patient_gender
    String? purpose_of_sampling
    String? purpose_of_sequencing
    String? submitter_email
    String? treatment
    #GenBank formatting
    Int minlen = 50
    Int maxlen = 30000
    #runtime
    String docker_image = "quay.io/staphb/vadr:1.3"
    Int  memory = 1
    Int cpu = 1
    Int disk_size = 25
    Int preemptible_tries = 0
  }
  command <<<
    #Check date format
    if [[ ~{collection_date} != [0-9][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9] ]]
    then 
      echo "Incorrect collection date format; collection date must be in  YYYY-MM-DD format." 1>&2
      exit 1
    else
      year=$(echo ~{collection_date} | cut -f 1 -d '-')
    fi
    
    # capture sample variables
    isolate="~{organism}/~{host}/~{country}/~{submission_id}/${year}"
    gisaid_virus_name="~{gisaid_organism}/~{country}/~{submission_id}/$year"
    
    #Format BioSample Attributes
    echo -e "*sample_name\tsample_title\tbioproject_accession\t*organism\t*collected_by\t*collection_date\t*geo_loc_name\t*host\t*host_disease\t*isolate\t*isolation_source\tantiviral_treatment_agent\tcollection_device\tcollection_method\tdate_of_prior_antiviral_treat\tdate_of_prior_sars_cov_2_infection\tdate_of_sars_cov_2_vaccination\texposure_event\tgeo_loc_exposure\tgisaid_accession\tgisaid_virus_name\thost_age\thost_anatomical_material\thost_anatomical_part\thost_body_product\thost_disease_outcome\thost_health_state\thost_recent_travel_loc\thost_recent_travel_return_date\thost_sex\thost_specimen_voucher\thost_subject_id\tlat_lon\tpassage_method\tpassage_number\tprior_sars_cov_2_antiviral_treat\tprior_sars_cov_2_infection\tprior_sars_cov_2_vaccination\tpurpose_of_sampling\tpurpose_of_sequencing\tsars_cov_2_diag_gene_name_1\tsars_cov_2_diag_gene_name_2\tsars_cov_2_diag_pcr_ct_value_1\tsars_cov_2_diag_pcr_ct_value_2\tsequenced_by\tvaccine_received\tvirus_isolate_of_prior_infection\tdescription" > ~{submission_id}_biosample_attributes.tsv
    
    echo -e "~{submission_id}\t\t~{bioproject_accession}\t~{organism}\t~{collecting_lab}\t~{collection_date}\t~{country}: ~{state}\t~{host_sci_name}\t~{host_disease}\t${isolate}\t~{isolation_source}\t~{treatment}\t\t\t\t\t\t\t\t~{gisaid_accession}\t${gisaid_virus_name}\t~{patient_age}\t\t\t\t\t\t\t\t~{patient_gender}\t\t\t\t\t\t\t\t\t~{purpose_of_sampling}\t~{purpose_of_sequencing}\t\t\t\t\t\t\t\t" >> ~{submission_id}_biosample_attributes.tsv
    
    #Format SRA Reads & Metadata
    cp ~{read1_dehosted} ~{submission_id}_R1.fastq.gz
    cp ~{read2_dehosted} ~{submission_id}_R2.fastq.gz

    echo -e "bioproject_accession\tsample_name\tlibrary_ID\ttitle\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout\tplatform\tinstrument_model\tdesign_description\tfiletype\tfilename\tfilename2\tfilename3\tfilename4\tamplicon _PCR_primer_scheme\tamplicon_size\tsequencing_protocol_name\traw_sequence_data_processing_method \tdehosting_method\tsequence_submitter_contact_email" > ~{submission_id}_sra_metadata.tsv
    
    echo -e "~{bioproject_accession}\t~{submission_id}\t~{library_id}\tGenomic sequencing of ~{organism}: ~{isolation_source}\t~{library_strategy}\t~{library_source}\t~{library_selection}\t~{library_layout}\t~{seq_platform}\t~{instrument_model}\t~{design_description}\t~{filetype}\t~{submission_id}_R1.fastq.gz\t~{submission_id}_R2.fastq.gz\t\t\t~{amplicon_primer_scheme}\t~{amplicon_size}\t\t~{assembly_method}\t~{dehosting_method}\t~{submitter_email}" >> ~{submission_id}_sra_metadata.tsv 
     
    #Format GenBank metadata and assembly    
    ##GenBank assembly
    ###Adding metadata for genbank submissions and removing line breaks
    echo ">"~{submission_id} > ~{submission_id}_genbank_untrimmed.fasta
    grep -v ">" ~{assembly_fasta} | tr -d '\n' >> ~{submission_id}_genbank_untrimmed.fasta
    
    # remove terminal ambiguous nucleotides
    /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
      "~{submission_id}_genbank_untrimmed.fasta" \
      --minlen ~{minlen} \
      --maxlen ~{maxlen} \
      > "~{submission_id}_genbank.fasta"
    
    ##GenBank modifier
    echo -e "Sequence_ID\tcountry\thost\tisolate\tcollection-date\tisolation-source\tBioSample\tBioProject\tnote" > ~{submission_id}_genbank_modifier.tsv
    echo -e "~{submission_id}\t~{country}\t~{host_sci_name}\t${isolate}\t~{collection_date}\t~{isolation_source}\t~{biosample_accession}\t~{bioproject_accession}"  >> ~{submission_id}_genbank_modifier.tsv
    
  >>>
  output {
    File biosample_attributes = "~{submission_id}_biosample_attributes.tsv"
    File sra_metadata = "~{submission_id}_sra_metadata.tsv"
    File genbank_assembly = "~{submission_id}_genbank.fasta"
    File genbank_modifier = "~{submission_id}_genbank_modifier.tsv"
    File sra_read1 = "~{submission_id}_R1.fastq.gz"
    File sra_read2 = "~{submission_id}_R2.fastq.gz"
    Array[File] sra_reads = ["~{submission_id}_R1.fastq.gz","~{submission_id}_R2.fastq.gz"]
  }
  runtime {
    docker: "~{docker_image}"
    memory: "~{memory} GB"
    cpu: cpu    
    disks: "local-disk ~{disk_size} SSD"
    disk: disk_size + " GB" # TES
    preemptible: preemptible_tries
    maxRetries: 3
  }
}

task ncbi_prep_one_sample_se {
  input {
    #required files
    File assembly_fasta
    File reads_dehosted
    #required metadata
    String assembly_method
    String bioproject_accession
    String collecting_lab
    String collection_date
    String country
    String design_description
    String dehosting_method
    String filetype
    String host_disease
    String host
    String host_sci_name
    String instrument_model
    String isolation_source
    String library_id
    String library_layout
    String library_selection
    String library_source
    String library_strategy
    String organism
    String seq_platform
    String state
    String submission_id
    #optional metadata
    String? amplicon_primer_scheme
    String? amplicon_size
    String? biosample_accession = "{populate_with_bioSample_accession}"
    String? gisaid_accession
    String? gisaid_organism="hCoV-2019"
    String? patient_age
    String? patient_gender
    String? purpose_of_sampling
    String? purpose_of_sequencing
    String? submitter_email
    String? treatment
    #GenBank formatting
    Int minlen = 50
    Int maxlen = 30000
    #runtime
    String docker = "quay.io/staphb/vadr:1.3"
    Int  memory = 1
    Int cpu = 1
    Int disk_size = 25
    Int preemptible_tries = 0
  }
  command <<<
    #Check date format
    if [[ ~{collection_date} != [0-9][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9] ]]
    then 
      echo "Incorrect collection date format; collection date must be in  YYYY-MM-DD format." 1>&2
      exit 1
    else
      year=$(echo ~{collection_date} | cut -f 1 -d '-')
    fi
    
    # capture sample variables
    isolate="~{organism}/~{host}/~{country}/~{submission_id}/${year}"
    gisaid_virus_name="~{gisaid_organism}/~{country}/~{submission_id}/$year"
    
    #Format BioSample Attributes
    echo -e "*sample_name\tsample_title\tbioproject_accession\t*organism\t*collected_by\t*collection_date\t*geo_loc_name\t*host\t*host_disease\t*isolate\t*isolation_source\tantiviral_treatment_agent\tcollection_device\tcollection_method\tdate_of_prior_antiviral_treat\tdate_of_prior_sars_cov_2_infection\tdate_of_sars_cov_2_vaccination\texposure_event\tgeo_loc_exposure\tgisaid_accession\tgisaid_virus_name\thost_age\thost_anatomical_material\thost_anatomical_part\thost_body_product\thost_disease_outcome\thost_health_state\thost_recent_travel_loc\thost_recent_travel_return_date\thost_sex\thost_specimen_voucher\thost_subject_id\tlat_lon\tpassage_method\tpassage_number\tprior_sars_cov_2_antiviral_treat\tprior_sars_cov_2_infection\tprior_sars_cov_2_vaccination\tpurpose_of_sampling\tpurpose_of_sequencing\tsars_cov_2_diag_gene_name_1\tsars_cov_2_diag_gene_name_2\tsars_cov_2_diag_pcr_ct_value_1\tsars_cov_2_diag_pcr_ct_value_2\tsequenced_by\tvaccine_received\tvirus_isolate_of_prior_infection\tdescription" > ~{submission_id}_biosample_attributes.tsv
    
    echo -e "~{submission_id}\t\t~{bioproject_accession}\t~{organism}\t~{collecting_lab}\t~{collection_date}\t~{country}: ~{state}\t~{host_sci_name}\t~{host_disease}\t${isolate}\t~{isolation_source}\t~{treatment}\t\t\t\t\t\t\t\t~{gisaid_accession}\t${gisaid_virus_name}\t~{patient_age}\t\t\t\t\t\t\t\t~{patient_gender}\t\t\t\t\t\t\t\t\t~{purpose_of_sampling}\t~{purpose_of_sequencing}\t\t\t\t\t\t\t\t" >> ~{submission_id}_biosample_attributes.tsv
    
    #Format SRA Reads & Metadata
    cp ~{reads_dehosted} ~{submission_id}_R1.fastq.gz

    echo -e "bioproject_accession\tsample_name\tlibrary_ID\ttitle\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout\tplatform\tinstrument_model\tdesign_description\tfiletype\tfilename\tfilename2\tfilename3\tfilename4\tamplicon _PCR_primer_scheme\tamplicon_size\tsequencing_protocol_name\traw_sequence_data_processing_method \tdehosting_method\tsequence_submitter_contact_email" > ~{submission_id}_sra_metadata.tsv
    
    echo -e "~{bioproject_accession}\t~{submission_id}\t~{library_id}\tGenomic sequencing of ~{organism}: ~{isolation_source}\t~{library_strategy}\t~{library_source}\t~{library_selection}\t~{library_layout}\t~{seq_platform}\t~{instrument_model}\t~{design_description}\t~{filetype}\t~{submission_id}_R1.fastq.gz\t\t\t\t~{amplicon_primer_scheme}\t~{amplicon_size}\t\t~{assembly_method}\t~{dehosting_method}\t~{submitter_email}" >> ~{submission_id}_sra_metadata.tsv 
     
    #Format GenBank metadata and assembly    
    ##GenBank assembly
    ###Adding metadata for genbank submissions and removing line breaks
    echo ">"~{submission_id} > ~{submission_id}_genbank_untrimmed.fasta
    grep -v ">" ~{assembly_fasta} | tr -d '\n' >> ~{submission_id}_genbank_untrimmed.fasta

    # remove terminal ambiguous nucleotides
    /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
      "~{submission_id}_genbank_untrimmed.fasta" \
      --minlen ~{minlen} \
      --maxlen ~{maxlen} \
      > "~{submission_id}_genbank.fasta"
      
    ##GenBank modifier
    echo -e "Sequence_ID\tcountry\thost\tisolate\tcollection-date\tisolation-source\tBioSample\tBioProject\tnote" > ~{submission_id}_genbank_modifier.tsv
    echo -e "~{submission_id}\t~{country}\t~{host_sci_name}\t${isolate}\t~{collection_date}\t~{isolation_source}\t~{biosample_accession}\t~{bioproject_accession}"  >> ~{submission_id}_genbank_modifier.tsv
  >>>
  output {
    File biosample_attributes = "~{submission_id}_biosample_attributes.tsv"
    File sra_metadata = "~{submission_id}_sra_metadata.tsv"
    File genbank_assembly = "~{submission_id}_genbank.fasta"
    File genbank_modifier = "~{submission_id}_genbank_modifier.tsv"
    File sra_reads = "~{submission_id}_R1.fastq.gz"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    disk: disk_size + " GB" # TES
    preemptible: preemptible_tries
    maxRetries: 3
  }
}

task gisaid_prep_one_sample {
  input {
    #required files
    File assembly_fasta
    #required metadata
    String authors
    String assembly_method
    String collecting_lab
    String collecting_lab_address
    String collection_date
    String continent
    Float assembly_mean_coverage
    String country
    String gisaid_submitter
    String host
    String organism = "hCoV-19"
    String seq_platform
    String state
    String submission_id
    String submitting_lab
    String submitting_lab_address
    String type="betacoronavirus"
    #optional metadata
    String? county
    String? patient_gender = "unknown"
    String? last_vaccinated
    String? passage_details = "original"
    String? patient_age = "unknown"
    String? patient_status = "unknown"
    String? purpose_of_sequencing
    String? outbreak
    String? specimen_source
    String? treatment
    #runtime
    String docker = "quay.io/theiagen/utility:1.1"
    Int  memory = 1
    Int cpu = 1
    Int disk_size = 25
    Int preemptible = 0
  }
  command <<<
    #Check date format
    if [[ ~{collection_date} != [0-9][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9] ]]
    then 
      echo "Incorrect collection date format; collection date must be in  YYYY-MM-DD format." 1>&2
      exit 1
    else
      year=$(echo ~{collection_date} | cut -f 1 -d '-')
    fi 
    
    #Format GISAID metadata and assembly
    ##GISAID assembly
    gisaid_virus_name="~{organism}/~{country}/~{submission_id}/$year"   
    echo ">"${gisaid_virus_name} > ~{submission_id}_gisaid.fasta
    grep -v ">" ~{assembly_fasta} >> ~{submission_id}_gisaid.fasta
    
    ##GISAID tMetadata
    echo "submitter,fn,covv_virus_name,covv_type,covv_passage,covv_collection_date,covv_location,covv_add_location,covv_host,covv_add_host_info,covv_sampling_strategy,covv_gender,covv_patient_age,covv_patient_status,covv_specimen,covv_outbreak,covv_last_vaccinated,covv_treatment,covv_seq_technology,covv_assembly_method,covv_coverage,covv_orig_lab,covv_orig_lab_addr,covv_provider_sample_id,covv_subm_lab,covv_subm_lab_addr,covv_subm_sample_id,covv_authors,covv_comment,comment_type" >  ~{submission_id}_gisaid_metadata.csv

    echo "Submitter,FASTA filename,Virus name,Type,Passage details/history,Collection date,Location,Additional location information,Host,Additional host information,Sampling Strategy,Gender,Patient age,Patient status,Specimen source,Outbreak,Last vaccinated,Treatment,Sequencing technology,Assembly method,Coverage,Originating lab,Address,Sample ID given by the sample provider,Submitting lab,Address,Sample ID given by the submitting laboratory,Authors,Comment,Comment Icon" >> ~{submission_id}_gisaid_metadata.csv

    echo "\"~{gisaid_submitter}\",\"~{submission_id}.gisaid.fa\",\"~{organism}/~{country}/~{submission_id}/$year\",\"~{type}\",\"~{passage_details}\",\"~{collection_date}\",\"~{continent}/~{country}/~{state}/~{county}\",,\"~{host}\",,\"~{purpose_of_sequencing}\",\"~{patient_gender}\",\"~{patient_age}\",\"~{patient_status}\",\"~{specimen_source}\",\"~{outbreak}\",\"~{last_vaccinated}\",\"~{treatment}\",\"~{seq_platform}\",\"~{assembly_method}\",\"~{assembly_mean_coverage}\",\"~{collecting_lab}\",\"~{collecting_lab_address}\",,\"~{submitting_lab}\",\"~{submitting_lab_address}\",,\"~{authors}\",," >> ~{submission_id}_gisaid_metadata.csv
  >>>
  output {
    File gisaid_assembly = "~{submission_id}_gisaid.fasta"
    File gisaid_metadata = "~{submission_id}_gisaid_metadata.csv"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    disk: disk_size + " GB" # TES
    preemptible: preemptible
    maxRetries:   3
  }
}

task compile_assembly_n_meta {
  input {
    Array[File] single_submission_fasta
    Array[File] single_submission_meta
    Array[String] samplename
    Array[String] submission_id
    Array[String] vadr_num_alerts
    Int vadr_threshold = 0
    String repository
    String file_ext
    String date
    String docker = "quay.io/theiagen/utility:1.1"
    Int memory = 8
    Int cpu = 4
    Int disk_size = 100
    Int preemptible = 0
  }
  command <<<
  assembly_array=(~{sep=' ' single_submission_fasta})
  assembly_array_len=$(echo "${#assembly_array[@]}")
  meta_array=(~{sep=' ' single_submission_meta})
  meta_array_len=$(echo "${#meta_array[@]}")
  vadr_string="~{sep=',' vadr_num_alerts}"
  IFS=',' read -r -a vadr_array <<< ${vadr_string}
  vadr_array_len=$(echo "${#vadr_array[@]}")
  samplename_array=(~{sep=' ' samplename})
  samplename_array_len=$(echo "${#samplename_array[@]}")
  submission_id_array=(~{sep=' ' submission_id})
  submission_id_array_len=$(echo "${#submission_id_array[@]}")
  vadr_array_len=$(echo "${#vadr_array[@]}")
  passed_assemblies=""
  passed_meta=""

  #Create files to capture batched and excluded samples
  echo -e "~{repository} Identifier\tSamplename\tNumber of Vadr Alerts\tNote" > ~{repository}_batched_samples_~{date}.tsv
  echo -e "~{repository} Identifier\tSamplename\tNumber of Vadr Alerts\tNote" > ~{repository}_excluded_samples_~{date}.tsv

  #Ensure assembly, meta, and vadr arrays are of equal length
  echo "Samples: $samplename_array_len, Assemblies: $assembly_array_len, Metadata: $meta_array_len, Vadr: $vadr_array_len, Submission_IDs: $submission_id_array_len"
  if [ "$submission_id_array_len" -ne "$vadr_array_len" ]; then
    echo "Submission ID array (length: $submission_id_array_len) and vadr array (length: $vadr_array_len) are of unequal length." >&2
    exit 1
  else 
    echo "Submission ID array (length: $submission_id_array_len) and vadr array (length: $vadr_array_len) are of equal length."
  fi

  #remove samples that excede vadr threshold
  for index in ${!submission_id_array[@]}; do
    submission_id=${submission_id_array[$index]}
    samplename=${samplename_array[$index]}
    vadr=${vadr_array[$index]} 
    batch_note=""
    
    # check if the sample has submittable assembly file; if so remove those that excede vadr thresholds
    assembly=$(printf '%s\n' "${assembly_array[@]}" | grep "${submission_id}")
    metadata=$(printf '%s\n' "${meta_array[@]}" | grep "${submission_id}")

    echo -e "Submission_ID: ${submission_id}\n\tAssembly: ${assembly}\n\tMetadata: ${metadata}\n\tVADR: ${vadr}"
    
    if [ \( ! -z "${assembly}" \) -a \( ! -z "{$metadata}" \) ]; then
      repository_identifier=$(grep -e ">" ${assembly} | sed 's/\s.*$//' | sed 's/>//g' )  
      re='^[0-9]+$'
      if ! [[ "${vadr}" =~ $re ]] ; then
        batch_note="No VADR value to evaluate"
        echo -e "\t$submission_id removed: ${batch_note}"
        echo -e "$repository_identifier\t$samplename\t$vadr\t$batch_note" >> ~{repository}_excluded_samples_~{date}.tsv
      elif [ "${vadr}" -le "~{vadr_threshold}" ] ; then
        passed_assemblies=( "${passed_assemblies[@]}" "${assembly}")
        passed_meta=( "${passed_meta[@]}" "${metadata}")
        echo -e "\t$submission_id added to batch"
        echo -e "$repository_identifier\t$samplename\t$vadr\t$batch_note" >> ~{repository}_batched_samples_~{date}.tsv
      else 
        batch_note="Number of vadr alerts (${vadr}) exceeds threshold ~{vadr_threshold}"
        echo -e "\t$submission_id removed: ${batch_note}"
        echo -e "$repository_identifier\t$samplename\t$vadr\t$batch_note" >> ~{repository}_excluded_samples_~{date}.tsv
      fi
    else 
      batch_note="Assembly or metadata file missing" 
      repository_identifier="NA"
      echo -e "\t$submission_id removed: ${batch_note}"
      echo -e "$repository_identifier\t$samplename\t$vadr\t$batch_note" >> ~{repository}_excluded_samples_~{date}.tsv
    fi

  done

  count=0
  for i in ${passed_meta[*]}; do
      # grab header from first sample in meta_array
      while [ "$count" -lt 1 ]; do
        head -n -1 $i > ~{repository}_upload_meta_~{date}.~{file_ext}
        count+=1
      done
      #populate csv with each samples metadata
      tail -n1 $i >> ~{repository}_upload_meta_~{date}.~{file_ext}
  done

  cat ${passed_assemblies[*]} > ~{repository}_upload_~{date}.fasta
  >>>
  output {
    File? upload_meta   = "${repository}_upload_meta_~{date}.~{file_ext}"
    File? upload_fasta  = "${repository}_upload_~{date}.fasta"
    File batched_samples = "${repository}_batched_samples_~{date}.tsv"
    File excluded_samples = "${repository}_excluded_samples_~{date}.tsv"
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    disk: disk_size + " GB" # TES
    preemptible: preemptible
    maxRetries: 3
  }
}

task compile_biosamp_n_sra {
input {
  Array[File] single_submission_biosample_attirbutes
  Array[File] single_submission_sra_metadata
  Array[String] single_submission_sra_reads
  String date
  String? gcp_bucket
  String docker = "quay.io/theiagen/utility:1.1"
  Int memory = 16
  Int cpu = 4
  Int disk_size = 100
  Int preemptible = 0
}
  command <<<
  biosample_attributes_array=(~{sep=' ' single_submission_biosample_attirbutes})
  biosample_attributes_array_len=$(echo "${#biosample_attributes_array[@]}")
  sra_metadata_array=(~{sep=' ' single_submission_sra_metadata})
  sra_metadata_array_len=$(echo "${#sra_metadata_array[@]}")
  sra_reads_array="~{sep=' ' single_submission_sra_reads}"
  sra_reads_arra_len=$(echo "${#sra_reads_arra[@]}")
  
  # Compile BioSample attributes
  biosamp_count=0
  for i in ${biosample_attributes_array[*]}; do
      # grab header from first sample in meta_array
      while [ "${biosamp_count}" -lt 1 ]; do
        head -n -1 $i > biosample_attributes_~{date}.tsv
        biosamp_count+=1
      done
      #populate csv with each samples metadata
      tail -n1 $i >> biosample_attributes_~{date}.tsv
  done
  
  # Compile SRA metadata
  sra_count=0
  for i in ${sra_metadata_array[*]}; do
      # grab header from first sample in meta_array
      while [ "${sra_count}" -lt 1 ]; do
        head -n -1 $i > sra_metadata_~{date}.tsv
        sra_count+=1
      done
      #populate csv with each samples metadata
      tail -n1 $i >> sra_metadata_~{date}.tsv
  done
  
  # move sra read data to gcp bucket if one is specified; zip into single file if not
  if [[ ! -z "~{gcp_bucket}" ]]
  then 
    echo "Moving read data to provided GCP Bucket ~{gcp_bucket}"
    echo "Running: gsutil -m cp -n ${sra_reads_array[@]} ~{gcp_bucket}"
    gsutil -m cp -n ${sra_reads_array[@]} ~{gcp_bucket}       
  else 
    echo "Preparing SRA read data into single zipped-file"
    mkdir sra_reads_~{date} 
    for i in ${sra_reads_array[*]}; do
      mv $i sra_reads_~{date}
    done  
    zip -r sra_reads_~{date}.zip sra_reads_~{date}
  fi
  >>>
  output {
    File biosample_attributes = "biosample_attributes_~{date}.tsv"
    File sra_metadata = "sra_metadata_~{date}.tsv"
    File? sra_zipped = "sra_reads_~{date}.zip"
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    disk: disk_size + " GB" # TES
    preemptible: preemptible
    maxRetries: 3
  }
}