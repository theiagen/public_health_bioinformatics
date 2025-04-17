version 1.0

task srst2_vibrio {
  meta {
    description: "Use of SRST2 to identify sequences of interest from a database of curated Vibrio sequences"
  }
  input {
    File read1
    File? read2
    String samplename
    Int min_percent_coverage
    Int max_divergence
    Int min_depth
    Int min_edge_depth
    Int gene_max_mismatch
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/srst2:0.2.0-vcholerae"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 8
  }
  command <<<
    if [ -z "~{read2}" ] ; then
      INPUT_READS="--input_se ~{read1}"
    else
      # This task expects/requires that input FASTQ files end in "_1.clean.fastq.gz" and "_2.clean.fastq.gz"
      # which is the syntax from TheiaProk read cleaning tasks
      INPUT_READS="--input_pe ~{read1} ~{read2} --forward _1.clean --reverse _2.clean"
    fi

    srst2 --version 2>&1 | tee VERSION
    srst2 \
      ${INPUT_READS} \
      --gene_db /vibrio-cholerae-db/vibrio_230224.fasta \
      --output ~{samplename} \
      --min_coverage ~{min_percent_coverage} \
      --max_divergence ~{max_divergence} \
      --min_depth ~{min_depth} \
      --min_edge_depth ~{min_edge_depth} \
      --gene_max_mismatch ~{gene_max_mismatch}
    
    # capture output TSV
    mv ~{samplename}__genes__*__results.txt ~{samplename}.tsv

    # capture detailed output TSV - not available if no results are outputed
    mv ~{samplename}__fullgenes__*__results.txt ~{samplename}.detailed.tsv || echo "No results" >  ~{samplename}.detailed.tsv

    # parsing block to account for when output columns do not exist
    python <<CODE
    import csv
    import re

    # Converting TSV file into list of dictionaries
    def tsv_to_dict(filename):
      result_list=[]
      with open(filename) as file_obj:
          reader = csv.DictReader(file_obj, delimiter='\t')
          for row in reader:
              result_list.append(dict(row))
      # only one sample is run, so there's only one row, flattening list
      return result_list[0]

    # Converting None to empty string
    conv = lambda i : i or '-'

    # Make characters human-readable 
    def translate_chars(string):
      translation = []
      if '?' in string:
        translation.append("low depth/uncertain")
      if '-' in string:
        translation.append("not detected")
      
      # in case we want to retrieve the allele information
      string = re.sub("\*|\?|-", "", string)

      if len(translation) > 0:
        return '(' + ';'.join(translation) + ')'
      return ""

    # load output TSV as dict 
    row = tsv_to_dict('~{samplename}.tsv')
  
    # presence or absence genes - ctxA, ompW and toxR
    with open("ctxA", "wb") as ctxA_fh:
      value = row.get("ctxA")
      presence = translate_chars(conv(value))
      if presence == "(not detected)":
        ctxA_fh.write(presence)
      else:
        result = "present" + ' ' + presence
        ctxA_fh.write(result.strip())
    
    with open("ompW", "wb") as ompW_fh:
      value = row.get("ompW")
      presence = translate_chars(conv(value))
      if presence == "(not detected)":
        ompW_fh.write(presence)
      else:
        result = "present" + ' ' + presence
        ompW_fh.write(result.strip())
    
    with open("toxR", "wb") as toxR_fh:
      value = row.get("toxR")
      presence = translate_chars(conv(value))
      if presence == "(not detected)":
        toxR_fh.write(presence)
      else:
        result = "present" + ' ' + presence
        toxR_fh.write(result.strip())
    
    # biotype - tcpA classical or tcpA ElTor
    with open("BIOTYPE", "wb") as biotype_fh:
      value_ElTor = translate_chars(conv(row.get("tcpA_ElTor")))
      value_classical = translate_chars(conv(row.get("tcpA_classical")))

      if value_ElTor == "(not detected)" and value_classical == "(not detected)":
        biotype_fh.write("(not detected)")
      else:
        if value_ElTor == "(not detected)":
          result = "tcpA_Classical" + ' ' + value_classical
          biotype_fh.write(result.strip())
        else:
          result = "tcpA_ElTor" + ' ' + value_ElTor
          biotype_fh.write(result.strip())
        
    # serogroup - O1 or O139
    with open("SEROGROUP", "wb") as serotype_fh:
      value_O1 = translate_chars(conv(row.get("wbeN_O1")))
      value_O139 = translate_chars(conv(row.get("wbfR_O139")))

      if value_O1 == "(not detected)" and value_O139 == "(not detected)":
        serotype_fh.write("(not detected)")
      else:
        if value_O1 == "(not detected)":
          result = "O139" + ' ' + value_O139
          serotype_fh.write(result.strip())
        else:
          result = "O1" + ' ' + value_O1
          serotype_fh.write(result.strip())
    CODE
  >>>
  output {
      File srst2_detailed_tsv = "~{samplename}.detailed.tsv"
      String srst2_docker = docker
      String srst2_database = "vibrio_230224"
      String srst2_version = read_string("VERSION")
      String srst2_vibrio_ctxA = read_string("ctxA")
      String srst2_vibrio_ompW = read_string("ompW")
      String srst2_vibrio_toxR = read_string("toxR")
      String srst2_vibrio_biotype = read_string("BIOTYPE")
      String srst2_vibrio_serogroup = read_string("SEROGROUP")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}