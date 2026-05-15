version 1.0

task fasta_to_ids {
  meta {
    description: "Return the headers only from a fasta file; copied from the Broad Institute"
  }
  input {
    File sequences_fasta
    Int disk_size = 375
    Int memory = 1
    Int cpu = 1
    String docker = "ubuntu"
  }
  String basename = basename(sequences_fasta, ".fasta")
  command <<<
    cat "~{sequences_fasta}" | grep \> | cut -c 2- > "~{basename}.txt"
  >>>
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File ids_txt = "~{basename}.txt"
  }
}

task tsv_join {
  meta {
    description: "Concatenates TSV files and checks if they have a collection date"
  }
  input {
    Array[File] input_tsvs
    Int memory = 4
    Int disk_size = 100
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0"
  }
  command <<<
    input_tsvs=(~{sep=' ' input_tsvs})

    for index in ${!input_tsvs[@]}; do
      file=${input_tsvs[$index]}
      if [ $index == 0 ]; then
        cat ${file} >> metadata-merged.tsv
      else
        tail -n +2 ${file} >> metadata-merged.tsv
      fi
    done

    if (head -n1 merged.tsv | grep -q "date"); then
      echo true | tee HAS_TIME
    else
      echo false | tee HAS_TIME
    fi
  >>>
  output {
    File out_tsv = "metadata-merged.tsv"
    Boolean has_time = read_boolean("HAS_TIME")
  }
  runtime {
    memory: memory + " GB"
    cpu: cpu
    docker: docker
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task filter_sequences_by_length {
  meta {
    description: "Filter sequences in a fasta file to enforce a minimum count of non-N bases; copied from the Broad Institute"
  }
  input {
    File sequences_fasta
    Int min_non_N = 1 # Minimum number of called bases (non-N, non-gap, A, T, C, G, and other non-N ambiguity codes)

    String docker = "us-docker.pkg.dev/general-theiagen/broadinstitute/viral-core:2.1.33"
    Int disk_size = 300
    Int memory = 1
    Int cpu = 1
  }
  String out_fname = sub(basename(sequences_fasta), ".fasta", ".filtered.fasta")
  command <<<
  python3 <<CODE
  import Bio.SeqIO
  import gzip

  n_total = 0
  n_kept = 0

  open_or_gzopen = lambda *args, **kwargs: gzip.open(*args, **kwargs) if args[0].endswith('.gz') else open(*args, **kwargs)

  with open_or_gzopen('~{sequences_fasta}', 'rt') as inf:
    with open_or_gzopen('~{out_fname}', 'wt') as outf:
      for seq in Bio.SeqIO.parse(inf, 'fasta'):
        n_total += 1
        ungapseq = seq.seq.ungap().upper()
        if (len(ungapseq) - ungapseq.count('N')) >= ~{min_non_N}:
          n_kept += 1
          Bio.SeqIO.write(seq, outf, 'fasta')
  
  n_dropped = n_total-n_kept
  
  with open('IN_COUNT', 'wt') as outf:
    outf.write(str(n_total)+'\n')
  with open('OUT_COUNT', 'wt') as outf:
    outf.write(str(n_kept)+'\n')
  with open('DROP_COUNT', 'wt') as outf:
    outf.write(str(n_dropped)+'\n')
  CODE
  >>>
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File filtered_fasta = out_fname
    Int sequences_in = read_int("IN_COUNT")
    Int sequences_dropped = read_int("DROP_COUNT")
    Int sequences_out = read_int("OUT_COUNT")
  }
}

task prep_augur_metadata {
  input {
    File assembly
    String? collection_date
    String? country
    String? division
    String? region

    String? location
    String? pango_lineage
    String? nextclade_clade
    String? organism = "sars-cov-2"

    Int disk_size = 10
    Int memory = 3
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
  }
  command <<<   
    # Set strain name by assembly header
    assembly_header=$(grep -e ">" ~{assembly} | sed 's/\s.*$//' |  sed 's/>//g' )

    pangolin_header=""
    nextclade_header=""

    # if pango_lineage defined, add to metadata
    if [[ "~{pango_lineage}" ]]; then 
      pangolin_header="pango_lineage"
    fi
    # if pango_lineage defined, add to metadata
    if [[ "~{nextclade_clade}" ]]; then 
      nextclade_header="clade_membership"
    fi

    if [[ "~{organism}" == "sars-cov-2" ]]; then
      virus="ncov"
    else
      virus="~{organism}"
    fi

    # write everything to a file
    echo -e "strain\tvirus\t~{if defined(collection_date) then 'date' else ''}\tregion\tcountry\tdivision\tlocation\t${pangolin_header}\t${nextclade_header}" > augur_metadata.tsv
    echo -e "${assembly_header}\t${virus}\t~{collection_date}\t~{region}\t~{country}\t~{division}\t~{location}\t~{pango_lineage}\t~{nextclade_clade}" >> augur_metadata.tsv
  >>>
  output {
    File augur_metadata = "augur_metadata.tsv"
  }
  runtime {
      docker: docker
      memory: memory + " GB"
      cpu: cpu
      disks: "local-disk ~{disk_size} SSD"
      disk: disk_size + " GB"
      preemptible: 0
      maxRetries: 3
  }
}
