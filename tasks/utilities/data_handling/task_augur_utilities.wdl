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
    description: "Perform a full left outer join on multiple TSV tables. Each input tsv must have a header row, and each must must contain the value of id_col in its header. Inputs may or may not be gzipped. Unix/Mac/Win line endings are tolerated on input, Unix line endings are emitted as output. Unicode text safe. Copied from the Broad Institute"
  }
  input {
    Array[File] input_tsvs
    String id_col
    String out_basename = "merged"
    String out_suffix = ".tsv"
    Int memory = 7
    Int disk_size = 100
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/broadinstitute/viral-core:2.1.33"
  }
  command <<< 
  python3<<CODE
  import collections
  import csv
  import os.path
  import gzip, lzma, bz2
  import lz4.frame # pypi library: lz4
  import zstandard as zstd # pypi library: zstandard
  import util.file # viral-core

  # magic bytes from here:
  # https://en.wikipedia.org/wiki/List_of_file_signatures
  magic_bytes_to_compressor = {
    b"\x1f\x8b\x08":             gzip.open,      # .gz
    b"\xfd\x37\x7a\x58\x5a\x00": lzma.open,      # .xz
    b"\x42\x5a\x68":             bz2.open,       # .bz2
    b"\x04\x22\x4d\x18":         lz4.frame.open, # .lz4
    b"\x28\xb5\x2f\xfd":         util.file.zstd_open   # .zst (open using function above rather than library function)
  }
  extension_to_compressor = {
    ".gz":   gzip.open,      # .gz
    ".gzip": gzip.open,      # .gz
    ".xz":   lzma.open,      # .xz
    ".bz2":  bz2.open,       # .bz2
    ".lz4":  lz4.frame.open, # .lz4
    ".zst":  util.file.zstd_open,  # .zst (open using function above rather than library function)
    ".zstd": util.file.zstd_open   # .zst (open using function above rather than library function)
  }

  # max number of bytes we need to identify one of the files listed above
  max_len = max(len(x) for x in magic_bytes_to_compressor.keys())

  def open_or_compressed_open(*args, **kwargs):
    input_file = args[0]

    # if the file exists, try to guess the (de) compressor based on "magic numbers"
    # at the very start of the file
    if os.path.isfile(input_file):
      with open(input_file, "rb") as f:
        file_start = f.read(max_len)
      for magic, compressor_open_fn in magic_bytes_to_compressor.items():
        if file_start.startswith(magic):
          print("opening via {}: {}".format(compressor_open_fn.__module__,input_file))
          return compressor_open_fn(*args, **kwargs)
      # fall back to generic open if compression type could not be determine from magic numbers
      return open(*args, **kwargs)
    else:
      # if this is a new file, try to choose the opener based on file extension
      for ext,compressor_open_fn in extension_to_compressor.items():
        if str(input_file).lower().endswith(ext):
          print("opening via {}: {}".format(compressor_open_fn.__module__,input_file))
          return compressor_open_fn(*args, **kwargs)
      # fall back to generic open if compression type could not be determine from magic numbers
    return open(*args, **kwargs)

  # prep input readers
  out_basename = '~{out_basename}'
  join_id = '~{id_col}'
  in_tsvs = '~{sep="*" input_tsvs}'.split('*')
  readers = list(
    csv.DictReader(open_or_compressed_open(fn, 'rt'), delimiter='\t')
    for fn in in_tsvs)

  # prep the output header
  header = []
  for reader in readers:
    header.extend(reader.fieldnames)
  header = list(collections.OrderedDict(((h,0) for h in header)).keys())
  if not join_id or join_id not in header:
    raise Exception()

  # merge everything in-memory
  out_ids = []
  out_row_by_id = {}
  for reader in readers:
    for row in reader:
      row_id = row[join_id]
      row_out = out_row_by_id.get(row_id, {})
      for h in header:
        # prefer non-empty values from earlier files in in_tsvs, populate from subsequent files only if missing
        if not row_out.get(h):
          row_out[h] = row.get(h, '')
      out_row_by_id[row_id] = row_out
      out_ids.append(row_id)
  out_ids = list(collections.OrderedDict(((i,0) for i in out_ids)).keys())

  # write output
  with open_or_compressed_open(out_basename+'~{out_suffix}', 'w', newline='') as outf:
    writer = csv.DictWriter(outf, header, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    writer.writerows(out_row_by_id[row_id] for row_id in out_ids)

  if "date" in header:
    has_time = "true"
  else:
    has_time = "false"
  with open('HAS_TIME', 'w') as out:
      out.write(has_time)
  CODE
  >>>
  output {
    File out_tsv = "~{out_basename}~{out_suffix}"
    Boolean has_time = read_boolean("HAS_TIME")
  }
  runtime {
    memory: "~{memory} GB"
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
    echo -e "strain\tvirus\tdate\tregion\tcountry\tdivision\tlocation\t${pangolin_header}\t${nextclade_header}" > augur_metadata.tsv
    echo -e "\"${assembly_header}\"\t\"${virus}\"\t\"~{collection_date}\"\t\"~{region}\"\t\"~{country}\"\t\"~{division}\"\t\"~{location}\"\t\"~{pango_lineage}\"\t\"~{nextclade_clade}\"" >> augur_metadata.tsv
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