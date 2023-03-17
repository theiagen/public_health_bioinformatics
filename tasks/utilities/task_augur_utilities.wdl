version 1.0

task zcat { # in use in old augur workflows
  meta {
    description: "Glue together a bunch of text files that may or may not be compressed (autodetect among gz,xz,bz2,lz4,zst or uncompressed inputs). Optionally compress the output (depending on requested file extension)"
  }
  input {
    Array[File] infiles
    String output_name
    Int cpus = 4
    Int disk_size = 375
  }
  command <<<
  python3 <<CODE
  import os.path
  import gzip, lzma, bz2
  import lz4.frame # pypi library: lz4
  import zstandard as zstd # pypi library: zstandard

  # zstd_open() from: https://github.com/broadinstitute/viral-core/blob/master/util/file.py
  def zstd_open(fname, mode='r', **kwargs):
    '''Handle both text and byte decompression of the file.'''
    if 'r' in mode:
      with open(fname, 'rb') as fh:
        dctx = zstd.ZstdDecompressor()
        stream_reader = dctx.stream_reader(fh)
        if 'b' not in mode:
          text_stream = io.TextIOWrapper(stream_reader, encoding='utf-8-sig')
          yield text_stream
          return
        yield stream_reader
    else:
      with open(fname, 'wb') as fh:
        cctx = zstd.ZstdCompressor(level=kwargs.get('level', 10), threads=util.misc.sanitize_thread_count(kwargs.get('threads', None)))
        stream_writer = cctx.stream_writer(fh)
        if 'b' not in mode:
          text_stream = io.TextIOWrapper(stream_reader, encoding='utf-8')
          yield text_stream
          return
        yield stream_writer

  # magic bytes from here:
  # https://en.wikipedia.org/wiki/List_of_file_signatures
  magic_bytes_to_compressor = {
    b"\x1f\x8b\x08":             gzip.open,      # .gz
    b"\xfd\x37\x7a\x58\x5a\x00": lzma.open,      # .xz
    b"\x42\x5a\x68":             bz2.open,       # .bz2 
    b"\x04\x22\x4d\x18":         lz4.frame.open, # .lz4
    b"\x28\xb5\x2f\xfd":         zstd_open       # .zst (open using function above rather than library function)
  }
  extension_to_compressor = {
    ".gz":   gzip.open,      # .gz
    ".gzip": gzip.open,      # .gz
    ".xz":   lzma.open,      # .xz
    ".bz2":  bz2.open,       # .bz2 
    ".lz4":  lz4.frame.open, # .lz4
    ".zst":  zstd_open,      # .zst (open using function above rather than library function)
    ".zstd": zstd_open       # .zst (open using function above rather than library function)
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

  with open_or_compressed_open("~{output_name}", 'wt') as outf:
    for infname in "~{sep='*' infiles}".split('*'):
      with open_or_compressed_open(infname, 'rt') as inf:
        for line in inf:
          outf.write(line)
  CODE

  # gather runtime metrics
  cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
  cat /proc/loadavg > CPU_LOAD
  { cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes || echo 0; } > MEM_BYTES
  >>>
  runtime {
    docker: "quay.io/broadinstitute/viral-core:2.1.33"
    memory: "1 GB"
    cpu:    cpus
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File combined = "${output_name}"
    Int max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
    Int runtime_sec = ceil(read_float("UPTIME_SEC"))
    String cpu_load = read_string("CPU_LOAD")
  }
}

task fasta_to_ids { ## in use
  meta {
    description: "Return the headers only from a fasta file"
  }
  input {
    File sequences_fasta
    Int disk_size = 375
  }
  String basename = basename(sequences_fasta, ".fasta")
  command <<<
    cat "~{sequences_fasta}" | grep \> | cut -c 2- > "~{basename}.txt"
  >>>
  runtime {
    docker: "ubuntu"
    memory: "1 GB"
    cpu:    1
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File ids_txt = "~{basename}.txt"
  }
}

task tsv_join { ## in use
  meta {
    description: "Perform a full left outer join on multiple TSV tables. Each input tsv must have a header row, and each must must contain the value of id_col in its header. Inputs may or may not be gzipped. Unix/Mac/Win line endings are tolerated on input, Unix line endings are emitted as output. Unicode text safe."
  }
  input {
    Array[File]+   input_tsvs
    String         id_col
    String         out_basename = "merged"
    String         out_suffix = ".txt"
    Int            machine_mem_gb = 7
    Int disk_size = 100
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
  CODE
  >>>
  output {
    File out_tsv = "~{out_basename}~{out_suffix}"
  }
  runtime {
    memory: "~{machine_mem_gb} GB"
    cpu: 2
    docker: "quay.io/broadinstitute/viral-core:2.1.33"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task filter_sequences_by_length { # in use
  meta {
    description: "Filter sequences in a fasta file to enforce a minimum count of non-N bases."
  }
  input {
    File   sequences_fasta
    Int    min_non_N = 1

    String docker = "quay.io/broadinstitute/viral-core:2.1.33"
    Int disk_size = 300
  }
  parameter_meta {
    sequences_fasta: {
      description: "Set of sequences in fasta format",
      patterns: ["*.fasta", "*.fa"]
    }
    min_non_N: {
      description: "Minimum number of called bases (non-N, non-gap, A, T, C, G, and other non-N ambiguity codes accepted)"
    }
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
    memory: "1 GB"
    cpu :   1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File filtered_fasta    = out_fname
    Int  sequences_in      = read_int("IN_COUNT")
    Int  sequences_dropped = read_int("DROP_COUNT")
    Int  sequences_out     = read_int("OUT_COUNT")
  }
}

task set_sc2_defaults { # establish sars-cov-2 default values for augur
  input {
    String nextstrain_ncov_repo_commit = "23d1243127e8838a61b7e5c1a72bc419bf8c5a0d" # last updated on 2023-03-07
    Int disk_size = 50
  }
  command <<<
    wget -q "https://github.com/nextstrain/ncov/archive/~{nextstrain_ncov_repo_commit}.tar.gz"
    tar -xf "~{nextstrain_ncov_repo_commit}.tar.gz" --strip-components=1
  >>>
  output {
    Int min_num_unambig = 27000
    File clades_tsv = "defaults/clades.tsv"
    File lat_longs_tsv = "defaults/lat_longs.tsv"
    File reference_fasta = "defaults/reference_seq.fasta"
    File reference_genbank = "defaults/reference_seq.gb"
    File auspice_config = "defaults/auspice_config.json"
    Float min_date = 2020.0
    Int pivot_interval = 1
    String pivot_interval_units = "weeks"
    Float narrow_bandwidth = 0.05
    Float proportion_wide = 0.0
  }
  runtime {
    docker: "quay.io/staphb/augur:16.0.3"
    memory: "1 GB"
    cpu:   1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
  }
}

task set_flu_defaults { # establish flu default values for augur
  input {
    String flu_segment
    String? flu_subtype

    File flu_lat_longs_tsv = "gs://theiagen-public-files-rp/terra/flu-references/lat_longs.tsv"

    Int disk_size = 50
  }
  command <<<
    # set h1n1 defaults
    if [ ~{flu_subtype} == "H1N1" ]; then
      if [ ~{flu_segment} == "HA" ]; then
        echo "gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_ha.gb" | tee FLU_REFERENCE_FASTA
        echo "gs://theiagen-public-files-rp/terra/flu-references/clades_h1n1pdm_ha.tsv" | tee FLU_CLADE_FILE
      elif [ ~{flu_segment} == "NA" ]; then
        echo "gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_na.gb" | tee FLU_REFERENCE_FASTA
      else 
        echo "Uh oh! Your flu segment was not recognized. The only accepted options are \"HA\" or \"NA\"! You provided: \"~{flu_segment}\"."
        exit 1
      fi
      echo "gs://theiagen-public-files-rp/terra/flu-references/auspice_config_h1n1pdm.json" | tee AUSPICE_CONFIG

    # set h3n2 defaults
    elif [ ~{flu_subtype} == "H3N2" ]; then
      if [ ~{flu_segment} == "HA" ]; then
        echo "gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_ha.gb" | tee FLU_REFERENCE_FASTA
        echo "gs://theiagen-public-files-rp/terra/flu-references/clades_h3n2_ha.tsv" | tee FLU_CLADE_FILE
      elif [ ~{flu_segment} == "NA" ]; then
        echo "gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_na.gb" | tee FLU_REFERENCE_FASTA
      else 
        echo "Uh oh! Your flu segment was not recognized. The only accepted options are \"HA\" or \"NA\"! You provided: \"~{flu_segment}\"."
        exit 1
      fi
      echo "gs://theiagen-public-files-rp/terra/flu-references/auspice_config_h3n2.json" | tee AUSPICE_CONFIG

    # set vic defaults
    elif [ ~{flu_subtype} == "Victoria" ]; then
      if [ ~{flu_segment} == "HA" ]; then
        echo "gs://theiagen-public-files-rp/terra/flu-references/reference_vic_ha.gb" | tee FLU_REFERENCE_FASTA
        echo "gs://theiagen-public-files-rp/terra/flu-references/clades_vic_ha.tsv" | tee FLU_CLADE_FILE
      elif [ ~{flu_segment} == "NA" ]; then
        echo "gs://theiagen-public-files-rp/terra/flu-references/reference_vic_na.gb" | tee FLU_REFERENCE_FASTA
      else 
        echo "Uh oh! Your flu segment was not recognized. The only accepted options are \"HA\" or \"NA\"! You provided: \"~{flu_segment}\"."
        exit 1
      fi
      echo "gs://theiagen-public-files-rp/terra/flu-references/auspice_config_vic.json" | tee AUSPICE_CONFIG

    # set yam defaults
    elif [ ~{flu_subtype} == "Yamagata" ]; then
      if [ ~{flu_segment} == "HA" ]; then
        echo "gs://theiagen-public-files-rp/terra/flu-references/reference_yam_ha.gb" | tee FLU_REFERENCE_FASTA
        echo "gs://theiagen-public-files-rp/terra/flu-references/clades_yam_ha.tsv" | tee FLU_CLADE_FILE
      elif [ ~{flu_segment} == "NA" ]; then
        echo "gs://theiagen-public-files-rp/terra/flu-references/reference_yam_na.gb" | tee FLU_REFERENCE_FASTA
      else 
        echo "Uh oh! Your flu segment was not recognized. The only accepted options are \"HA\" or \"NA\"! You provided: \"~{flu_segment}\"."
        exit 1
      fi
      echo "gs://theiagen-public-files-rp/terra/flu-references/auspice_config_yam.json" | tee AUSPICE_CONFIG

    else 
      echo "Uh oh! Your flu subtype was not recognized. The only accepted options are \"H1N1\", \"H3N2\", \"Victoria\", or \"Yamagata\"! You provided: \"~{flu_subtype}\"."
      exit 1
    fi

  >>>
  output {
    Int min_num_unambig = 900
    File? clades_tsv = read_string("FLU_CLADE_FILE")
    File lat_longs_tsv = flu_lat_longs_tsv
    File reference_fasta = read_string("FLU_REFERENCE_FASTA")
    File reference_genbank = read_string("FLU_REFERENCE_FASTA")
    File auspice_config = read_string("AUSPICE_CONFIG")
     # the following values were set to match the nextstrain build: https://github.com/nextstrain/seasonal-flu/blob/fa7a3d323a7e193f1ba1ae5b3e40249edb87c491/Snakefile_base#L713
    Float min_date = 2020.0
    Int pivot_interval = 1
    Float narrow_bandwidth = 0.1666667
    Float proportion_wide = 0.0
  }
  runtime {
    docker: "quay.io/staphb/augur:16.0.3"
    memory: "1 GB"
    cpu:   1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
  }
}

task prep_augur_metadata {
  input {
    File assembly
    String collection_date
    String country
    String state
    String continent

    String county = ""
    String? pango_lineage

    Int disk_size = 10
  }
  command <<<   
    # Set strain name by assembly header
    assembly_header=$(grep -e ">" ~{assembly} | sed 's/\s.*$//' |  sed 's/>//g' )

    pangolin_header=""

    # if pango_lineage defined, add to metadata
    if [[ "~{pango_lineage}" ]] ; then 
      pangolin_header="pango_lineage"
    fi

    echo -e "strain\tvirus\tdate\tregion\tcountry\tdivision\tlocation\t${pangolin_header}" > augur_metadata.tsv
    echo -e "\"${assembly_header}\"\t\"ncov\"\t\"~{collection_date}\"\t\"~{continent}\" \t\"~{country}\"\t\"~{state}\"\t\"~{county}\"\t\"~{pango_lineage}\"" >> augur_metadata.tsv

  >>>
  output {
    File augur_metadata = "augur_metadata.tsv"
  }
  runtime {
      docker: "quay.io/theiagen/utility:1.1"
      memory: "3 GB"
      cpu: 1
      disks: "local-disk ~{disk_size} SSD"
      disk: disk_size + " GB"
      preemptible: 0
      maxRetries: 3
  }
}