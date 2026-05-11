version 1.0

task freyja_long_format_single {
  meta {
    description: "Builds per-sample Freyja metadata and converts it to long format for downstream Microreact visualization."
  }
  input {
    String samplename
    String freyja_lineages
    String freyja_abundances
    Float freyja_coverage
    String? collection_date
    String? collection_site
    Float? latitude
    Float? longitude
    String? group_by
    Int mincov
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/freyja-microreact:1.0.1"
    Int disk_size = 50
    Int memory = 4
    Int cpu = 2
  }
  command <<<
    # build per-sample metadata file
    header="samplename\tfreyja_lineages\tfreyja_abundances\tfreyja_coverage\tcollection_date\tcollection_site"
    row="~{samplename}\t~{freyja_lineages}\t~{freyja_abundances}\t~{freyja_coverage}\t~{collection_date}\t~{collection_site}"

    # check if latitude and longitude are provided, add them to the header and row
    if [ -n "~{latitude}" ]; then
      header="${header}\tlatitude"
      row="${row}\t~{latitude}"
    fi

    if [ -n "~{longitude}" ]; then
      header="${header}\tlongitude"
      row="${row}\t~{longitude}"
    fi

    echo -e "${header}\n${row}" > freyja_metadata.tsv

    # if group_by is provided, pass it through to freyja_to_long.py
    if [ -n "~{group_by}" ]; then
      freyja_to_long.py freyja_metadata.tsv ~{samplename}_freyja_demixed_parsed_long.tsv \
        --sample-col samplename \
        --group-by ~{group_by} \
        --mincov ~{mincov}
    else
      freyja_to_long.py freyja_metadata.tsv ~{samplename}_freyja_demixed_parsed_long.tsv \
        --sample-col samplename \
        --mincov ~{mincov}
    fi
  >>>
  output {
    File freyja_parsed_format_tsv = "~{samplename}_freyja_demixed_parsed_long.tsv"
    String freyja_long_format_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
  }
}

task freyja_long_format_multi {
  meta {
    description: "Builds a multi-sample Freyja metadata table and converts it to long format for downstream Microreact visualization."
  }
  input {
    Array[String] samplenames
    Array[String] freyja_lineages
    Array[String] freyja_abundances
    Array[Float] freyja_coverages
    Array[String]? collection_dates
    Array[String]? collection_sites
    Array[Float]? latitudes
    Array[Float]? longitudes
    String? group_by
    Int mincov
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/freyja-microreact:1.0.1"
    Int disk_size = 50
    Int memory = 4
    Int cpu = 2
  }
  command <<<
    # assemble multi-sample metadata table from the per-array inputs; optional
    # collection arrays are materialized as empty files via write_lines([]) when undefined
    python3 <<CODE
    def read_lines(path):
      with open(path) as f:
        content = f.read().strip()
      return content.split('\n') if content else []

    samplenames = read_lines("~{write_lines(samplenames)}")
    lineages = read_lines("~{write_lines(freyja_lineages)}")
    abundances = read_lines("~{write_lines(freyja_abundances)}")
    dates = read_lines("~{if defined(collection_dates) then write_lines(select_first([collection_dates])) else write_lines([])}")
    sites = read_lines("~{if defined(collection_sites) then write_lines(select_first([collection_sites])) else write_lines([])}")
    lats = read_lines("~{if defined(latitudes) then write_lines(select_first([latitudes])) else write_lines([])}")
    lons = read_lines("~{if defined(longitudes) then write_lines(select_first([longitudes])) else write_lines([])}")
    coverages = read_lines("~{write_lines(freyja_coverages)}")

    headers = ['samplename', 'freyja_lineages', 'freyja_abundances', 'freyja_coverage', 'collection_date', 'collection_site']
    if lats:
      headers.append('latitude')
    if lons:
      headers.append('longitude')

    with open('freyja_metadata.tsv', 'w') as f:
      f.write('\t'.join(headers) + '\n')
      for i in range(len(samplenames)):
        row = [samplenames[i], lineages[i], abundances[i], coverages[i],
               dates[i] if dates else '',
               sites[i] if sites else '']
        if lats:
          row.append(lats[i])
        if lons:
          row.append(lons[i])
        f.write('\t'.join(row) + '\n')
    CODE

    if [ -n "~{group_by}" ]; then
      freyja_to_long.py freyja_metadata.tsv freyja_demixed_parsed_long.tsv \
        --sample-col samplename \
        --group-by ~{group_by} \
        --mincov ~{mincov}
    else
      freyja_to_long.py freyja_metadata.tsv freyja_demixed_parsed_long.tsv \
        --sample-col samplename \
        --mincov ~{mincov}
    fi
  >>>
  output {
    File freyja_parsed_format_tsv = "freyja_demixed_parsed_long.tsv"
    String freyja_long_format_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
  }
}
