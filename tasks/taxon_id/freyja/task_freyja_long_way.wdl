version 1.0

task freyja_long_format_single {
    input {
        String samplename
        String freyja_lineages
        String freyja_abundances
        String? collection_date
        String? collection_site
        Float? latitude
        Float? longitude
        String? group_by
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/freyja-microreact:1.0.0"
        Int disk_size = 50
        Int memory = 4
        Int cpu = 2
    }

    command <<<
        # create a metadata file for the sample
        header="samplename\tfreyja_lineages\tfreyja_abundances\tcollection_date\tcollection_site"
        row="~{samplename}\t~{freyja_lineages}\t~{freyja_abundances}\t~{collection_date}\t~{collection_site}"

        # cheb if latitude and longitude are provided, add them to the header and row
        if [ -n "~{latitude}" ]; then
            header="${header}\tlatitude"
            row="${row}\t~{latitude}"
        fi
        
        # long
        if [ -n "~{longitude}" ]; then
            header="${header}\tlongitude"
            row="${row}\t~{longitude}"
        fi

        # write the metadata to a tsv file
        echo -e "${header}\n${row}" > freyja_metadata.tsv

        # freyja the long way, if groupby is provided, set --group-by to the value of group_by
        if [ -n "~{group_by}" ]; then
            freyja_to_long.py freyja_metadata.tsv ~{samplename}_freyja_long_format.tsv --sample-col samplename --group-by ~{group_by}
        else
            freyja_to_long.py freyja_metadata.tsv ~{samplename}_freyja_long_format.tsv --sample-col samplename
    >>>

    output {
        File freyja_long_format_tsv = "~{samplename}_freyja_long_format.tsv"
    }
    runtime {
        docker: docker
        disk: disk_size
        memory: memory
        cpu: cpu
    }
}

task freyja_long_format_multi {
    input {
        Array[String] samplenames
        Array[String] freyja_lineages
        Array[String] freyja_abundances
        Array[String]? collection_dates
        Array[String]? collection_sites
        Array[Float]? latitudes
        Array[Float]? longitudes
        String? group_by
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/freyja-microreact:1.0.0"
        Int disk_size = 50
        Int memory = 4
        Int cpu = 2
    }

    command <<<
        python3 <<CODE
        def read_lines(path):
            with open(path) as f:
                content = f.read().strip()
            return content.split('\n') if content else []

        samplenames = read_lines("~{write_lines(samplenames)}")
        lineages    = read_lines("~{write_lines(freyja_lineages)}")
        abundances  = read_lines("~{write_lines(freyja_abundances)}")
        dates       = read_lines("~{if defined(collection_dates) then write_lines(select_first([collection_dates])) else write_lines([])}")
        sites       = read_lines("~{if defined(collection_sites) then write_lines(select_first([collection_sites])) else write_lines([])}")
        lats        = read_lines("~{if defined(latitudes) then write_lines(select_first([latitudes])) else write_lines([])}")
        lons        = read_lines("~{if defined(longitudes) then write_lines(select_first([longitudes])) else write_lines([])}")

        headers = ['samplename', 'freyja_lineages', 'freyja_abundances', 'collection_date', 'collection_site']
        if lats:
            headers.append('latitude')
        if lons:
            headers.append('longitude')

        with open('freyja_metadata.tsv', 'w') as f:
            f.write('\t'.join(headers) + '\n')
            for i in range(len(samplenames)):
                row = [samplenames[i], lineages[i], abundances[i],
                       dates[i] if dates else '',
                       sites[i] if sites else '']
                if lats:
                    row.append(lats[i])
                if lons:
                    row.append(lons[i])
                f.write('\t'.join(row) + '\n')
        CODE

        if [ -n "~{group_by}" ]; then
            freyja_to_long.py freyja_metadata.tsv freyja_long_format.tsv --sample-col samplename --group-by ~{group_by}
        else
            freyja_to_long.py freyja_metadata.tsv freyja_long_format.tsv --sample-col samplename
    >>>

    output {
        File freyja_long_format_tsv = "freyja_long_format.tsv"
    }

    runtime {
        docker: docker
        disk: disk_size
        memory: memory
        cpu: cpu
    }
}