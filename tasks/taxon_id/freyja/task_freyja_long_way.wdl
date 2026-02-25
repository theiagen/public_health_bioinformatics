version 1.0

task freyja_long_way_single {
    input {
        String samplename
        String freyja_lineages
        String freyja_abundances
        String? collection_date
        String? collection_site
        String? latitude
        String? longitude
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/freyja-microreact:0.1.0"
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

        # freyja the long way
        freyja_to_long.py freyja_metadata.tsv ~{samplename}_freyja_long_format.tsv --sample-col samplename
    >>>

    output {
        File freyja_long_format = "~{samplename}_freyja_long_format.tsv"
    }
    runtime {
        docker: docker
        disk: disk_size
        memory: memory
        cpu: cpu
    }
}