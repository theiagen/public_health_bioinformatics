version 1.0

task pmga {
    meta {
        description: "Serogrouping and serotyping of all Neisseria species and Haemophilus influenzae"
    }

    input {
        File assembly
        String samplename
        String docker = "quay.io/staphb/pmga:3.0.2"
        Int? cpu = 4
    }

    command <<<
        echo $(pmga --version 2>&1) | sed 's/.*pmga //; s/ .*\$//' | tee VERSION
        pmga \
            ~{assembly} \
            --blastdir /data/blastdbs \
            --threads ~{cpu} \
            --prefix ~{samplename}

        # Parse pmga TSV
        # https://github.com/rpetit3/pmga#pmga-output-files
        cut -f 2 pmga/~{samplename}.txt | tail -n 1 | tee PMGA_SPECIESDB
        cut -f 3 pmga/~{samplename}.txt | tail -n 1 | tee PMGA_SEROTYPE
        cut -f 4 pmga/~{samplename}.txt | tail -n 1 | tee PMGA_GENES
        cut -f 5 pmga/~{samplename}.txt | tail -n 1 | tee PMGA_NOTES
    >>>

    output {
        String version = read_string("VERSION")
        String pmga_docker = "~{docker}"
        String pmga_speciesdb = read_string("PMGA_SPECIESDB")
        String pmga_serotype = read_string("PMGA_SEROTYPE")
        String pmga_genes = read_string("PMGA_GENES")
        String pmga_notes = read_string("PMGA_NOTES")
        File pmga_results = "./pmga/~{samplename}.txt"
        File pmga_allele_matrix = "./pmga/~{samplename}-allele-matrix.txt"
        File pmga_blast_final = "./pmga/~{samplename}-blast-final-results.json.gz"
        File pmga_blast_raw = "./pmga/~{samplename}-blast-raw-results.json.gz"
        File pmga_loci_counts = "./pmga/~{samplename}-loci-counts.txt"
        File pmga_gff = "./pmga/~{samplename}.gff.gz"
    }

    runtime {
        docker: "~{docker}"
        memory: "8 GB"
        cpu: 4
        disks: "local-disk 50 SSD"
        preemptible: 0
    }
}
