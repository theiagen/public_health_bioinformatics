version 1.0

task kraken_output_parser {
  input {
    File kraken2_report
    Array[String] taxon_ids
    Int? read_count_threshold
    
    Int cpu = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/krakentools:d4a2fbe"
    Int memory = 4
  }
  meta {
    description: "Parse Kraken2 report to extract taxon IDs identified within the sample. Intended for use within TheiaViral_Panel."
  }
  command <<<
    set -euo pipefail

    python3 <<CODE
    taxon_array = "~{sep=' ' taxon_ids}".split(' ')
    parsed_array = []
    for line in open("~{kraken2_report}"):
        line = line.strip().split("\t")
        taxon_id = str(line[4])
        read_count= int(line[1])
        if (taxon_id in taxon_array):
          if read_count >= ~{read_count_threshold}:
            parsed_array.append(taxon_id)
          else:
            print(f"Taxon ID {taxon_id} has {read_count} reads, below threshold of ~{read_count_threshold}")
    with open("parsed_ids.txt", "w") as out:
        out.write("\\n".join(parsed_array) + "\\n")
        out.close()
    CODE

  >>>
  output {
    Array[String] parsed_taxon_ids = read_lines("parsed_ids.txt")
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
