version 1.0

task ete3_taxon_id {
  input {
    String taxon # can be taxon id (int) or organism name (string)
    String? rank # limit input taxon to the user specified rank
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/ete3:3.1.3"
    Int cpu = 1
    Int memory = 4
    Int disk_size = 50
  }
  command <<<
    # fail hard
    set -euo pipefail

    date | tee DATE
    ete3 version | sed 's| Tool path.*||' | tee ETE3_VERSION

    echo "DEBUG: Obtaining taxon report for taxon: ~{taxon} and rank: ~{rank}"

    python3 <<CODE
    from ete3 import NCBITaxa

    taxa = NCBITaxa()
    # check if the taxon is an integer taxid or a string name
    try: 
      query_taxid = int("~{taxon}")
    except ValueError:
      try:
        query_taxid = taxa.get_name_translator(["~{taxon}"]).get("~{taxon}")[0]
      # error out if the taxon was not recovered
      except TypeError:
        raise ValueError(f"ERROR: Invalid taxon name: '~{taxon}'")

    lineage = taxa.get_lineage(query_taxid)
    lineage2rank = taxa.get_rank(lineage)
    # will overwrite multiple "no rank" entries, but these won't be used later
    rank2lineage = {v: k for k, v in lineage2rank.items()}
    query_rank = list(taxa.get_rank([query_taxid]).keys())[0]
    query_taxon = taxa.get_taxid_translator([query_taxid])[query_taxid]
    print(f"DEBUG: Resolved input taxon '~{taxon}' to taxid: {query_taxid}, name: {query_taxon}, rank: {query_rank}")

    check_rank = "~{rank}".lower()

    # Reported taxon info is based on the 'rank' input (if provided). Otherwise, use raw rank.
    if not check_rank:
      print("DEBUG: No input 'rank' provided, defaulting to raw taxon ID, name, and rank.")
      reported_taxon_name = raw_taxon_name
      reported_taxon_rank = raw_taxon_rank
    elif check_rank not in rank2lineage:
      reported_taxon_name = reported_taxon_id = reported_taxon_rank = 'N/A'
      raise ValueError(f"ERROR: Input taxon rank '{check_rank}' is not valid for taxon: '{query_taxon}'.")
    else:
      reported_taxon_id = rank2lineage[check_rank]
      reported_taxon_name = taxa.get_taxid_translator([reported_taxon_id])[reported_taxon_id]
      reported_taxon_rank = check_rank

    outputs = {
      "TAXON_ID": str(reported_taxon_id),
      "TAXON_NAME": str(reported_taxon_name).lower(),
      "TAXON_RANK": str(reported_taxon_rank).lower(),
      "RAW_TAXON_ID": str(query_taxid),
      "RAW_TAXON_NAME": str(query_taxon).lower(),
      "RAW_TAXON_RANK": str(query_rank).lower()
    }
    for filename, value in outputs.items():
      with open(filename, "w") as f:
        f.write(value)
    CODE
    echo "DEBUG: Reported TAXON_ID: $(cat TAXON_ID), TAXON_NAME: $(cat TAXON_NAME), TAXON_RANK: $(cat TAXON_RANK)"
  >>>
  output {
    String taxon_id = read_string("TAXON_ID")
    String taxon_name = read_string("TAXON_NAME")
    String taxon_rank = read_string("TAXON_RANK")
    String raw_taxon_id = read_string("RAW_TAXON_ID")
    String raw_taxon_name = read_string("RAW_TAXON_NAME")
    String raw_taxon_rank = read_string("RAW_TAXON_RANK")
    String ete3_version = read_string("ETE3_VERSION")
    String ete3_docker = docker
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}