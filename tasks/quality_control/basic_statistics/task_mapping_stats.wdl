version 1.0

task mapping_stats {
  input {
    File bamfile # aligned reads
    String samplename
    File read1
    File? read2
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/pysam:1.23"
  }
  command <<<
    # fail hard
    set -euo pipefail

    date | tee DATE
    samtools --version | head -n1 | tee VERSION

    samtools stats ~{bamfile} > ~{samplename}.stats.txt
    samtools coverage ~{bamfile} -m -o ~{samplename}.cov.hist
    samtools coverage ~{bamfile} -o ~{samplename}.cov.txt
    samtools flagstat ~{bamfile} > ~{samplename}.flagstat.txt

     # Extracting coverage, depth, meanbaseq, and meanmapq
    python3 <<CODE
    import json
    with open("~{samplename}.cov.txt") as f:
      lines = f.readlines()

    # compile data for each sequence into a dictionary
    seq2data = {}
    for line in lines[1:]: # Skip header
      seq_data = line.strip().split("\t")
      seq, start, end, reads, cov_bases, coverage, meandepth, meanbaseq, meanmapq = seq_data
      seq2data[seq] = {
        "length": int(end) - int(start),
        "reads": int(reads),
        "cov_bases": int(cov_bases),
        "coverage": float(coverage),
        "meandepth": float(meandepth),
        "meanbaseq": float(meanbaseq),
        "meanmapq": float(meanmapq)
      }

    # sort by name
    seq2data = {k: v for k, v in sorted(seq2data.items(), key=lambda x: x[0])}

    # calculate averages across all sequences
    total_length = sum(part["length"] for part in seq2data.values())
    total_mapped_reads = sum(part["reads"] for part in seq2data.values())
    total_cov_bases = sum(part["cov_bases"] for part in seq2data.values())

    total_breadth = total_cov_bases / total_length
    total_depth = sum(part["meandepth"] * part["length"] for part in seq2data.values()) / total_length
    total_meanbaseq = sum(part["meanbaseq"] * part["cov_bases"] for part in seq2data.values()) / total_cov_bases
    total_meanmapq = sum(part["meanmapq"] * part["reads"] for part in seq2data.values()) / total_mapped_reads

    # report total stats
    with open("COVERAGE", "w") as f:
      f.write(str(total_breadth))
    with open("DEPTH", "w") as f:
      f.write(str(total_depth))
    with open("MEANBASEQ", "w") as f:
      f.write(str(total_meanbaseq))
    with open("MEANMAPQ", "w") as f:
      f.write(str(total_meanmapq))

    # report sequence specific mapping stats
    with open("SEQ2COVERAGE.json", "w") as f:
      json.dump({seq: part["coverage"] for seq, part in seq2data.items()}, f, indent=4)
    with open("SEQ2DEPTH.json", "w") as f:
      json.dump({seq: part["meandepth"] for seq, part in seq2data.items()}, f, indent=4)
    CODE

    # parse inputted reads for total read count
    read1_count=$(samtools view -c ~{read1})
    if [ ~{if defined(read2) then "true" else "false"} == "true" ]; then
      read2_count=$(samtools view -c ~{read2})
      total_reads=$(echo $(($read1_count + $read2_count)))
    else
      total_reads=$read1_count
    fi
    # exclude supplementary, unmapped, and secondary alignments from the mapped count
    mapped_reads=$(samtools view -c -F 0x904 ~{bamfile})

    # Check for empty values and set defaults to avoid errors
    if [ -z "$total_reads" ]; then total_reads="1"; fi  # Avoid division by zero
    if [ -z "$mapped_reads" ]; then mapped_reads="0"; fi

    # Calculate the percentage of mapped reads
    percentage_mapped_reads=$(awk "BEGIN {printf \"%.2f\", ($mapped_reads / $total_reads) * 100}")

    # If the percentage calculation fails, default to 0.0
    if [ -z "$percentage_mapped_reads" ]; then percentage_mapped_reads="0.0"; fi

    # Output the result
    echo $percentage_mapped_reads | tee PERCENTAGE_MAPPED_READS
  >>>
  output {
    String date = read_string("DATE")
    String samtools_version = read_string("VERSION")
    Map[String, Float] coverage_by_sequence = read_json("SEQ2COVERAGE.json")
    Map[String, Float] depth_by_sequence = read_json("SEQ2DEPTH.json")
    File stats = "~{samplename}.stats.txt"
    File cov_hist = "~{samplename}.cov.hist"
    File cov_stats = "~{samplename}.cov.txt"
    File flagstat = "~{samplename}.flagstat.txt"
    Float coverage = read_string("COVERAGE")
    Float depth = read_string("DEPTH")
    Float meanbaseq = read_string("MEANBASEQ")
    Float meanmapq = read_string("MEANMAPQ")
    Float percentage_mapped_reads = read_string("PERCENTAGE_MAPPED_READS")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}
