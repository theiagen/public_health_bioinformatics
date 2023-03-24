version 1.0

task theiacov_summary {
  input {
    String samplename
    String theiacov_workflow
    String theiacov_version
    String theiacov_analysis_date
    String seq_platform
    String primer_bed_name
    Float percent_reference_coverage
    Int? consensus_n_variant_min_depth
    Float? s_gene_mean_coverage
    Float? s_gene_percent_coverage
    Int number_N
    String pango_lineage
    String pangolin_conflicts
    String pangolin_notes
    String pangolin_assignment_version
    String pangolin_docker
    String pango_lineage_expanded
    String pangolin_versions
    String nextclade_clade
    String nextclade_aa_subs
    String nextclade_aa_dels
    String vadr_num_alerts
    Int assembly_length_unambiguous
    Float assembly_mean_coverage
    String assembly_method
    Int number_Degenerate
    Int number_Total
    Float meanbaseq_trim
    Float meanmapq_trim
    Int num_reads_clean1
    String num_reads_clean2 = ""
    String num_reads_clean_pairs = ""
    Int num_reads_raw1
    String num_reads_raw2 = ""
    String num_reads_raw_pairs = ""
    Float kraken_human
    String kraken_human_dehosted = ""
    Float kraken_sc2
    String kraken_sc2_dehosted = ""
    String? primer_trimmed_read_percent
    String? artic_version
    String? artic_docker
    String? medaka_reference
    String bbduk_docker = ""
    String bwa_version = ""
    String fastq_scan_version = ""
    String ivar_variant_version = ""
    String ivar_version_consensus = ""
    String ivar_version_primtrim = ""
    String kraken_version = ""
    String nextclade_version = ""
    String nextclade_docker = ""
    String nextclade_ds_tag = ""
    String samtools_version
    String samtools_version_consensus = ""
    String samtools_version_primtrim = ""
    String samtools_version_stats = ""
    String trimmomatic_version = ""
    String vadr_docker = ""
  }
  command <<<
      python3<<CODE
      from collections import OrderedDict
      import json
      data = OrderedDict((
          ('sample', '~{samplename}'),
          ('theiacov_workflow', '~{theiacov_workflow}'),
          ('theiacov_version', '~{theiacov_version}'),
          ('theiacov_analysis_date', '~{theiacov_analysis_date}'),
          ('seq_platform', '~{seq_platform}'),
          ('primer_bed_name', '~{primer_bed_name}'),
          ('percent_reference_coverage', '~{percent_reference_coverage}'),
          ('consensus_n_variant_min_depth', '~{consensus_n_variant_min_depth}'),
          ('s_gene_mean_coverage', '~{s_gene_mean_coverage}'),
          ('s_gene_percent_coverage', '~{s_gene_percent_coverage}'),
          ('number_n', '~{number_N}'),
          ('pangolin_lineage', '~{pango_lineage}'),
          ('pangolin_conflicts', '~{pangolin_conflicts}'),
          ('pangolin_notes', '~{pangolin_notes}'),
          ('pangolin_assignment_version', '~{pangolin_assignment_version}'),
          ('pango_lineage_expanded', '~{pango_lineage_expanded}'),
          ('pangolin_versions', '~{pangolin_versions}'),
          ('pangolin_docker', '~{pangolin_docker}'),
          ('nextclade_clade', '~{nextclade_clade}'),
          ('nextclade_aa_subs', '~{nextclade_aa_subs}'),
          ('nextclade_aa_dels', '~{nextclade_aa_dels}'),
          ('vadr_num_alerts', '~{vadr_num_alerts}'),
          ('assembly_length_unambiguous', '~{assembly_length_unambiguous}'),
          ('assembly_mean_coverage', '~{assembly_mean_coverage}'),
          ('assembly_method', '~{assembly_method}'),
          ('number_degenerate', '~{number_Degenerate}'),
          ('number_total', '~{number_Total}'),
          ('meanbaseq_trim', '~{meanbaseq_trim}'),
          ('meanmapq_trim', '~{meanmapq_trim}'),
          ('num_reads_clean1', '~{num_reads_clean1}'),
          ('num_reads_clean2', '~{num_reads_clean2}'),
          ('num_reads_clean_pairs', '~{num_reads_clean_pairs}'),
          ('num_reads_raw1', '~{num_reads_raw1}'),
          ('num_reads_raw2', '~{num_reads_raw2}'),
          ('num_reads_raw_pairs', '~{num_reads_raw_pairs}'),
          ('kraken_human', '~{kraken_human}'),
          ('kraken_human_dehosted', '~{kraken_human_dehosted}'),
          ('kraken_sc2', '~{kraken_sc2}'),
          ('kraken_sc2_dehosted', '~{kraken_sc2_dehosted}'),
          ('primer_trimmed_read_percent', '~{primer_trimmed_read_percent}'),
          ('artic_version', '~{artic_version}'),
          ('artic_docker', '~{artic_docker}'),
          ('medaka_reference', '~{medaka_reference}'),
          ('bbduk_docker', '~{bbduk_docker}'),
          ('bwa_version', '~{bwa_version}'),
          ('fastq_scan_version', '~{fastq_scan_version}'),
          ('ivar_variant_version', '~{ivar_variant_version}'),
          ('ivar_version_consensus', '~{ivar_version_consensus}'),
          ('ivar_version_primtrim', '~{ivar_version_primtrim}'),
          ('kraken_version', '~{kraken_version}'),
          ('nextclade_version', '~{nextclade_version}'),
          ('nextclade_docker', '~{nextclade_docker}'),
          ('nextclade_ds_tag', '~{nextclade_ds_tag}'),
          ('samtools_version', '~{samtools_version}'),
          ('samtools_version_consensus', '~{samtools_version_consensus}'),
          ('samtools_version_primtrim', '~{samtools_version_primtrim}'),
          ('samtools_version_stats', '~{samtools_version_stats}'),
          ('trimmomatic_version', '~{trimmomatic_version}'),
          ('vadr_docker', '~{vadr_docker}')
      ))
      with open(f'~{samplename}.results.json', "wt") as json_fh:
          json.dump(data, json_fh)
      CODE
  >>>
  output {
    File summary = '~{samplename}.results.json'
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}

task merge_theiacov_summary {
  input {
    Array[File?] clearlabs_summaries
    Array[File?] illumina_pe_summaries
    Array[File?] illumina_se_summaries
    Array[File?] ont_summaries
  }
  Array[File] clearlabs = select_all(clearlabs_summaries)
  Array[File] illumina_pe = select_all(illumina_pe_summaries)
  Array[File] illumina_se = select_all(illumina_se_summaries)
  Array[File] ont = select_all(ont_summaries)
  command <<<
      python3<<CODE
      import json
      header = None
      rows = []
      results = [
          *'~{sep=" " clearlabs}'.split(),
          *'~{sep=" " illumina_pe}'.split(),
          *'~{sep=" " illumina_se}'.split(),
          *'~{sep=" " ont}'.split()
      ]
      with open("theiacov-results.json", 'wt') as results_fh:
          for result in results:
              with open(result, 'rt') as json_fh:
                  row = []
                  json_data = json.load(json_fh)
                  results_fh.write(f'{json.dumps(json_data)}\n')
                  if not header:
                      header = json_data.keys()
                      rows.append("\t".join(header))
                  for col in header:
                      row.append(json_data[col])
                  rows.append("\t".join(row))

      with open("theiacov-results.tsv", "wt") as tsv_fh:
          for row in rows:
              tsv_fh.write(f'{row}\n')
      CODE
  >>>
  output {
    File summaries_tsv  = "theiacov-results.tsv"
    File summaries_json = "theiacov-results.json"
  }
  runtime {
    docker: "python:3.9.5-slim"
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}