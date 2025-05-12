version 1.0

import "../../tasks/assembly/task_spades.wdl" as task_spades
import "../../tasks/assembly/task_megahit.wdl" as task_megahit
import "../../tasks/assembly/task_skesa.wdl" as task_skesa
import "../../tasks/alignment/task_bwa.wdl" as task_bwa
import "../../tasks/quality_control/read_filtering/task_pilon.wdl" as task_pilon
import "../../tasks/quality_control/read_filtering/task_filter_contigs.wdl" as task_filter_contigs

workflow digger_denovo {
  input {
    File read1
    File? read2
    String samplename
    String assembler = "skesa" # Options: spades, skesa, megahit
    Int min_contig_length = 200
    String? kmers
    Boolean call_pilon = false
    String? assembler_options # Extra assembler options
    Boolean run_filter_contigs = true # Default: Filter contigs after assembly
    # Optional parameters for spades
    String? spades_type = "isolate"
    Int? spades_cpu
    Int? spades_memory
    Int? spades_disk_size
    String? spades_docker
    # Optional parameters for skesa
    Int? skesa_cpu
    Int? skesa_memory
    Int? skesa_disk_size
    String? skesa_docker
    # Optional parameters for megahit
    Int? megahit_cpu
    Int? megahit_memory
    Int? megahit_disk_size
    String? megahit_docker
    # Optional parameters for bwa
    Int? bwa_cpu
    Int? bwa_memory
    Int? bwa_disk_size
    String? bwa_docker
    # Optional parameters for pilon
    Int? pilon_cpu
    Int? pilon_memory
    Int? pilon_disk_size
    String? pilon_docker
    Int? pilon_min_mapping_quality = 60 # Shovill default
    Int? pilon_min_base_quality = 3 # Shovill default
    Float? pilon_min_depth = 0.25 # Shovill default
    String? pilon_fix = "bases" # Options: all, snps, indels, gaps
    # Optional parameters for filtering 
    Int filter_contigs_min_length = 200 # Default we set before
    Float filter_contigs_min_coverage = 2.0 # Default we set before
    Boolean filter_contigs_skip_length_filter = false
    Boolean filter_contigs_skip_coverage_filter = false
    Boolean filter_contigs_skip_homopolymer_filter = false
    Int? filter_contigs_cpu
    Int? filter_contigs_memory
    Int? filter_contigs_disk_size
    String? filter_contigs_docker
  }
  if (assembler == "spades") {
    call task_spades.spades {
      input:
        read1 = read1,
        read2 = read2,
        samplename = samplename,
        kmers = kmers,
        spades_type = spades_type,
        spades_opts = assembly_options,
        cpu = spades_cpu,
        memory = spades_memory,
        disk_size = spades_disk_size,
        docker = spades_docker
    }
  }
  if (assembler == "megahit") {
    call task_megahit.megahit {
      input:
        read1 = read1,
        read2 = read2,
        samplename = samplename,
        kmers = kmers,
        min_contig_length = min_contig_length,
        megahit_opts = assembly_options,
        cpu = megahit_cpu,
        memory = megahit_memory,
        disk_size = megahit_disk_size,
        docker = megahit_docker
    }
  }
  if (assembler == "skesa") {
    call task_skesa.skesa {
      input:
        read1 = read1,
        read2 = read2,
        samplename = samplename,
        min_contig_length = min_contig_length,
        skesa_opts = assembly_options,
        cpu = skesa_cpu,
        memory = skesa_memory,
        disk_size = skesa_disk_size,
        docker = skesa_docker
    }
  }
  if (call_pilon) {
    call task_bwa.bwa {
      input:
        read1 = read1,
        read2 = read2,
        samplename = samplename,
        reference_genome = select_first([spades.assembly_fasta, megahit.assembly_fasta, skesa.assembly_fasta]),
        cpu = bwa_cpu,
        memory = bwa_memory,
        disk_size = bwa_disk_size,
        docker = bwa_docker
    }
    call task_pilon.pilon {
      input:
        assembly = select_first([spades.assembly_fasta, megahit.assembly_fasta, skesa.assembly_fasta]),
        bam = bwa.sorted_bam,
        bai = bwa.sorted_bai,
        samplename = samplename,
        min_mapping_quality = pilon_min_mapping_quality,
        min_base_quality = pilon_min_base_quality,
        min_depth = pilon_min_depth,
        fix = pilon_fix,
        cpu = pilon_cpu,
        memory = pilon_memory,
        disk_size = pilon_disk_size,
        docker = pilon_docker
    }
  }
  if (run_filter_contigs) {
    call task_filter_contigs.filter_contigs {
      input:
        samplename = samplename,
        assembly_fasta = select_first([pilon.assembly_fasta, spades.assembly_fasta, megahit.assembly_fasta, skesa.assembly_fasta]),
        min_length = filter_contigs_min_length,
        min_coverage = filter_contigs_min_coverage,
        skip_length_filter = filter_contigs_skip_length_filter,
        skip_coverage_filter = filter_contigs_skip_coverage_filter,
        skip_homopolymer_filter = filter_contigs_skip_homopolymer_filter,
        cpu = filter_contigs_cpu,
        memory = filter_contigs_memory,
        disk_size = filter_contigs_disk_size,
        docker = filter_contigs_docker
      }
    }
  output {
    File assembly_fasta = select_first([filter_contigs.filtered_fasta, pilon.assembly_fasta, spades.assembly_fasta, megahit.assembly_fasta, skesa.assembly_fasta])
    File? contigs_gfa = spades.assembly_gfa
    File? filtered_contigs_metrics = filter_contigs.assembly_filtering_metrics
    File? pilon_changes = pilon.changes
    File? pilon_vcf = pilon.vcf
    String assembler_used = assembler
    String? assembler_version = select_first([spades.spades_version, megahit.megahit_version, skesa.skesa_version])
  }
}