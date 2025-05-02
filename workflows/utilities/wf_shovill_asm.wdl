version 1.0

import "../../tasks/assembly/task_spades.wdl" as task_spades
import "../../tasks/assembly/task_megahit.wdl" as task_megahit
import "../../tasks/assembly/task_skesa.wdl" as task_skesa
import "../../tasks/alignment/task_bwa.wdl" as task_bwa
import "../../tasks/quality_control/read_filtering/task_pilon.wdl" as task_pilon

workflow shovill_asm {
  input {
    File read1
    File? read2
    String samplename
    String assembler = "spades" # Options: spades, skesa, megahit
    String? kmers
    Boolean use_pilon = false
    String? opts # Extra assembler options

    # Optional parameters
    Int? spades_cpu
    Int? spades_memory
    Int? spades_disk_size
    String? spades_docker

    Int? skesa_cpu
    Int? skesa_memory
    Int? skesa_disk_size
    String? skesa_docker

    Int? megahit_cpu
    Int? megahit_memory
    Int? megahit_disk_size
    String? megahit_docker
    
    # bwa/pilon parameters
    Int? bwa_cpu
    Int? bwa_memory
    Int? bwa_disk_size
    String? bwa_docker

    Int? pilon_cpu
    Int? pilon_memory
    Int? pilon_disk_size
    String? pilon_docker
  }
  if (assembler == "spades") {
    call task_spades.spades {
      input:
        read1 = read1,
        read2 = read2,
        samplename = samplename,
        kmers = kmers,
        spades_opts = opts,
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
        megahit_opts = opts,
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
        skesa_opts = opts,
        cpu = skesa_cpu,
        memory = skesa_memory,
        disk_size = skesa_disk_size,
        docker = skesa_docker
    }
  }
  if (use_pilon) {
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
        cpu = pilon_cpu,
        memory = pilon_memory,
        disk_size = pilon_disk_size,
        docker = pilon_docker
    }
  }
  output {
    File final_assembly = select_first([pilon.assembly_fasta, spades.assembly_fasta, megahit.assembly_fasta, skesa.assembly_fasta])
    File? pilon_changes = pilon.changes
    File? pilon_vcf = pilon.vcf
    String assembler_used = assembler
    String? assembler_version = select_first([spades.spades_version, megahit.megahit_version, skesa.skesa_version])
  }
}