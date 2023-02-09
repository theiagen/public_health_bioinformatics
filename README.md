# Public Health Bioinformatics (PHB)
Bioinformatics workflows for characterization, epidemiology and sharing of pathogen genomes. More information about the steps undertaken in these workflows is available via the [Theiagen Public Resources Documentation](https://theiagen.notion.site/Theiagen-Public-Health-Resources-a4bd134b0c5c4fe39870e21029a30566). Support for running these workflows can be sought by raising a [GitHub issue](https://github.com/theiagen/public_health_bioinformatics/issues/new) or by contacting Theiagen at support@theiagen.com.

These workflows are written in [WDL](https://github.com/openwdl/wdl), a language for specifying data processing workflows with a human-readable and writeable syntax. They have been developed by [Theiagen genomics](https://theiagen.com/) to primarily run on the [Terra.bio](https://terra.bio/) platform but can be run locally or on an HPC system at the command-line with Cromwell or miniWDL. 

### Contributors & Influence
* Based on collaborative work with Andrew Lang, PhD & his [Genomic Analysis WDL workflows](https://github.com/AndrewLangvt/genomic_analyses)
* Workflows and task development influenced by The Broad's [Viral Pipes](https://github.com/broadinstitute/viral-pipelines)
* TheiaCoV workflows for genomic characterization influenced by UPHL's [Cecret](https://github.com/UPHL-BioNGS/Cecret) & StaPH-B's [Monroe](https://staph-b.github.io/staphb_toolkit/workflow_docs/monroe/)
* The TheiaCoV workflow user community. To provide feedback, please raise a GitHub issue [GitHub issue](https://github.com/theiagen/public_health_viral_genomics/issues/new).

### Contributing to the PHVG workflows
Contributions to the workflows contained in this repository are warmly welcomed. Our style guide may be found below for convenience of formatting.

2-space indents (no tabs), no line break for openning braces, single space when defining input/output variables & runtime attributes, single-line breaks between non-intended constructs, and task commands enclosed with triple braces (`<<< ... >>>`). 

<em>E.g.</em>:
```
workflow w {
  input {
    String input
  }
  call task_01 {
    input:
      input = input
  }
  call task_02 {
    input: 
      input = input
  }
  output {
    File task_01_out = task_01.output
    File task_02_out = task_02.output 
  }      
}

task task1 {
  input {
    String input
    String docker = "theiagen/utility:1.1"
  }
  command <<<
    echo '~{input}' > output.txt
  >>>
  output {
    File output = "output.txt"
  }
  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 2
    disks "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 0
  }
}

task task_02 {
  input{
    String input
    String docker = "theiagen/utility:1.1"
  }
  command <<<
    echo '~{input}' > output.txt
  >>>
  output {
    File output = "output.txt"
  }
  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 2
    disks "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 0
  }
}
```
Style guide inspired by [scottfrazer](https://gist.github.com/scottfrazer)'s [WDL Best Pratcices Style Guide](https://gist.github.com/scottfrazer/aa4ab1945a6a4c331211)
