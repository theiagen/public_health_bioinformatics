# PHB Code Contributions

Theiagen Genomics’ [**Public Health Bioinformatics (PHB)**](https://github.com/theiagen/public_health_bioinformatics) workflows are written in [WDL](https://github.com/openwdl/wdl), a language for specifying data processing workflows with a human-readable and writable syntax. Contributions to the workflows contained in the repository are warmly welcomed.

This document gives coding conventions for the WDL code comprising the workflow and task development for PHB. This style guide evolves over time as additional conventions are identified and past conventions are rendered obsolete by changes in the language itself.

Style guide inspired by Scott Frazer’s [WDL Best Practices Style Guide](https://gist.github.com/scottfrazer/aa4ab1945a6a4c331211).

## General Guidelines

***Modularity and Metadata***

- **Best Practice:** Place tasks and workflows in separate files to maintain modularity and clarity.
- **Add a `meta` block** to every task and workflow to provide a brief description of its purpose.

    ```bash
    meta {
      description: "This tool does X"
    }
    ```

***Docker Containers***

- Use a specific Docker container version instead of 'latest' to ensure reproducibility and prevent unexpected changes in container behavior.

    ```bash
    String docker = "quay.io/docker_image:version"
    ```

- Preferentially use containers [`Google's Artifact Registry`](https://console.cloud.google.com/artifacts/docker/general-theiagen/us) rather than those from [`quay.io`](http://quay.io) or [`dockerhub`](https://hub.docker.com/)

***Indentation and Whitespace***

- Use 2-space indentation for all blocks. Avoid using tabs to ensure uniform formatting across editors:

    ```bash
    # perform action
    if [ condition ]; then
      perform_action(variable)
    fi
    ```

- Use a single space when defining variables (`this = that` *not* `this= that` (unless a bash variable where `this=that` is required))

***Bracket and Spacing Conventions***

- Avoid line breaks for opening braces. Keep them on the same line as the declaration. i.e `input {` instead of `input\n{`
  
  ```bash
    # Correct
  input {
    String input_variable
  }

  # Incorrect
  input
  {
    String input_variable
  }
  ```

- Use single space when defining input/output variables & runtime attributes  (`output {` instead of `output{`)
- Separate non-indented constructs (like input and output sections) with a single-line break for readability.

***Command Block Syntax***

- Enclose command blocks in triple angle brackets (<<< ... >>>) for consistency and easier handling of multi-line scripts. It also avoids issues with unescaped special characters in the command block:

    ```bash
    command <<<
      tool --input ~{input} --output ~{output}
    >>>
    ```

## Task Blocks

A WDL task block defines a discrete, reusable step in a workflow. To ensure readability and consistency, follow these conventions when writing task blocks. Include single spaces between the input, command, output, and runtime sections and their enclosing curly brackets.

```bash
task example_task {
  input {

  }
  command <<<
  
  >>>
  output {

  }
  runtime {

  }
}
```

??? toggle "`input` block"
    - The following conventions are used to expose docker, CPU, memory, and disk size

      ```bash
      input {
        String docker = "quay.io/example:1.0.0"  # Docker container for the task
        Int cpu = 4                              # Number of CPUs
        Int memory = 16                          # Memory in GB
        Int disk_size = 100                      # Disk space in GB
      }
      ```

    - If the task accepts additional optional parameters, include an args field.
        ```bash
        input {
          String args = ""
        }
        ```
    - Input and output lists should not be formatted to have the equal sign aligned, but instead use a single space before and after the `=`
        
        ```bash
        correct_output = "output_file"
        long_variable_name = "long_file_name"
        ```
        
    - Expose Docker as an input and runtime variable:
        
        ```bash
        input {
          String docker = "quay.io/example:1.0.0"
        }
        ...
        output {
          String used_docker = docker
        }
        runtime {
          docker: docker
        }
        ```

??? toggle "`command` block"
    - Ensure use of line breaks between different sections of code to improve readability

        ```bash
        # Perform task step 1
        if [ condition ]; then
          action1(variable)
        fi
        
        # Perform task step 2
        if [ another_condition ]; then
          action2(variable)
        fi
        ```
        
    - Split command calls into multiple lines if they have user input variables and/or if the length of the command is very long to avoid text wrapping and/or side-scrolling, e.g.
        - Use backslashes for continuation and indentation to clarify structure:
        
        ```bash
         tool \
          --input ~{input_file} \
          --output ~{output_file} \
          --option1 ~{option1} \
          ...
          --optionN ~{optionN}
        ```
        
    - Add comments that
        - Explain what the optional parameters are
        - Provide links to the tool documentation so future readers of the code know where to find that information
        - Explain what non-intuitive bash/python text wrangling actions do, e.g.
            
            ```bash
            ## awk for gene column ($6) to grab subtype ($15)
            cat ~{file} | awk -F '\t' '{if ($6=="M1") print $15}' > FLU_TYPE
            ```

??? toggle "`output` block"
    - The output block specifies the files or variables produced by the task. Follow these conventions:

      ```bash
      output {
        File result_csv = "output.csv"  # CSV file generated
        File result_log = "log.txt"    # Log file
        String 
      }
      ```
        
    - Ensure the docker container is exposed as an output string, e.g.
        
        ```bash
        input {
          String docker = "us-docker.pkg.dev/general-theiagen/tool:version"
        }
        ...
        output {
          String XX_docker = docker
        }
        runtime {
          docker: docker
        }
        ```

??? toggle "`runtime` block"
    - The runtime block defines the compute resources and environment for the task.
    - Always specify a Docker:

        ```bash
        runtime {
          docker: docker
          cpu: cpu
          memory: memory
          disk: disk_size
        }
        ```

## Workflow Blocks

A WDL workflow block orchestrates the execution of tasks and subworkflows. It defines the inputs, calls tasks or subworkflows, and specifies the final outputs. To ensure readability and consistency, follow these conventions when writing workflow blocks:

### General Guidelines

- Include a block of `import` statements (sorted in alphabetical order).
  - When a workflow imports a task, ensure it is imported under a unique name to avoid conflicts.

```bash
import "../tasks/task_task1.wdl" as task1_task
import "../tasks/task_task2.wdl" as task2_task
```

- A `workflow` block with:

- An `input` section:

  ```bash
  input {
    String input
    String task1_docker = "us-docker.pkg.dev/general-theiagen/tool:version"
    String? task1_optional_argument
  }
  ```

- `call` sections for specified tasks:

  ```bash
  call task1_task.task1 {
    input:
      input = input,
      docker = task1_docker
  }
  ```

- An `output` section:

  - Define all workflow outputs in this section.
  - Use descriptive names for each output variable.

    ```bash
    output {
      # Task 1 outputs
      File task1_out_csv = task1.output_csv
      String task1_version = task1.version

      # Subworkflow outputs
      File subworkflow_out_tsv = subworkflow.task3_out_tsv
      String subworkflow_version = subworkflow.task3_version
    }
    ```


Example formatting is shown below.

??? toggle "wf_example_wf.wdl"

    ```bash
    import "../tasks/task_task1.wdl" as task1_task
    import "../tasks/task_task2.wdl" as task2_task

    import "../workflows/wf_subworkflow.wdl" as subworkflow

    workflow example_wf {
      input {
        String input
        String task1_docker = "us-docker.pkg.dev/general-theiagen/task_1:version"
        String task2_docker = "us-docker.pkg.dev/general-theiagen//task_2:version"
        String? hidden_task3_argument 
        String? hidden_task3_docker
        String? hidden_task4_argument
        String? hidden_task4_docker
      }
      call task1_task.task1 {
        input:
          input = input,
          docker = task1_docker
      }
      call task2_task.task2 {
        input: 
          input = input,
          docker = task2_docker
      }
      call subworkflow.subworkflow {
        input:
          input = input
      }
      output {
        # Task 1 outputs
        File task1_out_csv = task1.output_csv
        String task1_version = task1.version
        String task1_docker = task1.docker
        # Task 2 outputs
        File task2_out_tsv = task2.output_tsv
        String task2_version = task2.version
        String task2_docker = task2.docker
        # Subworkflow outputs
        File task3_out_tsv = subworkflow.task3_out_tsv
        String task3_version = subworkflow.task3_version
        String task3_docker = subworkflow.task3_docker
      }      
    }
    ```


??? toggle "wf_subworkflow.wdl"
    ```bash
    import "../tasks/task_task3.wdl" as task3_task
    import "../tasks/task_task4.wdl" as task4_task

    workflow subworkflow {
      input {
        String input
        
        # optional inputs for tasks inside subworkflows cannot
        #  be seen on Terra, so make them available at the subworkflow
        #  level so they can be modified by a Terra user
        String? task3_argument 
        String? task3_docker
      }
      call task3_task.task3 {
        input:
          input = input,
          args = task3_argument,
          docker = task3_docker
      }
      output {
        File task3_out_tsv = task3.output_tsv
        String task3_version = task3.version
        String task3_docker = task3.docker
      }
    }
    ```

---

??? toggle "`input` section"
    - Optional inputs that should be able to be edited by the user, such as docker containers should be exposed on the workflow level as in the example
    - In the case of subworkflows, all optional inputs should be exposed on the workflow level so that they can be modified by users on Terra

??? toggle "`call` task sections"
    - There should be no blank lines between tasks in workflows

        ```bash
        task A {
        }
        task B {
        }
        ```
        
    - Label a group of outputs by the source/species for organizational purposes when a workflow has many different outputs
        
        ```ebnf
        output {
          ...
          # task99 outputs
          String task99_ouput
          String task99_file
          ...
         }
        ```
