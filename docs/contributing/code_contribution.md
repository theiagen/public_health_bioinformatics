# PHB Code Contributions

Theiagen Genomics’ [**Public Health Bioinformatics (PHB)**](https://github.com/theiagen/public_health_bioinformatics) workflows are written in [WDL](https://github.com/openwdl/wdl), a language for specifying data processing workflows with a human-readable and writable syntax. Contributions to the workflows contained in the repository are warmly welcomed.

This document gives coding conventions for the WDL code comprising the workflow and task development for PHB. This style guide evolves over time as additional conventions are identified and past conventions are rendered obsolete by changes in the language itself.

Style guide inspired by Scott Frazer’s [WDL Best Practices Style Guide](https://gist.github.com/scottfrazer/aa4ab1945a6a4c331211).

## General Guidelines

- Put tasks and workflows in separate files in the appropriate folders.
- Always add a description as metadata

    ```bash
    meta {
      description: "This tool does X"
    }
    ```

- Ensure that the docker container is locked to a given version, not `latest`

    ```bash
    String docker = "quay.io/docker_image:version"
    ```

- Preferentially use containers [`Google's Artifact Registry`](https://console.cloud.google.com/artifacts/docker/general-theiagen/us) rather than those from [`quay.io`](http://quay.io) or [`dockerhub`](https://hub.docker.com/)
- Use 2-space indents (no tabs)

    ```bash
    # perform action
    if [ this ]; then
      action1(variable)
    fi
    ```

- Do not use line break for opening braces
- Use single space when defining input/output variables & runtime attributes  (`output {` instead of `output{`)
- Use single-line breaks between non-intended constructs
- Enclose task commands with triple angle brackets (`<<< ... >>>`)
- Consistently use white space with variables (`this = that` *not* `this= that` (unless a bash variable where `this=that` is required))

## Task Blocks

The task should contain the following sections. Include _single_ spaces between input, command, output, and runtime closing and opening curly brackets.

```bash
input {

}
command <<<

>>>
output {

}
runtime {

}
```

??? toggle "`input` block"
    - The following conventions are used to expose docker, CPU, memory, and disk size

        ```bash
        input {
          String docker = "..."
        	Int cpu = x
        	Int memory = y
        	Int disk_size = z
        }
        ```
        
    - If additional arguments should be allowed to be passed to the task, this input should follow the convention below:
        
        ```bash
        input {
          String args = ""
        }
        ```
        
    - Input and output lists should not be formatted to have the equal sign aligned, but instead use a single space before and after the `=`
        
        ```bash
        output1_x = string1
        output2_that_does_y = string2
        ```
        
    - Ensure the docker container is exposed as an input and as an output string
        
        ```bash
        input {
        	String docker = ""
        }
        ...
        output {
            String XX_docker = docker
        }
        runtime {
        	docker: docker
        }
        ```

??? toggle "`command` block"
    - Ensure use of line breaks between different sections of code to improve readability
 
        ```bash
        # if this, perform action 1
        if [ this ]; then
          action1(variable)
        fi
        
        # if that, perform action 2
        if [ that ]; then
          action2(variable)
        fi
        ```
        
    - Split command calls into multiple lines if they have user input variables and/or if the length of the command is very long to avoid text wrapping and/or side-scrolling, e.g.
        - Use indentation as appropriate
        
        ```bash
         tool \
          --option1 ~{option1} \
          --option2 ~{option2} \
          ...
          --option999 ~{option999}
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
    - File types should be clearly stated in the output name variables

        ```bash
        output1_csv = file1.csv
        output2_tsv = file2.tsv
        ```
        
    - Ensure the docker container is exposed as an output string, e.g.
        
        ```bash
        input {
        	String docker
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
    - Always use a docker container

## Workflow Blocks

The workflow/sub-workflow file should contain:

- a block of `import` statements (alphabetical order),
    - When a workflow imports a task, make sure that it is imported under a different name than the task it is calling
- a `workflow` block with
    - an `input` section
    - `call` sections for specified tasks
    - an `output` section

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
