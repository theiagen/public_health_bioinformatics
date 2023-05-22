<!--
Thank your for contributing to Theiagen's Public Health Bioinformatics repository!
Please fill in the appropriate checklist below and delete whatever is not relevant.

Documentation on how to contribute can be found at https://github.com/theiagen/public_health_bioinformatics#contributing-to-the-phb-workflows

Please replace all '[ ]' with '[X]' to demonstrate completion.

Please delete all text within '<>'. 
-->
Closes <Issue number>

## :hammer_and_wrench: Changes being made

<!--Here give examples of the changes you've made to this pull request. Include an itemized list if you can. It'll help the reviewer-->

## :brain: Context and Rationale

<!--What's the context of the changes? Where there any trade-offs you had to consider?-->

## :test_tube: Testing

### Locally

<!--Please show, with screenshots when possible, that your changes pass the local execution of the workflow-->

### Terra

<!--Please show, with screenshots when possible and/or a URL to the job execution, that your changes pass the execution of the workflow on Terra.bio-->

## :microscope: Quality checks

<!--Please ensure that your changes respect the following quality checks.-->

Pull Request (PR) checklist:
- [ ] Include a description of what is in this pull request in this message.
- [ ] If it's a task or a workflow, it includes a descriptive section

```
meta {
    description: "This tool does X"
  }
```

- [ ] Ensure all docker images are using quay.io if possible
- [ ] The docker image is modifiable, and saved as a String output

```
input {
	String docker = "quay.io/docker_image:version"
}
...
output {
    String XX_docker = docker
}
runtime {
	docker: docker
}
```

- [ ] Everything follows the style guide
  - [ ] Consistent usage of white space with variables (`this = that` *not* `this= that` (unless a bash variable where `this=that` is required))
  - [ ] Consistent use of brackets (`output {` instead of `output{`)
  - [ ] Ensure line breaks between different sections of code to improve readability

  ```
  # perform action
  if [ this ]; then
    action1(variable)
  fi

  # perform action 2
  if [ that ]; then
    action2(variable)
  fi
  ```
  
  - [] Indentation to be two spaces from the parent line

  ```
  # perform action
  if [ this ]; then
    action1(variable)
  fi
  ```

  - [ ] Add comments in tasks to explain what nonintuitive bash text wrangling actions do

  ```
  ## awk for gene column ($6) to grab subtype ($15)
  cat ~{file} | awk -F '\t' '{if ($6=="M1") print $15}' > FLU_TYPE
  ```

  - [ ] Add comments in tasks in the `command {}` block that either (a) explain what the optional parameters are or (b) provide links to the tool documentation so future readers of the code know where to find that information.
  - [ ] Split command calls into multiple lines if has user input variables and/or if the length of the command is very long to avoid text wrapping and/or side-scrolling.

  ```
  tool \
  --option1 ~{option1} \
  --option2 ~{option2} \
  ...
  --option999 ~{option999}
  ```

  NOT

  ```
  tool --option1 ~{option1} --option2 ~{option2} ... --option999 ~{option999}
  ```

  - [ ] No blank lines between tasks in sub-workflows

  ```
  task A {
  }
  task B {
  }
  ```

  - [ ] Input and output list should not be formatted to have the equal sign aligned

  ```
  output1 = file1
  output2_that_does_x = string2 
  ```

  - [ ] Include spaces between input, command, output and runtime closing and opening curly brackets

  ```
  input {
  
  }
  command <<<
  
  >>>
  output {

  }
  runtime {

  }
  ```