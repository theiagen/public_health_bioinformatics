<!--
Thank you for contributing to Theiagen's Public Health Bioinformatics repository! 

Please ensure your contributions are formatted in line with or style guide found here: https://github.com/theiagen/public_health_bioinformatics#contributing-to-the-phb-workflows and follow the instructions with <>  to complete this PR.

As you create the PR, please provide all information that could be required for validation of the workflow.
-->

This PR closes issue #<Issue number>.

🗑️ This dev branch should <NOT> be deleted after merging to main.

### :brain: Aim, Context and Functionality
<!--Please describe the aim of this PR, why the changes were made, and how the workflow should now function -->

### :hammer_and_wrench:  Impacted Workflows/Tasks & Changes Being Made
<!--
-Please use bullet points or headings to describe what is being added or modified to each impacted workflow or task, and the reasoning for those choices. 
-Consider inserting before and after tables or pictures to demonstrate the consequences of the changes on files etc.
-->

### :clipboard: Workflow/Task Steps
This will affect users of the workflow(s) even if they don’t change any workflow inputs <!--  Delete as appropriate -->: Yes/No

Running this workflow on different occasions could result in different results, e.g. due to use of a live database or "latest" docker image <!--  Delete as appropriate -->: Yes/No

#### 🔄 Data Processing Changes
<!-- How are data processed differently through the steps of the task/workflow? Please fill out the sections below.-->
Software or software versions changed: 

Databases or database versions changed:

Data processing/commands changed:

#### ➡️ Modified Inputs
<!--Which inputs of the workflow/task have been added/removed/modified? How have these been modified, e.g input name, type, default parameters, acceptable input ranges etc?-->

#### ⬅️ Modified Outputs
<!--Which outputs of the workflow/task have been added/removed/modified? How have these been modified, e.g. output variable name, output content, output type, file changes?-->

### :test_tube: Testing
#### Test Dataset
<!--Briefly describe what samples were used for testing, e.g. what organism/s, pathogen diversity, etc-->

#### Local Testing
<!--Please show, with screenshots if possible, that your changes pass the local execution of the workflow.
If the whole test dataset was not used, please specify which samples were tested and verify the results were as anticipated.-->

#### Terra Testing
<!--Please show, with screenshots if possible and/or a URL to the job execution, that your changes pass the execution of the workflow on Terra and that all results were as anticipated (including outputs you didn't expect to change!)-->

#### Suggested Scenarios for Reviewer to Test

<!--Please list any potential scenarios that the reviewer should test, including edge cases or data types-->

#### Version Release Testing
<!-- 
-Will changes require functional or validation testing during the release?
-Do new samples need to be added to validation datasets? If so, upload these as a new data table in the validation template workspace with the anticipated results (https://app.terra.bio/#workspaces/theiagen-validations/PHB_Validation_TEMPLATE?) Please describe the new samples here and why these have been chosen.
-Are there any output files that should be checked after running the version release testing?
-->

## :microscope: Final checks
<!--Please delete Yes/No as appropriate to confirm that your changes meet the following quality checks.-->
- The workflow/task has been tested locally and results, including file contents, are as anticipated: Yes/No
- The workflow/task has been tested on Terra and results, including file contents, are as anticipated: Yes/No
- The CI/CD has been adjusted and tests are passing: Yes/No
- Code changes follow the [style guide](https://theiagen.notion.site/Style-Guide-WDL-Workflow-Development-bb456f34322d4f4db699d4029050481c): Yes/No
- Documentation on the "main" version of the Public Health Resources Page has been updated: Yes/No
