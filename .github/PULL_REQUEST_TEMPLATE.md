<!--
Thank you for contributing to Theiagen's Public Health Bioinformatics repository!
Please fill in the appropriate checklist below and delete whatever is not relevant.

Please replace all '[ ]' with '[X]' to demonstrate completion.

Please amend all text within '<>' as appropriate. 

Documentation on how to contribute can be found at https://github.com/theiagen/public_health_bioinformatics#contributing-to-the-phb-workflows
-->



This PR closes issue <Issue #>.

🗑️ This dev branch should <NOT> be deleted after merging to main.

### :brain: Context, Functionality and Rationale
<!--Please describe the aim of this PR, what changes have been made to workflow functionality, and the rationale for this-->

### :hammer_and_wrench:  Impacted Workflows/Tasks & Changes Being Made
<!--
-Please use bullet points or headings to describe what is being added or modified to each impacted workflow and task. 
-Consider inserting before and after pictures or tables to demonstrate the consequences of the changes on files etc
-->

### 🧑‍🔬 Expected Workflow Changes for Users
<!--
-Will this affect users of the workflow(s) even if they don’t change any workflow inputs?
-->

### :clipboard: Workflow/Task Steps
<!--What are the main steps of your workflow/task? 
Any trade-offs for decisions made?-->

#### ➡️ Modified Inputs
<!--Which inputs of the workflow/task have been added/removed/modified? How have these been modified, e.g input name, type, default parameters, acceptable input ranges etc?-->

#### ⬅️ Modified Outputs
<!--Which outputs of the workflow/task have been added/removed/modified? How have these been modified, e.g output variable name, output content, output type, file changes?-->

#### 🔄 Modified  Actions
<!-- 
-How are data processed (differently) through the steps of the task/workflow? Make it explicit enough so that someone who doesn't have deep knowledge of the workflow/task can understand how the rationale was implemented
-List any changes being made to tool or database workflow components, including version changes. -->

### :test_tube: Testing
#### Test Dataset
<!--Briefly describe what samples were used for testing, e.g. what organism/s, pathogen diversity, etc-->

#### Local Testing
<!--Please show, with screenshots if possible, that your changes pass the local execution of the workflow.
If the whole test dataset was not used, please specify which samples were tested and verify the results were as anticipated.-->

#### Terra Testing
<!--Please show, with screenshots if possible and/or a URL to the job execution, that your changes pass the execution of the workflow on Terra and that the results were as anticipated-->

#### Suggested Scenarios for Reviewer to Test

<!--Please list any potential scenarios that the reviewer should test, including edge cases or data types-->

#### Version Release Testing
<!-- 
-Will changes require functional or validation testing during the release?
-Do new samples need to be added to validation datasets? If so, please specify which samples should be added
-Are there any output files that should be checked after running the version release testing?
-->

## :microscope: Quality checks
<!--Please check the boxes [X] to confirm that your changes meet the following quality checks.-->
- [ ] The workflow/task has been tested locally and results, including file contents, are as anticipated.
- [ ] The workflow/task has been tested on Terra and results, including file contents, are as anticipated.
- [ ] The CI/CD has been adjusted and tests are passing
- [ ] Code changes follow the [style guide](https://theiagen.notion.site/Style-Guide-WDL-Workflow-Development-bb456f34322d4f4db699d4029050481c)
- [ ] Documentation on the "main" version of the Public Health Resources Page has been updated
