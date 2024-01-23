<!--
Thank you for contributing to Theiagen's Public Health Bioinformatics repository! 

Please ensure your contributions are formatted in line with or style guide found here: https://github.com/theiagen/public_health_bioinformatics#contributing-to-the-phb-workflows and follow the instructions with <>  to complete this PR.

As you create the PR, please provide all information required to validate the workflow.
-->

This PR closes #<Issue number>.

🗑️ This dev branch should <NOT> be deleted after merging to main.

## :brain: Aim, Context and Functionality
<!--Please describe the aim of this PR, why the changes were made, and how the workflow should now function -->

## :hammer_and_wrench:  Impacted Workflows/Tasks & Changes Being Made
This will affect the behavior of the workflow(s) even if users don’t change any workflow inputs relative to the last version <!--  Delete as appropriate -->: Yes/No

Running this workflow on different occasions could result in different results, e.g. due to use of a live database, "latest" docker image, or stochastic data processing <!--  Delete as appropriate. If yes, please describe. -->: Yes/No 
<!--
-Please use bullet points or headings to describe what is being added or modified to each impacted workflow or task, and the reasoning for those choices. 
-Consider inserting before and after tables or pictures to demonstrate the consequences of the changes on files etc.
-->

## :clipboard: Workflow/Task Step Changes

#### 🔄 Data Processing 
<!-- How are data processed differently through the steps of the task/workflow? 
Please describe in the sections below. 
If nothing has changed, please explicitly say so.-->

Docker/software or software versions changed: 

Databases or database versions changed:

Data processing/commands changed:

File processing changed:

Compute resources changed:

#### ➡️ Inputs 
<!--Which inputs of the workflow/task have been added/removed/modified? 
How have these been modified, e.g input name, type, default parameters, acceptable input ranges etc? 
If nothing has changed, please explicitly say so.-->

#### ⬅️ Outputs 
<!--Which outputs of the workflow/task have been added/removed/modified? 
How have these been modified, e.g. output variable name, output content, output type, file changes? 
If nothing has changed, please explicitly say so.-->

## :test_tube: Testing 
#### Test Dataset
<!--Briefly describe what samples were used for testing, e.g. what organism/s, pathogen diversity, etc. -->

#### Commandline Testing with MiniWDL or Cromwell (optional)
<!--
Please show, with screenshots if possible, that your changes pass the local execution of the workflow.
If the whole test dataset was not used, please specify which samples were tested and verify the results were as anticipated. 
If local testing was not undertaken/possible, please explicitly state this.-->

#### Terra Testing
<!--Please show, with screenshots if possible and/or a URL to the job execution, that your changes pass the execution of the workflow on Terra and that all results were as anticipated (including outputs you didn't expect to change!)-->

#### Suggested Scenarios for Reviewer to Test
<!--Please list any potential scenarios that the reviewer should test, including edge cases or data types-->

#### Theiagen Version Release Testing (optional)
<!-- 
-Will changes require functional or validation testing (checking outputs etc) during the release?
-Do new samples need to be added to validation datasets? If so, upload these to the appropriate validation workspace Google bucket (). Please describe the new samples here and why these have been chosen.
-Are there any output files that should be checked after running the version release testing?
-->

## :microscope: Final Developer Checklist
<!--Please mark boxes [X] -->
- [ ] The workflow/task has been tested locally and results, including file contents, are as anticipated
- [ ] The workflow/task has been tested on Terra and results, including file contents, are as anticipated
- [ ] The CI/CD has been adjusted and tests are passing (to be completed by Theiagen developer)
- [ ] Code changes follow the [style guide](https://theiagen.notion.site/Style-Guide-WDL-Workflow-Development-bb456f34322d4f4db699d4029050481c)


## 🎯 Reviewer Checklist 
<!--  Indicate NA when not applicable  -->
- [ ] All impacted workflows/tasks have been tested on Terra with a different dataset than used for development
- [ ] All reviewer-suggested scenarios have been tested and any additional
- [ ] All changed results have been confirmed to be accurate
- [ ] All workflows/tasks impacted by change/s have been tested using a standard validation dataset to ensure no unintended change of functionality
- [ ] All code adheres to the style guide
- [ ] MD5 sums have been updated
- [ ] The PR author has addressed all comments

## 🗂️ Associated Documentation (to be completed by Theiagen developer)
<!--  Indicate NA when not applicable -->
- [ ] Relevant documentation on the Public Health Resources "PHB Main" has been updated
- [ ] Workflow diagrams have been updated to reflect changes
