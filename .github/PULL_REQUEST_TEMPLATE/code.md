<!--
Thank you for contributing to Theiagen's Public Health Bioinformatics repository! 

Please ensure your contributions are formatted following our style guide, which can be found here: <https://theiagen.github.io/public_health_bioinformatics/latest/contributing/code_contribution/>.

As you create the PR, please provide any necessary information as suggested in the comments that will help us test your PR.
-->

<!-- Indicate the issue number if applicable; otherwise, delete -->
This PR closes #

🗑️ This dev branch should <NOT> be deleted after merging to main.

## :brain: Summary
<!-- Please summarize what this PR does -->

## :zap: Impacted Workflows/Tasks
<!-- Please list what workflows and/or tasks are impacted by this change -->

This PR may lead to different results in pre-existing outputs: **Yes/No**

This PR uses an element that could cause duplicate runs to have different results: **Yes/No**
<!-- This may be due to using a live database or stochastic data processing. If yes, please describe. -->

## :hammer_and_wrench: Changes
<!-- Describe your changes. -->

### :gear: Algorithm
<!-- Have any changes been made to the algorithm or processing changes under the hood? This can include any changes to the task/workflow algorithm; Docker, software, or database versions; compute resources; etc. If so, please explain. -->

### ➡️ Inputs
<!-- Have any inputs been added or altered? If so, list out the changes. -->

### ⬅️ Outputs
<!-- Have any outputs been added or altered? If so, list out the changes. -->

## :test_tube: Testing
<!-- Please describe how you tested this PR. -->

### Suggested Scenarios for Reviewer to Test
<!-- Please list any potential scenarios that the reviewer should test, including edge cases or data types -->

## :microscope: Final Developer Checklist
<!-- Please mark boxes [X] -->
- [ ] The workflow/task has been tested and results, including file contents, are as anticipated
- [ ] The CI/CD has been adjusted and tests are passing (Theiagen developers)
- [ ] Code changes follow the [style guide](https://theiagen.github.io/public_health_bioinformatics/main/contributing/code_contribution/)
- [ ] Documentation and/or workflow diagrams have been updated if applicable and follow the [documentation style guide](https://theiagen.github.io/public_health_bioinformatics/main/contributing/doc_contribution/)
  - [ ] You have updated the "Last Known Changes" field for any affected workflows in the respective workflow documentation page and for the entry in the `docs/assets/tables/all_workflows.tsv` table to be the tag for the next upcoming release. If you do not know the tag, please put "vX.X.X"

## 🎯 Reviewer Checklist
<!--  Indicate NA when not applicable  -->
- [ ] All changed results have been confirmed
- [ ] You have tested the PR appropriately (see the [testing guide](https://theiagen.notion.site/PR-Testing-Guide-Determining-Appropriate-Levels-of-Testing-4764e98a6aeb460185039c0896714590) for more information)
- [ ] All code adheres to the [style guide](https://theiagen.github.io/public_health_bioinformatics/main/contributing/code_contribution/)
- [ ] MD5 sums have been updated
- [ ] The PR author has addressed all comments
- [ ] The documentation has been updated and adheres to the [documentation style guide](https://theiagen.github.io/public_health_bioinformatics/main/contributing/doc_contribution/)
