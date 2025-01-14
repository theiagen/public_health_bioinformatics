# Workflow Failure Philosophy

## Our Approach to Workflow Failures

At Theiagen, **we believe our workflows should only fail because of technical issues, not because of poor quality data**. Our goal is to create workflows that can handle data in any condition and still provide meaningful results, _especially_ if that data isn’t perfect.

### What You Can Expect

- **Your workflow will keep running, even with imperfect data**

    _**Data quality shouldn't cause failures**_. Poor or incomplete data should never stop your workflow. Instead, the workflow will process it and (hopefully!) provide meaningful outputs.

- **Your workflow will provide useful feedback, not errors**

    If an issue arises in your data, such as missing or invalid data in a template control, you can know that any _**workflow failures are due to underlying programmatic issues**_, not your data.

- **You’ll gain a better understanding of your data**

    Since poor-quality data will not cause workflow failures, the _**relevant QC results will be available as output**_, so you can understand what's happened and make any needed adjustments moving forward.

### Ongoing Improvements

While we’ve made a lot of progress, we’re still working on fully implementing this philosophy across all of our workflows. If you encounter an issue where poor data quality leads to a failure, please let us know. Your feedback helps us make continuous improvements.

---

**Thanks for being part of the process!** We’re always working to improve and your feedback plays a huge role in making that happen. Together, we’ll keep making things run smoother and easier for everyone.

If you experience a workflow failure related to data quality, we want to hear from you! Please reach out to us at <support@theiagen.com> with the following details:

- The type of data involved
- The error messages or failures encountered
- The steps that led to the issue
    - if this error was generated on the command-line, please include the full command used
    - if this error was generated when running the workflow with Terra.bio, please provide a link to the specific workflow's job history page
