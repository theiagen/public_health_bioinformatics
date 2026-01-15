# Getting Started with the Command-Line

!!! dna "What is WDL?"
    Running workflows on the command-line requires the direct use of the WDL (Workflow Development Language). As the name suggests, this is the workflow management language that is used to write and execute workflows. Frank has put together a great video describing 📺 [**WDL Task and Workflow Files**](https://www.youtube.com/watch?v=DNtdra59Y6o) and you can find full instructions below on running these WDL workflows.

## Step 1: Obtain the Workflow and Data

You will need to have access to the WDL workflow file (.wdl) and any associated input files (such as reference genomes, input data files, etc.). To do this, complete the following steps:

### 1. Install Git (if not already installed)

If you don't already have Git installed on your system, you will need to install it. Here's how you can install Git on some common operating systems:

??? toggle "Linux (Ubuntu/Debian)"

    ```bash
    sudo apt update
    sudo apt install git
    ```

??? toggle "macOS"

    Git is usually pre-installed on macOS. However, you can install or update it using Homebrew:
    
    ```bash
    brew install git
    ```

??? toggle "Windows"

    Download and install Git from the official website: <https://git-scm.com/download/win>

### 2. Clone the Repository

1. Open your terminal.
2. Create a directory where you want to store the cloned repository and navigate to it.

    ```bash
    mkdir /path/to/your/desired/new/directory
    cd /path/to/your/desired/new/directory
    ```

3. Clone the <https://github.com/theiagen/public_health_bioinformatics> repository from GitHub using the following command:

    ```bash
    git clone https://github.com/theiagen/public_health_bioinformatics.git
    ```

4. After running the command, Git will download all the repository files and set up a local copy in the directory you specified.

### 3. Navigate to the Cloned Repository

1. Change your working directory to the newly cloned repository:

    ```bash
    cd public_health_bioinformatics
    ```

2. You're now inside the cloned repository's directory. Here, you should find all the files and directories from the GitHub repository.

### 4. Verify the Cloned Repository

You can verify that the repository has been cloned successfully by listing the contents of the current directory using the `ls` (on Linux/macOS) or `dir` (on Windows) command:

```bash
ls
```

This should display the files and directories within the <https://github.com/theiagen/public_health_bioinformatics.git> repository.

Congratulations! You've successfully cloned the <https://github.com/theiagen/public_health_bioinformatics.git> repository from GitHub to your local command-line environment. You're now ready to proceed with running the bioinformatics analysis workflows using WDL as described in subsequent steps.

## Step 2: Install docker and miniWDL

Docker and miniwdl will be required for command-line execution. We will check if these are installed on your system and if not, install them now.

1. Open your terminal.
2. Navigate to the directory where your workflow and input files are located using the `cd` command:
 
    ```bash
    cd /path/to/your/workflow/directory
    ```

3. Check if Docker is installed:

    ```bash
    docker --version
    ```

    If Docker is not installed, follow the official installation guide for your operating system: [**https://docs.docker.com/get-docker/**](https://docs.docker.com/get-docker/)

4. Check if **`miniwdl`** is installed:

    ```bash
    miniwdl --version
    ```

    If **`miniwdl`** is not installed, you can install it using pip:

    ```bash
    pip install miniwdl
    ```

## Step 3: Set up the input.json file for your WDL workflow

In a WDL (Workflow Description Language) workflow, an input JSON file is used to provide attributes (values/files etc) for input variables into the workflow. The names of the input variables must match the names of inputs specified in the workflow file. The workflow files can be found within the git repository that you cloned. Each input variable can have a specific type of attribute, such as String, File, Int, Boolean, Array, etc. Here's a detailed outline of how to specify different types of input variables in an input JSON file:

??? toggle "String Input"
    To specify a string input, use the name of the input variable as the key and provide the corresponding string value. Example:

    ```json
    {
      "sampleName": "VirusSample1",
      "primerSequence": "ACGTGTCAG"
    }
    ```

??? toggle "File Input"
    To specify a file input, provide the path to the input file relative to the directory where you run the `miniwdl` command. Example:

    ```json
    {
      "inputFastq": "data/sample.fastq",
      "referenceGenome": "reference/genome.fasta"
    }
    
    ```

??? toggle "Int Input"
    To specify an integer input, provide the integer value. These do not require quotation marks. Example:

    ```json
    {
      "minReadLength": 50,
      "maxThreads": 8
    }
    ```

??? toggle "Boolean Input"
    To specify a boolean input, use `true` or `false` (lowercase). Example:

    ```json
    {
      "useQualityFiltering": true,
      "useDuplicateRemoval": false
    }
    ```

??? toggle "Array Input"

    To specify an array input, provide the values as an array. Example:

    ```json
    {
      "sampleList": ["Sample1", "Sample2", "Sample3"],
      "thresholds": [0.1, 0.05, 0.01]
    }
    ```

## Step 4: Execute the Workflow

Run the workflow using `miniwdl` with the following command, replacing `your_workflow.wdl` with the actual filename of your WDL workflow and `input.json` with the filename of your input JSON file.

```bash
miniwdl run your_workflow.wdl --input input.json
```

## Step 5: Monitor Workflow Progress

You can monitor the progress of the workflow by checking the console output for updates and log messages. This can help you identify any potential issues or errors during execution.

??? tip "Tips for monitoring your workflow"

    ##### Tips for monitoring workflow progress {% raw %} {#tips-for-monitoring} {% endraw %}

    After you've started the workflow using the **`miniwdl run`** command, you'll see various messages appearing in the terminal. These messages provide information about the various steps of the workflow as they are executed. Monitoring this output is crucial for ensuring that the workflow is progressing as expected.

    The console output will typically show:

    1. **Task Execution:** You will see messages related to the execution of individual tasks defined in your workflow. These messages will include details about the task's name, input values, and progress.
    2. **Logging Information:** Workflow tasks often generate log messages to provide information about what they are doing. These logs might include details about software versions, input data, intermediate results, and more.
    3. **Execution Progress:** The output will indicate which tasks have completed and which ones are currently running. This helps you track the overall progress of the workflow.
    4. **Error Messages:** If there are any errors or issues during task execution, they will be displayed in the console output. These error messages can help you identify problems and troubleshoot them.
    5. **Timing Information:** You might also see timing information for each task, indicating how long they took to execute. This can help you identify tasks that might be taking longer than expected.

    **Example Console Output:**

    Here's an example of what the console output might look like while the workflow is running:

    ```bash

    Running: task1
    Running: task2
    Completed: task1 (Duration: 5s)
    Running: task3
    Error: task2 (Exit Code: 1)
    Running: task4
    ...

    ```

    In this example, you can see that **`task1`** completed successfully in 5 seconds, but **`task2`** encountered an error and exited with a non-zero exit code. This kind of output provides insight into the progress and status of the workflow.

    **What to Look For:**

    As you monitor the console output, pay attention to:

    - **Successful Task Completion:** Look for messages indicating tasks that have completed successfully. This ensures that the workflow is progressing as intended.
    - **Error Messages:** Keep an eye out for any error messages or tasks that exit with non-zero exit codes. These indicate issues that need attention.
    - **Task Order:** The order of task messages can provide insights into the workflow's logic and execution flow.
    - **Timing:** Notice how long each task takes to complete. If a task takes significantly longer than expected, it might indicate a problem.

    **Early Troubleshooting:**

    If you encounter errors or unexpected behavior, the console output can provide valuable information for troubleshooting. You can search for the specific error messages to understand the problem and take appropriate action, such as correcting input values, adjusting parameters, or addressing software dependencies. 

    Monitoring the workflow progress through the console output is an essential practice for successful execution. It allows you to track the status of individual tasks, identify errors, and ensure that your analysis is proceeding as planned. Regularly reviewing the output will help you address any issues and improve the efficiency of your bioinformatics workflow.

??? tip "What to do if you need to cancel a run"
    ##### Canceling a Running Workflow {% raw %} {#canceling-a-run} {% endraw %}

    Canceling a running workflow is an important step in case you need to stop the execution due to errors, unexpected behavior, or any other reason. If you're using `miniwdl` to run your workflow, here's how you can cancel a workflow run while it's in progress:

    1. **Ctrl + C**: The simplest way to cancel a running command in the terminal is to press `Ctrl + C`. This sends an interrupt signal to the running process, which should gracefully terminate it. However, keep in mind that this might not work for all scenarios, and some tasks might not be able to cleanly terminate.
    2. **Terminate Docker Containers**: If your workflow involves Docker containers, you might need to ensure that any Docker containers launched by the workflow are also terminated. To do this, you can manually stop the Docker containers associated with the workflow. You can use the `docker ps` command to list running containers and `docker stop <container_id>` to stop a specific container.
    3. **Kill the miniwdl Process**: If the `Ctrl + C` approach doesn't work, you might need to explicitly kill the `miniwdl` process running in the terminal. To do this, you can use the `kill` command. First, find the process ID (PID) of the `miniwdl` process by running:
        
        ```bash
        ps aux | grep miniwdl
        
        ```
        
        Identify the PID in the output and then run:
        
        ```bash
        kill -9 <PID>
        
        ```
        
        This forcefully terminates the process.
        
    4. **Clean Up Intermediate Files**: Depending on the workflow and how tasks are structured, there might be intermediate files or resources that were generated before the cancellation. You might need to manually clean up these files to free up disk space.
    5. **Check for Workflow-Specific Cancellation**: Some workflows might have specific mechanisms to handle cancellation. Refer to the workflow documentation or user guide to understand if there's a recommended way to cancel the workflow gracefully.
    6. **Check for Any Remaining Resources**: After canceling the workflow, it's a good practice to check for any remaining resources that might need to be cleaned up. This could include temporary files, Docker images, or other resources that were created during the workflow's execution.

    Remember that canceling a workflow might leave the system in an inconsistent state, especially if some tasks were partially executed. After canceling, it's a good idea to review the output and logs to identify any cleanup actions you might need to take.

    It's important to approach workflow cancellation carefully, as abruptly terminating processes can potentially lead to data loss or other unintended consequences. Always make sure you understand the workflow's behavior and any potential side effects of cancellation before proceeding.

## Step 6: Review Output

Once the workflow completes successfully, you will find the output files and results in the designated output directory as defined in your WDL workflow.

??? toggle "Substep 1: Locate the Output Directory"

    Before you begin reviewing outputs, make sure you know where the output directory of your workflow is located. This is typically specified in the workflow configuration or input JSON file. Navigate to this directory using the **`cd`** command in your terminal.
    
    ```bash
    cd /path/to/your/output/directory
    ```

??? toggle "Substep 2: Logs"

    Logs are a valuable source of information about what happened during each step of the workflow. Each task in the workflow might generate its own log file. Here's how to review logs:
    
    1. Use the **`ls`** command to list the files in the output directory:
        
        ```bash
        
        ls
        
        ```
        
    2. Look for log files with names that correspond to the tasks in your workflow. These files often have a **`.log`** extension.
    3. Open a log file using a text editor like **`less`** or **`cat`**:
        
        ```bash
        
        less task_name.log
        
        ```
        
        Use the arrow keys to navigate through the log, and press **`q`** to exit.
        
    4. Inspect the log for messages related to the task's execution, input values, software versions, and any errors or warnings that might have occurred.

??? toggle "Substep 3: stderr (Standard Error) and stdout (Standard Output)"

    stderr and stdout are streams where processes write error messages and standard output, respectively. These are often redirected to files during workflow execution. Here's how to review them:
    
    1. Use the **`ls`** command to list the files in the output directory.
    2. Look for files with names like **`task_name.err`** (for stderr) and **`task_name.out`** (for stdout).
    3. Open the files using a text editor:
        
        ```bash
        
        less task_name.err
        less task_name.out
        
        ```
        
        These files might contain additional information about the task's execution, errors, and output generated during the analysis.

??? toggle "Substep 4: Reviewing Output Files"

    Workflow tasks might generate various types of output files, such as plots, reports, or data files. Here's how to review them:
    
    1. Use the **`ls`** command to list the files in the output directory.
    2. Identify the files generated by your workflow tasks.
    3. Depending on the file type, you can use different tools to open and view them. For example, you might use **`less`** or a text editor for text-based files, or an image viewer for image files.

??? toggle "Substep 5: Interpretation and Troubleshooting"

    As you review the outputs, keep these points in mind:
    
    - **Successful Execution:** Look for indicators of successful task execution, such as expected messages, correct output files, and absence of error messages.
    - **Errors and Warnings:** Pay close attention to any error or warning messages in logs, stderr, or stdout. These can help you identify issues that need troubleshooting.
    - **Input Values and Parameters:** Verify that input values and parameters were correctly passed to tasks. Incorrect input can lead to unexpected behavior.
    - **Software Versions:** Check if the versions of the tools and software used in the workflow match what you expected.
    - **Intermediate Outputs:** Review intermediate outputs generated by tasks. These might provide insights into the workflow's progress and results.

??? toggle "Substep 6: Make Notes and Take Action"

    As you review the outputs, make notes of any issues, errors, or unexpected behavior you encounter. Depending on the severity of the issues, you might need to:
    
    - Adjust input parameters.
    - Re-run specific tasks.
    - Debug and troubleshoot errors.
    - Consult the workflow documentation.
    - Reach out to the Theiagen Genomics bioinformatics experts for assistance. (<support@theiagen.com>)

!!! info "Output Review Conclusion"

    Reviewing the outputs of your bioinformatics workflow is a critical step to ensure the quality of your analysis. Logs, stderr, stdout, and generated output files provide valuable insights into the execution process and results. By carefully reviewing these outputs and addressing any issues, you can enhance the reliability and accuracy of your bioinformatics analysis.

## Step 7: Troubleshooting and Debugging

1. If the workflow encounters errors or fails to execute properly, review the error messages in the terminal.
2. Check for any missing input files, incorrect paths, or issues related to software dependencies.
3. Double-check your input JSON file to ensure that all required inputs are correctly specified.

Congratulations! You have successfully executed a bioinformatics analysis workflow using WDL on the command-line. This tutorial covered the basic steps to run a WDL workflow using the `miniwdl` command-line tool.

Remember that the specific steps and commands might vary depending on the details of your workflow, software versions, and environment. Be sure to consult the documentation for `miniwdl`, WDL, and any other tools you're using for more advanced usage and troubleshooting.

Happy analyzing!
