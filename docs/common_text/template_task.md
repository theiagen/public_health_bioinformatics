??? task "`<tool_name>`: <brief description\>"

    _Provide sufficient information about your task here so that users can understand how to use the task, what actions it performs, and any important considerations or requirements. Use admonitions as appropriate to highlight important information, tips, or technical details._

    _Be sure any links within this file are written relative to the destination file. For example, if this file is in `docs/common_text/` but is being included in `docs/workflows/standalone/` and you want to link to a file in `docs/assets/figures/`, **the link should be written as `[link_name](../../assets/figures/link_destination.md)`**._

    _If your task is used slighly differently in different contexts, you can provide conditionals inside of comments._

    _To include any conditional information in the destination page, add the `condition="condition_name"` parameter to the `include_md` macro call on the destination page. Please note that the `if: <condition_name>` and `endif` syntax within the markdown comment (`<!-- -->`) is required for a conditional statement to work._ 

    _Here is an example:_
    
    <!-- if: condition -->

    `<!-- if: <condition_name> -->`
    
    !!! tip "Conditional Content"
        _This is content that is only shown if <condition_name\> is provided. See the `kraken_task.md` common_text file for examples on conditional usage._ 
    
    `<!-- endif -->`
        
    <!-- endif -->
  
    !!! techdetails "<tool_name> Technical Details"        
        _This section is required for all tasks. If the Software Source Code, Software Documentation, and/or Original Publication(s) fields are not applicable, they may be removed._

        |  | Links |
        | --- | --- |
        | Task | [<link_to_task_wdl>.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/<link_to_task_wdl>.wdl) |
        | Software Source Code | [<tool name>](<link to source code>) |
        | Software Documentation |  [<tool name>](<link to source code documentation>) |
        | Original Publication(s) | [<publication name>](<link to publication>) |
