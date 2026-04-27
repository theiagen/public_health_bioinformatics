??? task "`microreact_export`: Data Export"
    This task creates new Microreact projects or updates existing projects using Terra tables as metadata inputs. This is done by accessing Microreact using a user's access token and interacting with Microreact's API to submit created project files directly to Microreact. 

    If a user does not have an access token they will still receive a finalized project file in JSON format. If a submission to Microreact is attempted a response JSON will be returned to the user for review. 

    !!! techdetails "`microreact_export` Technical Details"
        |  | Links |
        | --- | --- |
        | WDL Task | [task_microreact_export.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_export/task_microreact_export.wdl) |