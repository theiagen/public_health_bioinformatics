??? task "`spaTyper`: Sequence Typing"
    spaTyper works by identifying the number and order of repeats in the _Staphylococcus_ protein A gene (also known as spa) of a _S. aureus_ assembly. The repeats are assigned a numerical code that is used to assign a spa-type. The tool uses the [Ridom SpaServer Database](http://spaserver2.ridom.de/) in order to assign the appropriate spa-type (please note that this link is typically considered unsafe by modern web browsers).

    Spa-types can be used for accurate and reliable typing of MRSA (methicillin-resistant _Staphylococcus aureus_) strains, and are often used in hospital infection control and epidemiological studies.

    !!! techdetails "spatyper Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_spatyper.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/staphylococcus/task_spatyper.wdl) |
        | Software Source Code | [spaTyper on GitHub](https://github.com/HCGB-IGTP/spaTyper) |
        | Software Documentation | [spaTyper on GitHub](https://github.com/HCGB-IGTP/spaTyper) |
        | Original Publication(s) | _Ridom SpaServer database_: [Typing of methicillin-resistant _Staphylococcus aureus_ in a university hospital setting by using novel software for spa repeat determination and database management](https://doi.org/10.1128/jcm.41.12.5442-5448.2003) |
