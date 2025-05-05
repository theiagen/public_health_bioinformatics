# PHB Documentation Contribution Guide

The documentation for PHB is hosted in the `docs/` directory. This documentation is written in Markdown and is built using [MkDocs](https://www.mkdocs.org/) and the [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) theme.

This guide is intended to provide a brief overview of the documentation structure and how to contribute to the documentation, including standard language and formatting conventions.

## Local Installation & Live Previews

Since the documentation is built off of the `main` branch, it is highly recommended to preview your changes before making a PR. You can do this by installing the necessary packages and previewing ("serving") the documentation locally.

To test your documentation changes, you will need to have the following packages installed on your local VM:

```bash
pip install mkdocs-material mkdocs-material-extensions mkdocs-git-revision-date-localized-plugin mike mkdocs-glightbox mkdocs-macros-plugin
```

Once installed, navigate to the top directory in PHB. The live preview server can be activated by running the following command:

```bash
mkdocs serve
```

This will prompt you to open your browser to the appropriate local host address (by default, localhost:8000). Every time you save a change, the documentation will automatically update in the browser.

### VSCode Extensions

Here are some VSCode Extensions can help you write and edit your markdown files (and allow you preview changes without running the server, though formatting will suffer):

- [Markdown Preview Enhanced (Yiyi Wang)](https://marketplace.visualstudio.com/items?itemName=shd101wyy.markdown-preview-enhanced) - This extension is good for previewing markdown files in VSCode, but is **not** good at rendering any of the more advanced features such as callouts or tables.
- [Markdown All in One (Yu Zhang)](https://marketplace.visualstudio.com/items?itemName=yzhang.markdown-all-in-one) - This extension allows you to use regular word-processing short-cuts to format your markdown files, like Ctrl-B to bold text, Ctrl-I for italics without having to manually type the `**` or `_` characters.
- [markdownlint (David Anson)](https://marketplace.visualstudio.com/items?itemName=DavidAnson.vscode-markdownlint) - This extension will help you catch any formatting errors in your markdown files.

### Helpful Websites

- [Excel to Markdown Table](https://tableconvert.com/excel-to-markdown) - This website will convert an Excel table into markdown format, which can be copied and pasted into your markdown file.
- [Material for MkDocs Reference](https://squidfunk.github.io/mkdocs-material/reference/) - This is the official reference for the Material for MkDocs theme, which will help you understand how to use the theme's features.
- [Dead Link Check](https://www.deadlinkchecker.com/) - This website will scan your website to ensure that all links are working correctly. This will only work on the deployed version of the documentation, not the local version.

## Standard Language & Formatting Conventions

In order to maintain cohesive documentation, the following language and formatting conventions should be followed:

### Language Conventions

The following language conventions should be followed when writing documentation:

- The documentation should be written in American English (sorry to our friends across the pond!)
- **The following variables should receive the following descriptions**:
    - `cpu` - Number of CPUs to allocate to the task
    - `disk_size` - Amount of storage (in GB) to allocate to the task
    - `docker` or `docker_image` - The Docker container to use for the task
    - `memory` - Amount of memory/RAM (in GB) to allocate to the task

### Formatting Conventions

- **Bold Text** - Use `**bold text**` to indicate text that should be bolded.
- _Italicized Text_ - Use `_italicized text_` to indicate text that should be italicized.
- ==Highlighted Text== - Use `==highlighted text==` to indicate text that should be highlighted.
- `Code` - Use ````code` ``` (backticks) to indicate text that should be formatted as code.
- ^^Underlined Text^^ - Use `^^underlined text^^` to indicate text that should be underlined (works with our theme; not all Markdown renderers support this).
- > Citations
    - Use a `>` to activate quote formatting for a citation. Make sure to separate multiple citations with a comment line (`<!-- -->`) to prevent the citations from running together.
    - Use a reputable citation style (e.g., Vancouver, Nature, etc.) for all citations.
- Callouts/Admonitions - These features are called "call-outs" in Notion, but are "Admonitions" in MkDocs. [I highly recommend referring to the Material for MkDocs documentation page on Admonitions to learn how best to use this feature](https://squidfunk.github.io/mkdocs-material/reference/admonitions/). Use the following syntax to create a callout:

    ```markdown
    !!! note
        This is a note. Observe I am indented with four spaces.
    ```

    Please see the [Admonition documentation](https://squidfunk.github.io/mkdocs-material/reference/admonitions/) for more information on how to change the title, enable toggles, and more.

    The following custom callout types are supported _in addition to the standard admonitions supported by our theme_ [more information here](https://squidfunk.github.io/mkdocs-material/reference/admonitions/#supported-types):
  
    !!! dna
        This is a DNA admonition. Admire the cute green DNA emoji. You can create this with the `!!! dna` syntax.

        Use this admonition when wanting to convey general information or highlight specific facts.

    ???+ toggle
        This is a toggle-able section. The emoji is an arrow pointing to the right downward. You can create this with the `??? toggle` syntax. I have added a `+` at the end of the question marks to make it open by default.

        Use this admonition when wanting to provide additional _optional_ information or details that are not strictly necessary, or take up a lot of space.

    ???+ task
        This is a toggle-able section **for a workflow task**. The emoji is a gear. Use the `??? task` syntax to create this admonition. Use `!!! task` if you want to have it be permanently expanded. I have add a `+` at the end of the question marks to make this admonition open by default and still enable its collapse.

        Use this admonition when providing details on a workflow, task, or tool.

    !!! caption
        This is a caption. The emoji is a painting. You can create this with the `!!! caption` syntax. A caption can be added beneath the picture and will also look nice.

        Use this admonition when including images or diagrams in the documentation.

    !!! techdetails
        This is where you will put technical details for a workflow task. You can create this by `!!! techdetails` syntax.

        Use this admonition when providing technical details for a workflow task or tool. These admonitions should include the following table:

        |  | Links |
        | --- | --- |
        | Task | [link to the task file in the PHB repository on GitHub] |
        | Software Source Code | [link to tool's source code] |
        | Software Documentation | [link to tool's documentation] |
        | Original Publication(s) | [link to tool's publication] |

        If any of these items are unfillable, delete the row.

- Images - Use the following syntax to insert an image:

    ```markdown
    !!! caption "Image Title"
        ![Alt Text](/path/to/image.png)
    ```

- Indentation - **_FOUR_** spaces are required instead of the typical two. This is a side effect of using this theme. If you use two spaces, the list and/or indentations will not render correctly. This will make your linter sad :(

    ```markdown
    - first item
        - second item
            - third item
    ```

- Tables - Use the following syntax to create a table

    ```markdown
    | Header 1 | Header 2 | Header 3 |
    |---|---|---|
    | value 1 | value2 | value3 |
    ```

    Note that this is not a "pretty" markdown table. This is because the spacing would be crazy in the markdown file, especially for tables with a lot of text and/or columns. The table will render correctly in the documentation.

- Links - Use the following syntax to create a link. This is works for both files and websites. If linking a file, use the relative path.

    ```markdown
    [Link Text](https://www.example.com)
    ```

- End all pages with an empty line

## Documentation Structure

A brief description of the documentation structure is as follows:

- `docs/` - Contains the Markdown files for the documentation.
    - `assets/` - Contains images and other files used in the documentation.
        - `figures/` - Contains images, figures, and workflow diagrams used in the documentation. For workflows that contain many images (such as BaseSpace_Fetch), it is recommended to create a subdirectory for the workflow.
        - `files/` - Contains files that are used in the documentation. This may include example outputs or templates. For workflows that contain many files (such as TheiaValidate), it is recommended to create a subdirectory for the workflow.
        - `logos/` - Contains Theiagen logos and symbols used in the documentation.
        - `metadata_formatters/` - Contains the most up-to-date metadata formatters for our submission workflows.
        - `new_workflow_template.md` - A template for adding a new workflow page to the documentation. You can see this template [here](../assets/new_workflow_template.md)
    - `contributing/` - Contains the Markdown files for our contribution guides, such as this file
    - `javascripts/` - Contains JavaScript files used in the documentation.
        - `tablesort.js` - A JavaScript file used to enable table sorting in the documentation.
    - `overrides/` - Contains HTMLs used to override theme defaults
        - `main.html` - Contains the HTML used to display a warning when the latest version is not selected
    - `stylesheets/` - Contains CSS files used in the documentation.
        - `extra.css` - A custom CSS file used to style the documentation; contains all custom theme elements (scrollable tables, resizable columns, Theiagen colors), and custom admonitions.
    - `workflows/` - Contains the Markdown files for each workflow, organized into subdirectories by workflow category
    - `workflows_overview/` - Contains the Markdown files for the overview tables for each display type: alphabetically, by applicable kingdom, and by workflow type.
    - `index.md` - The home/landing page for our documentation.

### Adding a Page for a New Workflow {% raw %} {#new-page} {% endraw %}

!!! tip "Hey, we've got a template for that!"
    Please see our template [here](../assets/new_workflow_template.md) for ease of use. Please remove all italicized text and replace with the appropriate information. If in doubt, please refer to existing documentation.

If you are adding a new workflow, there are a number of things to do in order to include the page in the documentation:

1. Add a page with the title of the workflow to appropriate subdirectory in `docs/workflows/`. Please use [the template](../assets/new_workflow_template.md) found in the `assets/` folder.
2. Collect the following information for your new workflow:
     - Workflow Name - Link the name with a relative path to the workflow page in appropriate `docs/workflows/` subdirectory
     - Workflow Description - Brief description of the workflow
     - Applicable Kingdom - Options: "Any taxa", "Bacteria", "Mycotics", "Viral"
     - Workflow Level (_on Terra_) - Options: "Sample-level", "Set-level", or ""
     - Command-line compatibility - Options: "Yes", "No", and/or "Some optional features incompatible"
     - The version where the last known changes occurred (likely the upcoming version if it is a new workflow -- if the upcoming version number is currently unknown, please use **vX.X.X**)
     - Link to the workflow on Dockstore (if applicable) - Workflow name linked to the information tab on Dockstore.
3. Format this information in a table (see template).
4. Copy the previously gathered information to ==**ALL THREE**== overview tables in `docs/workflows_overview/`:
     - `workflows_alphabetically.md` - Add the workflow in the appropriate spot based on the workflow name.
     - `workflows_kingdom.md` - Add the workflow in the appropriate spot(s) based on the kingdom(s) the workflow is applicable to. Make sure it is added alphabetically within the appropriate subsection(s).
     - `workflows_type.md` - Add the workflow in the appropriate spot based on the workflow type. Make sure it is added alphabetically within the appropriate subsection.
5. Copy the path to the workflow to ==**ALL**== of the appropriate locations in the `mkdocs.yml` file (under the `nav:` section) in the main directory of this repository. These should be the exact same spots as in the overview tables but without additional information. This ensures the workflow can be accessed from the navigation sidebar.
