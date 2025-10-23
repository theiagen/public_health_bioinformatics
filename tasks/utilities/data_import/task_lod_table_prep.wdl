version 1.0

task lod_table_prep {
  input {
    File lod_config_yaml
    String? input_table_name
    String? workspace_name
    String? project_name
    String? output_table_name
    String read1_column_name = "read1"
    String read2_column_name = "read2"
    Array[Int]? downsampling_levels

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/bioforklift:0.2.5-dev"

    Int memory = 4
    Int cpu = 1
    Int disk_size = 100
  }
  command <<<
    set -euo pipefail

    python3 <<CODE
    from bioforklift.terra import Terra, WorkflowConfig
    import pandas as pd
    import yaml

    # load config yaml
    with open("~{lod_config_yaml}", 'r') as f:
      config = yaml.safe_load(f)
    print("Loaded LOD config YAML file.")

    required_params = {
      'workspace_name': "~{workspace_name}",
      'project_name': "~{project_name}",
      'input_table_name': "~{input_table_name}",
      'output_table_name': "~{output_table_name}",
      'read1_column_name': "~{read1_column_name}",
      'read2_column_name': "~{read2_column_name}",
      'downsampling_levels': "~{sep=',' downsampling_levels}"
    }

    # check for required input parameters
    for param in required_params:
      if not config['workflow'].get(param) and not required_params[param]:
        raise ValueError(f"ERROR: Missing required input parameter: '{param}'. Not found in config or optional WDL inputs")
    print("All required parameters found.")

    # set parameters, giving precedence to WDL inputs if provided
    workspace_name = "~{workspace_name}" if "~{workspace_name}" else config['workflow']['workspace_name']
    project_name = "~{project_name}" if "~{project_name}" else config['workflow']['project_name']
    input_table_name = "~{input_table_name}" if "~{input_table_name}" else config['workflow']['input_table_name']
    output_table_name = "~{output_table_name}" if "~{output_table_name}" else config['workflow']['output_table_name']
    read1_column_name = "~{read1_column_name}" if "~{read1_column_name}" else config['workflow']['read1_column_name']
    read2_column_name = "~{read2_column_name}" if "~{read2_column_name}" else config['workflow']['read2_column_name']
    downsampling_levels = sorted([int(x) for x in "~{sep=',' downsampling_levels}"]) if "~{sep=',' downsampling_levels}" else sorted([int(x) for x in config['workflow']['downsampling_levels']])

    # Initialize Terra client
    terra = Terra(
      source_workspace = workspace_name,
      source_project = project_name,
      destination_workspace = workspace_name,
      destination_project = project_name,
    )

    # download the input table as a pandas dataframe
    input_df = terra.entities.download_table(input_table_name)
    input_header = list(input_df.columns)

    # fill in initial/expected values into the output LOD table
    output_df = pd.DataFrame(
      {
        f"entity:{output_table_name}_id" : input_df[f"entity:{input_table_name}_id"],
        "read1": input_df[f"{read1_column_name}"],
        "read2": input_df[f"{read2_column_name}"],
      }
    )
    num_rows = len(output_df)
    print(f"Preparing LOD table for {num_rows} samples...")

    # load reportables, expected values, and make sure columns to compare exist
    for target in config['reportables']:
      # copy the column to compare from input table to output table if it exists
      if target.get('column_to_compare') not in input_header:
        raise ValueError(f"ERROR: Column '{target['column_to_compare']}' not found in input table")
      else:
        output_df[target['column_to_compare']] = input_df[target['column_to_compare']]

      # creating columns of 'reportable_{target_name}' values
      output_df[f"reportable_{target['target_name']}"] = [target['target_values']] * num_rows

      # creating columns of 'expected_{target_name}' values
      for sample, value in target['expected_values'].items():
        # value is only assigned if sample from expected_values matches the sample_id in the input table
        mask = output_df[f"entity:{output_table_name}_id"] == sample
        if not mask.any():
          raise ValueError(f"No sample_id matches '{sample}' in expected_values for target '{target['target_name']}'")
        output_df.loc[mask, f"expected_{target['target_name']}"] = [value] * num_rows
    print("Populated all reportable and expected values for all targets.")

    # create downsampled sample for each level across all samples
    downsampled_df = pd.DataFrame({"downsampling_level": downsampling_levels})
    merged_df = output_df.merge(downsampled_df, how="cross")

    # modify entity ids (sample_id names) to reflect downsampling level in output table
    merged_df["entity:~{output_table_name}_id"] = merged_df["entity:~{output_table_name}_id"] + "_" + merged_df["downsampling_level"].astype(str) + "x"

    # combine original samples with all downsampled samples
    output_df = pd.concat([output_df, merged_df], ignore_index=True)
    output_df.to_csv("~{output_table_name}-data.tsv", sep="\t", index=False)
    print(f"Generated LOD table with {len(output_df)} total entries including downsampled levels.")

    # upload the output table to Terra
    updated_df = terra.entities.upload_entities(data=output_df, target="~{output_table_name}")
    result = terra.entities.create_entity_set("~{output_table_name}_set", "~{output_table_name}", updated_df)
    if result.ok:
        print("Entity set created successfully")

    CODE

  >>>
  output {
    File lod_table_tsv = "~{output_table_name}-data.tsv"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}