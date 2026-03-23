version 1.0

task lod_table_prep {
  input {
    File lod_config_yaml
    String workspace_name
    String project_name
    String input_table_name
    String output_table_name

    String read1_column_name
    String read2_column_name
    String taxon_column_name
    Array[Int] downsampling_levels

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/bioforklift:0.3.9-dev"

    Int memory = 4
    Int cpu = 1
    Int disk_size = 100
  }
  command <<<
    set -euo pipefail

    python3 <<CODE
    from bioforklift.terra import Terra
    import pandas as pd
    import yaml
    import json

    # load config yaml
    with open("~{lod_config_yaml}", 'r') as f:
      config = yaml.safe_load(f)
    print("Loaded LOD config YAML file.")

    # set parameters
    workspace_name = "~{workspace_name}"
    project_name = "~{project_name}"
    input_table_name = "~{input_table_name}"
    output_table_name = "~{output_table_name}"
    read1_column_name = "~{read1_column_name}"
    read2_column_name = "~{read2_column_name}"
    taxon_column_name = "~{taxon_column_name}"
    downsampling_levels = sorted([int(x) for x in "~{sep=',' downsampling_levels}".split(",")])
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
        f"{read1_column_name}": input_df[f"{read1_column_name}"],
        f"{read2_column_name}": input_df[f"{read2_column_name}"],
        f"{taxon_column_name}": input_df[f"{taxon_column_name}"],
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
      output_df[f"reportable_{target['target_name']}"] = [','.join(sorted(target['target_values']))] * num_rows

      # creating columns of 'expected_{target_name}' values
      # First initialize the column with None since not all samples will have the same expected values
      output_df[f"expected_{target['target_name']}"] = None
      for sample, value in target['expected_values'].items():
        # value is only assigned if sample from expected_values matches the sample_id in the input table
        mask = output_df[f"entity:{output_table_name}_id"] == sample
        if not mask.any():
          raise ValueError(f"No sample_id matches '{sample}' in expected_values for target '{target['target_name']}'")
        row_index = output_df.index[mask][0]
        output_df.at[row_index, f"expected_{target['target_name']}"] = ','.join(sorted(value))
    print("Populated all reportable and expected values for all targets.")

    # create downsampled sample for each level across all samples
    downsampled_df = pd.DataFrame({"downsampling_level": downsampling_levels})
    merged_df = output_df.merge(downsampled_df, how="cross")

    # modify entity ids (sample_id names) to reflect downsampling level in output table
    merged_df[f"entity:{output_table_name}_id"] = merged_df[f"entity:{output_table_name}_id"] + "_" + merged_df["downsampling_level"].astype(str) + "x"

    # combine original samples with all downsampled samples
    output_df = pd.concat([output_df, merged_df], ignore_index=True)
    output_df.to_csv(f"{output_table_name}_datatable.tsv", sep="\t", index=False)
    print(f"Generated LOD table with {len(output_df)} total entries including downsampled levels.")

    # upload the output table to Terra
    updated_df = terra.entities.upload_entities(data=output_df, target=f"{output_table_name}")
    result = terra.entities.create_entity_set(f"{output_table_name}_set", f"{output_table_name}", updated_df)
    if result.ok:
      print("Entity set created successfully")

    # write input and output json for launching the wf_lod_table_process workflow
    input_json = {
      f"lod_table_process.taxon": f"this.{taxon_column_name}",
      f"lod_table_process.read1": f"this.{read1_column_name}",
      f"lod_table_process.read2": f"this.{read2_column_name}",
      f"lod_table_process.samplename": f"this.{output_table_name}_id",
      f"lod_table_process.downsampling_level": "this.downsampling_level",
    }
    with open("wf_lod_table_process.input.json", 'w') as f:
      json.dump(input_json, f, indent=2)

    output_json = {
      f"lod_table_process.ncbi_taxon_id": "this.ncbi_taxon_id",
      f"lod_table_process.ncbi_taxon_name": "this.ncbi_taxon_name",
      f"lod_table_process.ncbi_read_extraction_rank": "this.ncbi_read_extraction_rank",
      f"lod_table_process.ete4_status": "this.ete4_status",
      f"lod_table_process.taxon_avg_genome_length": "this.taxon_avg_genome_length",
      f"lod_table_process.cg_pipeline_est_coverage": "this.cg_pipeline_est_coverage",
      f"lod_table_process.cg_pipeline_report": "this.cg_pipeline_report",
      f"lod_table_process.rasusa_status": "this.rasusa_status",
      f"lod_table_process.read1_subsampled": "this.read1_subsampled",
      f"lod_table_process.read2_subsampled": "this.read2_subsampled",
    }
    with open("wf_lod_table_process.output.json", 'w') as f:
      json.dump(output_json, f, indent=2)

    CODE

    # SUBMISSION TO BIOFORKLIFT
    DATE=$(date +"%Y-%m-%d %T")

    # wf_lod_table_process SUBMISSION
    bioforklift launch \
      --workspace "~{workspace_name}" \
      --project "~{project_name}" \
      --workflow_name "LOD_Table_Process_PHB" \
      --branch "tj-lod-lift" \
      --table "~{output_table_name}" \
      --comment "${DATE} Automated submission of LOD_Table_Process_PHB" \
      --input_json "wf_lod_table_process.input.json" \
      --output_json "wf_lod_table_process.output.json" \
      --verbose

  >>>
  output {
    File lod_table_tsv = "~{output_table_name}_datatable.tsv"
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