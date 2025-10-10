version 1.0

task lod_table_prep {
  input {
    String input_table_name
    String workspace_name
    String project_name

    String output_table_name

    String read1_column_name = "read1"
    String read2_column_name = "read2"
    String taxon_column_name = "gambit_predicted_taxon"

    Array[Int] downsampling_levels
    Array[String] expected_genes
    Array[String] expected_alleles
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/bioforklift:0.2.5"

    Int memory = 4
    Int cpu = 1
    Int disk_size = 100
  }
  command <<<
    set -euo pipefail

    python3 <<CODE
    from bioforklift.terra import Terra, WorkflowConfig
    import pandas as pd

    # Initialize Terra client
    terra = Terra(
      source_workspace="~{workspace_name}",
      source_project="~{project_name}",
      destination_workspace="~{workspace_name}",
      destination_project="~{project_name}",
    )

    # download the input table as a pandas dataframe
    input_df = terra.entities.download_table("~{input_table_name}")
    input_header = list(input_df.columns)

    required_columns = [
      "~{read1_column_name}",
      "~{read2_column_name}",
      "~{taxon_column_name}",
    ]

    # check that all required columns are present in the input table
    for col in required_columns:
        if col not in input_header:
          raise ValueError(f"ERROR: Column '{col}' not found in input table")

    # write header to output LOD table
    output_df = pd.DataFrame(
      {
        "entity:~{output_table_name}_id" : input_df["entity:~{input_table_name}_id"],
        "taxon": input_df["~{taxon_column_name}"],
        "expected_genes": ",".join(["~{sep=',' expected_genes}"]),
        "expected_alleles": ",".join(["~{sep=',' expected_alleles}"]),
        "read1": input_df["~{read1_column_name}"],
        "read2": input_df["~{read2_column_name}"],
      }
    )

    downsampled_df = pd.DataFrame({"downsampling_level": sorted([~{sep=',' downsampling_levels}])})
    merged_df = output_df.merge(downsampled_df, how="cross")
    merged_df["entity:~{output_table_name}_id"] = merged_df["entity:~{output_table_name}_id"] + "_" + merged_df["downsampling_level"].astype(str) + "x"
    output_df = pd.concat([output_df, merged_df], ignore_index=True)
    output_df.to_csv("~{output_table_name}-data.tsv", sep="\t", index=False)

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