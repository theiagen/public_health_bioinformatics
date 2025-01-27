version 1.0

task parse_amr_json {
  input {
    File input_json  # Input JSON file
    String output_csv_name = "amr_table.csv"  # Name for the output CSV
    String output_png_name = "amr_table.png"  # Name for the output PNG
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/amrsearch:0.1.0"
    Int cpu = 2
    Int memory = 4
    Int disk_size = 20
  }

  command <<< 
    # Create a Python script to parse the JSON and generate a PNG
    cat << EOF > parse_amr_json.py
import json
import csv
import matplotlib.pyplot as plt

def parse_amr_json(input_json, output_csv, output_png):
    with open(input_json, 'r') as f:
        data = json.load(f)

    library_version = data['library']['version']
    library_label = data['library']['label']
    resistance_profiles = data['resistanceProfile']

    table_data = []
    for profile in resistance_profiles:
        agent = profile['agent']['name']
        inferred_resistance = profile['state'].lower()
        determinants = []

        acquired = profile['determinants'].get('acquired', [])
        for item in acquired:
            determinants.append(item.get('gene', ''))

        variants = profile['determinants'].get('variants', [])
        for item in variants:
            gene = item.get('gene', '')
            variant = item.get('variant', '')
            determinants.append(f"{gene}_{variant}")

        determinants_str = '; '.join(determinants) if determinants else 'none'
        table_data.append([agent, inferred_resistance, determinants_str])

    # Write to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Agent", "Inferred Resistance", "Known Determinants"])
        csvwriter.writerows(table_data)

    # Generate PNG
    agents = [row[0] for row in table_data]
    resistances = [row[1] for row in table_data]
    determinants = [row[2] for row in table_data]

    plt.figure(figsize=(10, len(agents) * 0.5))
    plt.table(cellText=table_data, colLabels=["Agent", "Inferred Resistance", "Known Determinants"], loc='center')
    plt.axis('off')
    plt.title(f"AMR - Antimicrobial resistance (Library {library_label}, Version {library_version})")
    plt.savefig(output_png, bbox_inches='tight')

parse_amr_json("~{input_json}", "~{output_csv_name}", "~{output_png_name}")
EOF

    # Run the Python script
    python3 parse_amr_json.py
  >>>

  output {
    File output_csv = "~{output_csv_name}" 
    File output_png = "~{output_png_name}" 
  }

  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 1
  }
}
