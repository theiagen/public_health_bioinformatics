version 1.0

task parse_amr_json {
  input {
    File input_json  # Input JSON file
    String samplename
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

def parse_amr_json(input_json, samplename):
    with open(input_json, 'r') as f:
        data = json.load(f)

    library_version = data['library']['version']
    library_label = data['library']['label']
    resistance_profiles = data['resistanceProfile']

    # Define output filenames using samplename
    output_csv = f"{samplename}_amr_results.csv"
    output_png = f"{samplename}_amr_results.png"
    output_version_txt = "output_amr_version.txt"

    # Save the AMR search version to output
    with open(output_version_txt, 'w') as f:
        f.write(library_version + "\n")

    table_data = []
    for profile in resistance_profiles:
        agent = profile['agent']['name']
        inferred_resistance = profile['state'].capitalize()  # Capitalizing for consistency
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

    # Generate PNG with title
    fig, ax = plt.subplots(figsize=(10, len(table_data) * 0.5))
    ax.set_title(f"AMR - Antimicrobial Resistance (Library {library_label}, Version {library_version})", fontsize=12, weight='bold')
    ax.axis('off')

    # Create table
    table = plt.table(cellText=table_data, colLabels=["Agent", "Inferred Resistance", "Known Determinants"],
                      cellLoc='center', loc='center', bbox=[0, 0, 1, 1])

    # Style headers
    for col, label in enumerate(["Agent", "Inferred Resistance", "Known Determinants"]):
        cell = table[0, col]
        cell.set_text_props(weight="bold")
        cell.set_facecolor("lightgrey")

    # Style rows and highlight "Resistant" rows
    for row_idx, row in enumerate(table_data, start=1):
        for col_idx, value in enumerate(row):
            cell = table[row_idx, col_idx]
            if row[1] == "Resistant":  # Highlight row if "Resistant"
                cell.set_text_props(weight="bold")

    plt.savefig(output_png, bbox_inches='tight', dpi=300)

# Run the function
parse_amr_json("~{input_json}", "~{samplename}")
EOF

    # Run the Python script
    python3 parse_amr_json.py
  >>>

  output {
    File output_csv = "~{samplename}_amr_results.csv"
    File output_png = "~{samplename}_amr_results.png"
    File output_version = "output_amr_version.txt" 
  }

  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 3
  }
}
