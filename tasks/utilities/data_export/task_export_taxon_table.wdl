version 1.0

task export_taxon_table {
  input {
    File? taxon_table
    String? gambit_predicted_taxon
    String? terra_project
    String? terra_workspace
    String? samplename
    Boolean theiaviral_panel = false

    Map[String, String?] columns_to_export

    Int cpu = 1
    Int disk_size = 25
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-06-21"
    Int memory = 2
  }
  File columns_to_export_json = write_json(columns_to_export)
  meta {
    volatile: true
  }
  command <<<
    set -euo pipefail
    touch STATUS

    # capture taxon and corresponding table names from input taxon_table
    taxon_array=($(cut -f1 ~{taxon_table} | tail +2))
    table_array=($(cut -f2 ~{taxon_table} | tail +2))

    # replace whitespace from gambit_predicted_taxon with an underscore
    sample_taxon=$(echo ~{gambit_predicted_taxon} | tr ' ' '_')
  
    # prevent failures
    sample_table=""
    # set taxon and table vars
    echo "Checking if sample taxon should be exported to user-specified taxon table..."
    for index in "${!taxon_array[@]}"; do
      taxon=${taxon_array[$index]}
      table=${table_array[$index]}
      if [[ "${sample_taxon}" == *"${taxon}"* ]]; then
        sample_table=${table}
        break
      else 
        echo "${sample_taxon} does not match ${taxon}."
      fi
    done

    if [ "~{theiaviral_panel}" == "true" ]; then
      # check for taxon "other" in taxon table
      other_species=$(awk -F'\t' '$1=="other" {print $2}' ~{taxon_table})

      if [[ -z "${sample_table}" && -n "${other_species}" ]]; then
        echo "Assigning Other_Species"
        sample_table=${other_species}
      fi
    fi

    if [ -n "${sample_table}" ]; then

      jq -r '[.[] | .left], [.[] | .right] | @tsv' ~{columns_to_export_json} > exported_columns.tsv
     
      UPLOAD_DATE=$(date -I)
    
      echo -e "entity:${sample_table}_id\tupload_date\ttable_created_by\t$(head -n1 exported_columns.tsv)" > terra_table_to_upload.tsv
      awk -v date="$UPLOAD_DATE" -v samplename="~{samplename}" 'NR > 1 {print samplename"\t"date"\texport_taxon_table\t" $0}' exported_columns.tsv >> terra_table_to_upload.tsv

      python3 /scripts/import_large_tsv/import_large_tsv.py --project "~{terra_project}" --workspace "~{terra_workspace}" --tsv terra_table_to_upload.tsv
      echo "~{samplename} was added to the ${sample_table} table" > STATUS
    else
      echo "Table not defined for ~{gambit_predicted_taxon}" > STATUS
    fi
  >>>
  output {
    File? terra_table_to_upload = "terra_table_to_upload.tsv"
    String status = read_string("STATUS")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}