version 1.0

task freyja_microreact {
    input {
        File freyja_parsed_format_tsv
        String freyja_plot_name
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/freyja-microreact:1.0.1"
        Int disk_size = 50
        Int memory = 4
        Int cpu = 2
    }

    command <<<
        # check if all samples failed coverage threshold, if so, create an empty microreact file
        # step before writes "all samples are below coverage" to long format tsv, so check for that string in the tsv file
        if grep -q "all samples are below coverage" ~{freyja_parsed_format_tsv}; then
            echo "All samples are below coverage threshold, creating empty microreact file"
            touch ~{freyja_plot_name}.microreact
        else
            # micro react the long format tsv file
            freyja_microreact.py ~{freyja_parsed_format_tsv} --output ~{freyja_plot_name}.microreact
        fi
    >>>

    output {
        File freyja_microreact_output = "~{freyja_plot_name}.microreact"
    }
    runtime {
        docker: docker
        disk: disk_size
        memory: memory
        cpu: cpu
    }
}