version 1.0

task freyja_microreact {
    input {
        File freyja_long_format_tsv
        String freyja_plot_name
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/freyja-microreact:1.0.0"
        Int disk_size = 50
        Int memory = 4
        Int cpu = 2
    }

    command <<<
        # micro react the long format tsv file
        freyja_microreact.py ~{freyja_long_format_tsv} --output ~{freyja_plot_name}.microreact
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