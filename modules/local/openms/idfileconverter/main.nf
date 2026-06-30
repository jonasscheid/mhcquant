process OPENMS_IDFILECONVERTER {
    tag "$meta.id"
    label 'process_single'
    label 'openms'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'quay.io/biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(id_file), val(out_type)

    output:
    tuple val(meta), path("${prefix}.${out_type}"), emit: converted
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The only tool that bridges idXML <-> idparquet (regular ID tools preserve format).
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    if ("$id_file" == "${prefix}.${out_type}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    IDFileConverter \\
        -in $id_file \\
        -out ${prefix}.${out_type} \\
        -threads $task.cpus \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${out_type}
    """
}
