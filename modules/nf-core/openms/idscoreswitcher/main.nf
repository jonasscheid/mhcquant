process OPENMS_IDSCORESWITCHER {
    tag "$meta.id"
    label 'process_single'
    label 'openms'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0':
        'quay.io/biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(idxml)

    output:
    tuple val(meta), path("*.idparquet"), emit: idxml
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$idxml" == "${prefix}.idparquet") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    IDScoreSwitcher \\
        -in $idxml \\
        -out ${prefix}.idparquet \\
        -threads $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$idxml" == "${prefix}.idparquet") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    mkdir ${prefix}.idparquet
    """
}
