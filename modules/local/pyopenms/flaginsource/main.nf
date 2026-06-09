process PYOPENMS_FLAGINSOURCE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyopenms:3.4.1--py312h6b06db6_2' :
        'biocontainers/pyopenms:3.4.1--py312h6b06db6_2' }"

    input:
    tuple val(meta), path(fdr_filtered), path(rescored)

    output:
    tuple val(meta), path("*_insource_filtered.idXML"), emit: idxml
    tuple val(meta), path("*_insource_rescored.idXML"), emit: rescored_idxml
    tuple val("${task.process}"), val('pyopenms'), eval("pip show pyopenms | grep Version | sed 's/Version: //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    flag_in_source_fragments.py \\
        --filtered $fdr_filtered \\
        --rescored $rescored \\
        --out ${prefix}_insource_filtered.idXML \\
        --out-rescored ${prefix}_insource_rescored.idXML \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_insource_filtered.idXML
    touch ${prefix}_insource_rescored.idXML
    """
}
