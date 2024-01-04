process OPENMS_DECOYDATABASE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::openms=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.1.0--h8964181_3' :
        'biocontainers/openms:3.1.0--h8964181_3' }"

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("*.fasta"), emit: decoy
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${fasta.baseName}_decoy"

        """
        DecoyDatabase -in $fasta \\
            -out ${prefix}.fasta \\
            -decoy_string DECOY_ \\
            -decoy_string_position prefix \\
            -enzyme 'no cleavage'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
