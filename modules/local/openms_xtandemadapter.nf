process OPENMS_XTANDEMADAPTER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.8.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:2.8.0--h9ee0642_2' :
        'quay.io/biocontainers/openms-thirdparty:2.8.0--h9ee0642_2' }"

    input:
        tuple val(meta), path(mzml), path(fasta)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix                 = task.ext.prefix ?: "${mzml.baseName}"
        def fixed_modifications    = params.fixed_mods != " " ? "-fixed_modifications ${params.fixed_mods.tokenize(',').collect { "'${it}'"}.join(" ")}" : "-fixed_modifications"
        def variable_modifications = params.variable_mods != " " ? "-variable_modifications ${params.variable_mods.tokenize(',').collect { "'${it}'"}.join(" ")}" : "-variable_modifications"

        """
        XTandemAdapter -in $mzml \\
            -out ${prefix}.idXML \\
            -database $fasta \\
            -threads $task.cpus \\
            -precursor_mass_tolerance ${params.precursor_mass_tolerance} \\
            -fragment_mass_tolerance ${params.fragment_mass_tolerance} \\
            -missed_cleavages 0 \\
            -precursor_charge ${params.prec_charge} \\
            -activation_method ${params.activation_method} \\
            -enzyme '${params.enzyme}' \\
            -max_precursor_charge 3 \\
            -variable_modifications 'Oxidation (M)'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms-thirdparty: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
