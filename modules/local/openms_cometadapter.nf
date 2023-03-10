process OPENMS_COMETADAPTER {
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
        tuple val(meta), path("*.tsv")  , emit: tsv
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix                 = task.ext.prefix ?: "${mzml.baseName}"
        def fixed_modifications    = params.fixed_mods != " " ? "-fixed_modifications ${params.fixed_mods.tokenize(',').collect { "'${it}'"}.join(" ")}" : "-fixed_modifications"
        def variable_modifications = params.fixed_mods != " " ? "-variable_modifications ${params.fixed_mods.tokenize(',').collect { "'${it}'"}.join(" ")}" : "-fixed_modifications"
        def remove_precursor       = params.remove_precursor_peak ? "-remove_precursor_peak yes" : ""

        """
        rm /usr/local/bin/comet.exe
        CometAdapter -in $mzml \\
            -out ${prefix}.idXML \\
            -database $fasta \\
            -threads $task.cpus \\
            -pin_out ${prefix}.tsv \\
            -comet_executable comet.linux.exe \\
            -precursor_mass_tolerance ${params.precursor_mass_tolerance} \\
            -fragment_mass_tolerance ${params.fragment_mass_tolerance} \\
            -fragment_bin_offset ${params.fragment_bin_offset} \\
            -num_hits ${params.num_hits} \\
            -digest_mass_range ${params.digest_mass_range} \\
            -max_variable_mods_in_peptide ${params.number_mods} \\
            -missed_cleavages 0 \\
            -precursor_charge ${params.prec_charge} \\
            -activation_method ${params.activation_method} \\
            -enzyme '${params.enzyme}' \\
            -spectrum_batch_size ${params.spectrum_batch_size} \\
            -use_X_ions ${params.use_x_ions} \\
            -use_Z_ions ${params.use_z_ions} \\
            -use_A_ions ${params.use_a_ions} \\
            -use_C_ions ${params.use_c_ions} \\
            -use_NL_ions ${params.use_NL_ions} \\
            $remove_precursor \\
            $fixed_modifications \\
            $variable_modifications

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms-thirdparty: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
