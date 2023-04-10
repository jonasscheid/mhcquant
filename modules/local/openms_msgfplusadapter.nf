process OPENMS_MSGFPLUSADAPTER {
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
        def prefix           = task.ext.prefix ?: "${mzml.baseName}"

        """
        MSGFPlusAdapter -in $mzml \\
            -out ${prefix}.idXML \\
            -database $fasta \\
            -threads $task.cpus \\
            -precursor_mass_tolerance ${params.precursor_mass_tolerance} \\
            -fragment_method 'CID' \\
            -instrument 'high-res' \\
            -enzyme unspecific cleavage \\
            -protocol none \\
            -min_precursor_charge 2 \\
            -max_precursor_charge 3 \\
            -min_peptide_length 8 \\
            -max_peptide_length 12 \\
            -matches_per_spec ${params.num_hits} \\
            -max_mods ${params.number_mods} \\
            -variable_modifications 'Oxidation (M)'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms-thirdparty: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
