process PREDICT_POSSIBLE_CLASS1_NEOEPITOPES {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0' :
        'quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0' }"

    input:
        tuple val(meta), val(alleles), path(vcf)

    output:
        tuple val(meta), path("*.csv"), emit: csv
        tuple val(meta), path("*.txt"), emit: txt
        path "versions.yml"           , emit: versions

    script:
        def prefix           = task.ext.prefix ?: "${meta}_vcf_neoepitopes_class1"

        """
        vcf_neoepitope_predictor.py \\
            -t ${params.variant_annotation_style} \\
            -r ${params.variant_reference} \\
            -a '$alleles' -minl ${params.peptide_min_length} \\
            -maxl ${params.peptide_max_length} \\
            -v $vcf \\
            -o ${prefix}.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mhcflurry: \$(echo \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//') )
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
            fred2: \$(echo \$(python -c "import pkg_resources; print('fred2' + pkg_resources.get_distribution('Fred2').version)" | sed 's/^fred2//; s/ .*\$//'))
        END_VERSIONS
        """
}
