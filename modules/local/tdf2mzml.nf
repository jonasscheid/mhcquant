process TDF2MZML {
    tag "$meta.id"
    label 'process_high'


    if (params.enable_conda && meta.ext == 'd') {
            { exit 1, "Converting bruker tdf file format to mzml is only supported using docker. Aborting." }
    }


    container "mfreitas/tdf2mzml"

    input:
        tuple val(meta), path(brukerfile)

    output:
        tuple val(meta), path("*.mzML"), emit: mzml
        path "versions.yml"            , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${brukerfile.simpleName}"


        """
        tar -xzf $brukerfile -C .
        
        tdf2mzml.py -i ${prefix}.d -o ${prefix}.mzML


        echo 1.0.0 > versions.yml
        """
}
