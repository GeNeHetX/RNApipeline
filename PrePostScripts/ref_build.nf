nextflow.enable.dsl=2

params.ref_root = '/ref'
params.work_root = '/srv/slurm/scratch'
params.force = false

process BUILD_REFERENCE {
    tag 'RNApipeline reference'
    cpus 16
    memory '64 GB'
    time '48h'

    input:
    path builder, stageAs: 'ref_build_tod.sh'

    output:
    path 'reference-build.done'

    script:
    def force_arg = params.force ? '--force' : ''
    """
    bash ref_build_tod.sh \
        --ref-root '${params.ref_root}' \
        --work-root '${params.work_root}' \
        ${force_arg}
    touch reference-build.done
    """
}

workflow {
    builder = Channel.fromPath(
        "${projectDir}/ref_build_tod.sh",
        checkIfExists: true
    )
    BUILD_REFERENCE(builder)
}
