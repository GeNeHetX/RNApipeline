nextflow.enable.dsl=2

process STAGE_WORKFLOW_RESOURCES {
    tag 'RNApipeline workflow resources'
    executor 'local'
    scratch false
    stageInMode 'copy'

    input:
    path resource_source, name: 'resources.source'

    output:
    path 'resources', emit: resources

    script:
    """
    mkdir resources
    cp -a resources.source/. resources/
    """
}
