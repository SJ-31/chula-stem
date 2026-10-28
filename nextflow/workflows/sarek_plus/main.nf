workflow {

    main:
    input = channel.fromPath(workflow.outputDir)
    // TODO: build the manifest file from the outputDir, which is the
    // same output directory of sarek.
    // Throw an error if there aren't enough samples for it to work
}
