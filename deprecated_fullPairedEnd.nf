nextflow.enable.dsl=2

// Backward-compatible entrypoint; new runs should use mainRNAseq.nf.
include { Main } from './mainRNAseq.nf'

workflow {
    Main()
}
