nextflow.enable.dsl=2

include { DownloadReference; buildref; ValidateHotspots; PublishReference } from './modules/index.nf'

params.reference_name = params.reference_name ?: "ensembl_v${params.ensembl_release}_GRCh38"

workflow {
    if (!(params.ensembl_release.toString() ==~ /\d+/)) error 'Set --ensembl_release to a numeric Ensembl release'
    if (params.vep_cache_version && params.vep_cache_version.toString() != params.ensembl_release.toString())
        error 'VEP cache must use --ensembl_release; remove the separate --vep_cache_version override'
    if (!(params.reference_name ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Invalid --reference_name'
    reference_root = new File(params.reference_build_root.toString()).canonicalFile
    if (!new File(params.reference_build_root.toString()).isAbsolute() || !reference_root.isDirectory() || !reference_root.canWrite())
        error 'Set --reference_build_root to an existing, writable absolute directory'
    if (new File(reference_root, params.reference_name).exists())
        error "Reference destination already exists: ${reference_root}/${params.reference_name}; choose a new name"

    release = params.ensembl_release
    sources = [
        fasta: params.fasta_ref ?: "https://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        gtf: params.GTF ?: "https://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release}.chr.gtf.gz",
        cdna: params.cdna ?: "https://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        vep: params.vep_cache_url ?: "https://ftp.ensembl.org/pub/release-${release}/variation/indexed_vep_cache/homo_sapiens_vep_${release}_GRCh38.tar.gz",
        somalier: params.somalier_sites_url
    ]
    DownloadReference(sources)
    buildref(DownloadReference.out.files, file(params.hotspot_catalogue, checkIfExists: true),
        params.transcript_overrides ? file(params.transcript_overrides, checkIfExists: true) : [],
        file("${projectDir}/reference/reference_metadata.py"), file("${projectDir}/reference/hotspot_reference.py"))
    ValidateHotspots(buildref.out.reference, file("${projectDir}/reference/hotspot_reference.py"))
    PublishReference(buildref.out.reference, ValidateHotspots.out.validation,
        file("${projectDir}/reference/reference_metadata.py"))
}
