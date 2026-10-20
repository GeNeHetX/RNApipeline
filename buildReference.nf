nextflow.enable.dsl=2

include { DownloadReference; buildref } from './modules/index.nf'

params.reference_build_root = params.reference_build_root ?: '/ref'
params.reference_name = params.reference_name ?: "ensembl_v${params.ensembl_release}_GRCh38"

workflow {
    if (!(params.ensembl_release.toString() ==~ /\d+/))
        error 'Set --ensembl_release to a numeric Ensembl release'
    if (!(params.vep_cache_version.toString() ==~ /\d+/))
        error 'Set --vep_cache_version to a numeric VEP cache release'
    if (!(params.reference_name ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
        error 'Invalid --reference_name'

    reference_root = new File(params.reference_build_root.toString())
    if (!reference_root.isAbsolute())
        error "Reference output root must be an absolute path: ${reference_root}"
    reference_root = reference_root.canonicalFile
    if (!reference_root.isDirectory())
        error "Reference output root must be an existing directory: ${reference_root}"
    if (!reference_root.canWrite())
        error "No write permission for reference output root: ${reference_root}"

    reference_destination = new File(reference_root, params.reference_name)
    if (reference_destination.exists())
        error "Reference destination already exists: ${reference_destination}. Use it with mainRNAseq or choose a new name."

    release = params.ensembl_release
    vep_release = params.vep_cache_version
    fasta_url = params.fasta_ref ?: "https://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    gtf_url = params.GTF ?: "https://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release}.chr.gtf.gz"
    cdna_url = params.cdna ?: "https://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    known_vcf_url = params.known_vcf ?: "https://ftp.ensembl.org/pub/release-${release}/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz"
    vep_cache_url = params.vep_cache_url ?: "https://ftp.ensembl.org/pub/release-${vep_release}/variation/indexed_vep_cache/homo_sapiens_vep_${vep_release}_GRCh38.tar.gz"

    DownloadReference(
        fasta_url,
        gtf_url,
        cdna_url,
        known_vcf_url,
        vep_cache_url,
        vep_release,
        true
    )
    buildref(DownloadReference.out.files)
}
