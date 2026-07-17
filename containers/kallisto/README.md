# Kallisto SIFs

Build the two architecture-specific Kallisto 0.51.1 SIFs before building the
reference. The destination argument is a directory; each script writes
`kallisto-0.51.1.sif` there.

Run the ARM64 script on a PAM worker. Run the x86_64 script on a TOD worker
(Intel and AMD are both x86_64). The scripts require Apptainer and write the
SIF atomically, so a failed build does not replace an existing file.

Example destination layout:

```text
/ref/tools/rnapipeline/v1.7.0/kallisto/0.51.1/arm64/
/ref/tools/rnapipeline/v1.7.0/kallisto/0.51.1/amd64/
```

Direct execution on a matching worker:

```bash
containers/kallisto/build-arm64-sif.sh \
  /ref/tools/rnapipeline/v1.7.0/kallisto/0.51.1/arm64

containers/kallisto/build-x86_64-sif.sh \
  /ref/tools/rnapipeline/v1.7.0/kallisto/0.51.1/amd64
```

From a controller, submit the same scripts to the matching Slurm partition.
Use `pam_cpu` for ARM64. Use the enabled x86/TOD burst partition for x86_64;
the burst partition must be enabled before submission.

After both files exist, run `PrePostScripts/ref_build.nf`. It validates and
uses the matching SIF but never creates or publishes SIFs.
