# COSAP
Comparative Sequencing Analysis Platform

This repository hosts the development of the COSAP library.

## About COSAP

COSAP is an easy yet comprehensive pipeline creation tool for NGS data.
It provides reproducibility and aims to give deeper insight about the powers and limitations of the current tools by allowing users to compare results of different pipelines.

COSAP owns pipeline and bioinformatics-tool orchestration. Generic format parsing is
provided by [`cosap-formats`](https://github.com/marifergun/cosap_formats), while
pipeline-specific result parsing and GATK, tabix, comparator, and annotation
processors remain in this package. Historical parser imports continue to work as
compatibility APIs.

A typical variant calling pipeline consists of Read Trimming, Read Mapping, Duplicate Removal&Base Calibration, Variant Calling and Variant Annotation steps.
COSAP provides several tool options for each of these steps.

For local installation and container build instructions that use a sibling
`cosap_formats` checkout, see [Source builds](docs/source-builds.md).

## Usage

The documentation for this project has moved. Please visit [docs.cosap.bio](https://docs.cosap.bio) for the latest information.

### Mitochondrial Variant Calling

COSAP supports GATK's hg38 mitochondrial SNV/indel workflow with dual standard and shifted `chrM` references, mitochondrial-mode Mutect2, blacklist filtering, and optional Haplochecker contamination reporting. It accepts either an analysis-ready indexed BAM or paired FASTQs; FASTQ input is trimmed, mapped against the full hg38 reference, duplicate-marked, and recalibrated before mitochondrial calling.

The hg38 resource bundle must provide the GATK `chrM` standard and shifted references, BWA sidecars, dictionaries, shift-back chain, and blacklist under `chrM/`. The workflow currently runs on CPU and accepts BAM, not CRAM, input.

```bash
cosap mitochondrial \
  --workdir /workdir/mitochondrial \
  --sample-name sample1 \
  --bam /data/sample1.bam \
  --bam-index /data/sample1.bam.bai \
  --median-autosomal-coverage 30
```

### Apptainer GPU Usage

COSAP can launch child Apptainer containers for GPU-backed tools such as Parabricks.

- Launch the parent COSAP container with `apptainer exec --nv ...` so GPU devices and NVIDIA libraries are available inside COSAP.
- In `auto` mode, COSAP now detects that it is already running inside Apptainer and selects the Apptainer backend.
- GPU child containers are launched with `apptainer exec --nv ...` when `--device gpu` is used.
- COSAP bundles an unprivileged Apptainer install in the Docker image by default.
- `COSAP_APPTAINER_BIN` remains available as an override if you want to point COSAP at a host-provided Apptainer installation instead.

Example:

```bash
apptainer exec --nv \
  --bind hg38_bundle/:/cosap_data \
  --bind cosap_test_data/:/workdir \
  docker://itubioinformatics/cosap:<version> \
  cosap pipeline \
    --analysis_type germline \
    --workdir /workdir/test_gpu \
    --normal_sample /workdir/SRR31264074_1.fastq.gz,/workdir/SRR31264074_2.fastq.gz \
    --mapper bwa \
    --variant_caller deepvariant \
    --device gpu \
    --model_type WES
```

If you want to override the bundled Apptainer with a host installation, mount a full host Apptainer install and set `COSAP_APPTAINER_BIN` to the mounted executable. Mounting only `/usr/bin/apptainer` is not sufficient because Apptainer also needs its config and support files such as `apptainer.conf`.

## Licence

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Opening an issue

Please post bug reports
in [GitHub issues](https://github.com/MBaysanLab/cosap/issues).


---
