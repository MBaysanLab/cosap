# COSAP
Comparative Sequencing Analysis Platform

This repository hosts the development of the COSAP library.

## About COSAP

COSAP is an easy yet comprehensive pipeline creation tool for NGS data.
It provides reproducibility and aims to give deeper insight about the powers and limitations of the current tools by allowing users to compare results of different pipelines.

A typical variant calling pipeline consists of Read Trimming, Read Mapping, Duplicate Removal&Base Calibration, Variant Calling and Variant Annotation steps.
COSAP provides several tool options for each of these steps.

## Usage

The documentation for this project has moved. Please visit [docs.cosap.bio](https://docs.cosap.bio) for the latest information.

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
