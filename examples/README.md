# This folder contains pipeline examples created with COSAP

# Data

Fastq files used in example pipeline are from [Seqc2](https://www.nature.com/articles/s41587-021-00993-6) study. You can download them from following links:

Tumor: https://www.ebi.ac.uk/ena/browser/view/SRR7890850

Normal: https://www.ebi.ac.uk/ena/browser/view/SRR7890851

Once you have downloaded the data, place it in a directory named `data/` in the same folder as `example_run.py`.

Activate the conda environment
```bash
conda activate cosap
```
Run the example pipeline and that is it. 
```bash
python example_run.py
```

# Apptainer GPU note

When running COSAP itself inside Apptainer and using `--device gpu`, launch the parent container with `--nv`.

```bash
apptainer exec --nv \
  --bind /path/to/reference:/cosap_data \
  --bind /path/to/workdir:/workdir \
  docker://itubioinformatics/cosap:<version> \
  cosap pipeline ... --device gpu
```

COSAP will select the Apptainer backend automatically when already inside an Apptainer container, and GPU child containers such as Parabricks are launched with `--nv` as well. The COSAP Docker image bundles Apptainer by default. If you override it with a host install via `COSAP_APPTAINER_BIN`, mount a full Apptainer install, not just `/usr/bin/apptainer`, because the child runtime also needs its config and helper files.
