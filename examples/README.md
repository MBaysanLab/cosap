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

