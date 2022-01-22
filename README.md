# CoSAP
Comparative Sequencing Analysis Platform

This repository hosts the development of the CoSAP library.

## About CoSAP

CoSAP is an easy yet comprehensive pipeline creation tool for NGS data. It provides reproducibility and aims to give deeper insight about the powers and limitations of the current tools by
allowing users to compare results of different pipelines.

A typical variant calling pipeline consists of Read Trimming, Read Mapping, Duplicate Removal&Base Calibration, Variant Calling and Variant Annotation steps.
CoSAP provides several tool options for each of these steps.


## Prerequisites and Installation

CoSAP requires [Miniconda](https://docs.conda.io/en/latest/miniconda.html) to install required packages.
Once you install the miniconda,
clone the CoSAP repository:
```bash
git@github.com:MBaysanLab/cosap.git
```
cd into CoSAP:
```bash
cd cosap
```
and run: 
```bash
make install
```
---

## Usage

CoSAP consists of two main modules which are config creators and pipeline runner.

### Config Creation for Pipeline Steps

#### Fastq Reader
File handler for fastq files. This is generally where the pipeline starts.
```python 
from cosap import FastqReader

sample_fastq = FastqReader("/path/to/fastq.fastq", name="normal_sample")
```

Several fastq handler for a sample can be created as:
```python
germline_fastqs = [
    FastqReader("/path/to/fastq_1.fastq", name="normal_sample", read=1)
    FastqReader("/path/to/fastq_2.fastq", name="normal_sample", read=2)
]
```
The read argument is to mark the pairs.

#### Trimmer
Trimmer builder for adaptor trimming and quality control.
Takes the fastq file handler as input.

```python
from cosap import Trimmer

trimmer_germline = Trimmer(reads=germline_fastqs)
```
Here the germline_fastqs is the list of fastq file handler.

#### Mapper
Mapper builder for read mapping. Takes the Trimmer or FastqReader as input.
Currently following libraries are supported:
- [BWA](https://github.com/lh3/bwa)
- [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
- [Bowtie2](https://github.com/BenLangmead/bowtie2)



```python
from cosap import Mapper

germline_params = {
    "read_groups": {
    "ID": "H0164.2",
    "SM": "Pt28N",
    "PU": "0",
    "PL": "illumina",
    "LB": "Solexa-272222"
}

mapper_germline_bwa = Mapper(
    library="bwa",
    reads=trimmer_germline,
    params=germline_params
    }
)

mapper_germline_bowtie = Mapper(
    library="bowtie",
    reads=trimmer_germline,
    params=germline_params
    }
)

```

#### Duplicate Remover
Duplicate read tagger and remover builder. Takes Mapper as input.
```python
from cosap import MDUP

mdup_germline = MDUP(input_step=mapper_germline_bwa)
```

#### Base Recalibrator
GATK BaseRecalibrator builder. Takes Mapper or MDUP as input.
```python 
from cosap import Recalibrator

recalibrator_germline = Recalibrator(input_step=mdup_germline)
```

#### Elprep Preprocessing Tool
[Elprep](https://github.com/ExaScience/elprep) is a high performance tool for preprocessing.
Its functionality is same as duplicate remover&base recalibrator combined.
This tool requires up to 200GB of memory therefore is only recommended to be used on capable workstations and servers.

```python
from cosap import Elprep

elprep_recalibrator_germline = Elprep(input_step=mapper_germline_bwa)
```

#### Variant Caller
Variant caller builder for variant detection tools. Takes Mapper, MDUP, Recalibrator or Elprep of both normal and tumor samples as input.
Currently following libraries are supported:
 - [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360046788432-Mutect2)
 - [Varscan2](http://varscan.sourceforge.net/)
 - [Strelka2](https://github.com/Illumina/strelka)
 - [Octopus](https://github.com/luntergroup/octopus)
 - [MuSe](https://github.com/danielfan/MuSE)
 - [VarDict](https://github.com/AstraZeneca-NGS/VarDict)
 - [SomaticSniper](https://github.com/genome/somatic-sniper)

```python
from cosap import VariantCaller

sample_params = {"germline_sample_name":"Pt28N"}

mutect_caller = VariantCaller(
    library="mutect", 
    germline=recalibrator_germline, 
    tumor=recalibrator_tumor, 
    params=sample_params
)
strelka_caller = VariantCaller(
    library="strelka", 
    germline=recalibrator_germline, 
    tumor=recalibrator_tumor, 
    params=sample_params
)
```
If sample name is provied in the Mapper as read group, it must be provided in the VariantCaller params as well.

#### Variant Annotation
Variand annotator builder. Currently only supports [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)
```python
from cosap import Annotator

annotator = Annotator(library="vep", input_step=mutect_caller)
```

### Pipeline Creation
Here is the `Pipeline`:
```python
from cosap import Pipeline

pipeline = Pipeline()
```

Stacking steps into pipeline is easy as `.add()`:
```python
pipeline.add(trimmer_germline)
pipeline.add(trimmer_tumor)
pipeline.add(mapper_germline_bwa)
pipeline.add(mapper_tumor_bwa)
pipeline.add(mdup_germline)
pipeline.add(mdup_tumor)
pipeline.add(recalibrator_germline)
pipeline.add(recalibrator_tumor)
pipeline.add(mutect_caller)
pipeline.add(annotator)
```
> :warning: **Warning**: You need to add every pipeline step you have created!

`Config` can be created with `.build()`:
```python
pipeline_config = pipeline.build(workdir="/path/to/pipeline/workdir")
```
The config will be saved in the workdir when pipeline is run.

### Pipeline Runner
Runner module takes the pipeline config from previous states and run the pipeline using several steps. Currently it only supports [Snakemake](https://snakemake.readthedocs.io/en/stable/) backend.

```python
from cosap import PipelineRunner

runner = PipelineRunner()

runner.run_pipeline(pipeline_config=pipeline_config,backend="snakemake")
```
---
## Road Map
 - Extend tools for Trimmer, Mapper and preprocessing steps.
 - Support tumor only variant calling.
 - Add comprehensive quality control step.
 - Extend analysis steps with CNV and Fusion Detection.
 - Extend annotator with Annovar and Ella Anno.
 - Support more tools as backend.


## Opening an issue

Please post bug reports
in [GitHub issues](https://github.com/MBaysanLab/cosap/issues).


---
