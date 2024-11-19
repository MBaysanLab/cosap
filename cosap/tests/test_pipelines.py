import os

import pandas as pd
import pytest

from cosap import (
    MDUP,
    AnnotatorFactory,
    BamReader,
    CNVCallerFactory,
    FastqReader,
    Mapper,
    MapperFactory,
    MSICallerFactory,
    Pipeline,
    PipelineRunner,
    PreprocessorFactory,
    VariantCaller,
    VariantCallerFactory,
)
from cosap._pipeline_config import MappingKeys, VariantCallingKeys
from cosap._utils import join_paths
from cosap.tests import create_test_fastqs_from_chr1_ref_with_1_snp

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


@pytest.mark.integration
def test_bowtie_somaticsniper_pipeline(tmp_path):
    workdir = str(tmp_path)

    # Change read/write permissions for the test
    os.chmod(workdir, 0o777)

    (
        normal_fastqs,
        tumor_fastqs,
        single_SNP_tuple,
    ) = create_test_fastqs_from_chr1_ref_with_1_snp(workdir=workdir)
    snp_chrom, snp_pos, snp_ref, snp_alt = single_SNP_tuple

    run_basic_bowtie_somaticsniper_pipeline(
        workdir=workdir, normal_fastqs=normal_fastqs, tumor_fastqs=tumor_fastqs
    )

    vcf_path = os.path.join(workdir, "VCF", "somaticsniper", "all_TEST_VCF.vcf")
    assert os.path.isfile(vcf_path)

    somaticsniper_headers = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "NORMAL",
        "TUMOR",
    ]
    vcf_body = pd.read_csv(vcf_path, sep="\t", comment="#", names=somaticsniper_headers)
    print(vcf_body[["CHROM", "POS", "REF", "ALT"]])
    assert vcf_body[["CHROM", "POS", "REF", "ALT"]].equals(
        pd.DataFrame(
            {
                "CHROM": [snp_chrom],
                "POS": [snp_pos],
                "REF": [snp_ref],
                "ALT": [snp_alt],
            }
        )
    )


def run_basic_bowtie_somaticsniper_pipeline(
    workdir: str, normal_fastqs: list, tumor_fastqs: list
):
    normal_fastqs = [
        FastqReader(os.path.join(workdir, fastq), read=i)
        for i, fastq in enumerate(normal_fastqs, start=1)
    ]
    tumor_fastqs = [
        FastqReader(os.path.join(workdir, fastq), read=i)
        for i, fastq in enumerate(tumor_fastqs, start=1)
    ]

    mapper_params = {
        "read_groups": {
            "ID": "0",
            "SM": "TEMP_VAL",
            "PU": "0",
            "PL": "il",
            "LB": "0",
        }
    }

    mapper_normal = Mapper(
        input_step=normal_fastqs,
        library="bowtie",
        params={**mapper_params, "SM": "normal"},
    )
    mapper_tumor = Mapper(
        input_step=tumor_fastqs,
        library="bowtie",
        params={**mapper_params, "SM": "tumor"},
    )

    mdup_normal = MDUP(input_step=mapper_normal)
    mdup_tumor = MDUP(input_step=mapper_tumor)

    caller = VariantCaller(
        library="somaticsniper",
        normal_sample=mdup_normal,
        tumor_sample=mdup_tumor,
        name="TEST_VCF",
    )

    pipeline = (
        Pipeline()
        .add(mapper_normal)
        .add(mapper_tumor)
        .add(mdup_normal)
        .add(mdup_tumor)
        .add(caller)
    )

    config = pipeline.build(workdir=workdir)

    pipeline_runner = PipelineRunner()
    pipeline_runner.run_pipeline(config)


@pytest.fixture
def mapper(library):
    return MapperFactory.create(library)


@pytest.fixture
def mapper_builder(library):
    return MapperFactory.create(library, builder=True)


@pytest.fixture
def preprocessor(library):
    return PreprocessorFactory.create(library)


@pytest.fixture
def variant_caller(library):
    return VariantCallerFactory.create(library)


@pytest.fixture
def annotator(library):
    return AnnotatorFactory.create(library)


@pytest.fixture
def cnv_caller(library):
    return CNVCallerFactory.create(library)


@pytest.fixture
def msi_caller(library):
    return MSICallerFactory.create(library)


@pytest.mark.parametrize(
    "library",
    ["bowtie", "bwa2", "bwa"],
)
def test_mapper(mapper, library, tmp_path):
    assert mapper is not None

    # Create mapper object with builder pattern
    fastqs = [
        FastqReader(
            join_paths(TEST_DATA_DIR, "test_germline_1.fastq.gz"), read=1, name="test"
        ),
        FastqReader(
            join_paths(TEST_DATA_DIR, "test_germline_2.fastq.gz"), read=2, name="test"
        ),
    ]
    mapper_builder = Mapper(
        input_step=fastqs,
        library=library,
        params={
            "read_groups": {
                "ID": "0",
                "SM": "test",
                "PU": "0",
                "PL": "ILLUMINA",
                "LB": "0",
            }
        },
        output_dir=str(tmp_path),
    )
    mapper_config = mapper_builder.get_config()
    mapper.map(mapper_config[mapper_builder.key][mapper_builder.name])

    # Check if output file exists
    assert os.path.isfile(mapper_builder.get_output())
    # Check if file is not empty
    assert os.path.getsize(mapper_builder.get_output()) > 0


@pytest.mark.parametrize("library", ["mutect2"])
def test_tumor_only_variantcaller(variant_caller, library, tmp_path):
    assert variant_caller is not None

    # Create variant caller object with builder pattern
    bam = BamReader(join_paths(TEST_DATA_DIR, "tumor_test.bam"), name="tumor")

    variant_caller_builder = VariantCaller(
        library=library,
        tumor_sample=bam,
        params={
            "tumor_sample_name": "tumor",
        },
        output_dir=str(tmp_path),
    )
    variant_caller_config = variant_caller_builder.get_config()
    variant_caller.call_variants(
        variant_caller_config[variant_caller_builder.key][variant_caller_builder.name]
    )

    # Check if output file exists
    assert os.path.isfile(variant_caller_builder.get_output())
    # Check if file is not empty
    assert os.path.getsize(variant_caller_builder.get_output()) > 0


@pytest.mark.parametrize(
    "library",
    [
        "mutect2",
        # "somaticsniper",
        # "varscan",
        # "vardict",
        # "varnet",
    ],
)
def test_tumor_normal_variantcaller(variant_caller, library, tmp_path):
    assert variant_caller is not None

    # Create variant caller object with builder pattern
    tumor_bam = BamReader(
        join_paths(TEST_DATA_DIR, "tumor_test_small.bam"), name="test_tumor"
    )
    normal_bam = BamReader(
        join_paths(TEST_DATA_DIR, "normal_test_small.bam"), name="test_normal"
    )

    variant_caller_builder = VariantCaller(
        library=library,
        tumor_sample=tumor_bam,
        normal_sample=normal_bam,
        params={
            "tumor_sample_name": "test_tumor",
        },
        output_dir=str(tmp_path),
        bed_file=join_paths(TEST_DATA_DIR, "test.bed"),
    )
    variant_caller_config = variant_caller_builder.get_config()
    variant_caller.call_variants(
        variant_caller_config[variant_caller_builder.key][variant_caller_builder.name]
    )

    # Check if output file exists
    assert os.path.isfile(variant_caller_builder.get_output())
    # Check if file is not empty
    assert os.path.getsize(variant_caller_builder.get_output()) > 0


@pytest.mark.parametrize(
    "library",
    [
        #"haplotypecaller",
        "deepvariant",
    ],
)
def test_germline_variantcaller_gvcf(variant_caller, library, tmp_path):
    assert variant_caller is not None

    # Create variant caller object with builder pattern
    bam = BamReader(join_paths(TEST_DATA_DIR, "normal_test_small.bam"), name="normal")

    variant_caller_builder = VariantCaller(
        library=library, normal_sample=bam, output_dir=str(tmp_path), gvcf=True, bed_file="chr1:100000-200000"
    )
    variant_caller_config = variant_caller_builder.get_config()
    variant_caller.call_variants(
        variant_caller_config[variant_caller_builder.key][variant_caller_builder.name]
    )

    # Check if output file exists
    assert os.path.isfile(variant_caller_builder.get_output())
    # Check if file is not empty
    assert os.path.getsize(variant_caller_builder.get_output()) > 0
