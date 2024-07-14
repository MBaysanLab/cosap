import os
import pandas as pd
import pytest
from cosap import Pipeline, PipelineRunner, FastqReader, Mapper, VariantCaller, MDUP
from cosap.tests.pipeline_integration import create_test_fastqs_from_chr1_ref_with_1_snp


@pytest.mark.pipeline
@pytest.mark.takes_minutes
def test_bowtie_somaticsniper_pipeline(tmp_path):
    workdir = str(tmp_path)

    normal_fastqs, tumor_fastqs, single_SNP_tuple = create_test_fastqs_from_chr1_ref_with_1_snp(workdir=workdir)
    snp_chrom, snp_pos, snp_ref, snp_alt = single_SNP_tuple

    run_basic_bowtie_somaticsniper_pipeline(workdir=workdir, normal_fastqs=normal_fastqs, tumor_fastqs=tumor_fastqs)

    vcf_path = os.path.join(workdir, "VCF", "somaticsniper", "all_TEST_VCF.vcf")
    assert os.path.isfile(vcf_path)

    somaticsniper_headers = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR"]
    vcf_body = pd.read_csv(vcf_path, sep="\t", comment="#", names=somaticsniper_headers)
    print(vcf_body[["CHROM", "POS", "REF", "ALT"]])
    assert vcf_body[["CHROM", "POS", "REF", "ALT"]].equals(pd.DataFrame({
        "CHROM": [snp_chrom],
        "POS": [snp_pos],
        "REF": [snp_ref],
        "ALT": [snp_alt],
    }))


def run_basic_bowtie_somaticsniper_pipeline(workdir: str, normal_fastqs: list, tumor_fastqs: list):
    normal_fastqs = [FastqReader(os.path.join(workdir, fastq), read=i) for i, fastq in enumerate(normal_fastqs, start=1)]
    tumor_fastqs = [FastqReader(os.path.join(workdir, fastq), read=i) for i, fastq in enumerate(tumor_fastqs, start=1)]

    mapper_params = {"read_groups": {
                        "ID": "0",
                        "SM": "TEMP_VAL",
                        "PU": "0",
                        "PL": "il",
                        "LB": "0",
                    }}

    mapper_normal = Mapper(input_step=normal_fastqs, library="bowtie", params={**mapper_params, "SM": "normal"})
    mapper_tumor = Mapper(input_step=tumor_fastqs, library="bowtie", params={**mapper_params, "SM": "tumor"})

    mdup_normal = MDUP(input_step=mapper_normal)
    mdup_tumor = MDUP(input_step=mapper_tumor)

    caller = VariantCaller(library="somaticsniper", germline=mdup_normal, tumor=mdup_tumor, name="TEST_VCF")

    pipeline = (Pipeline()
        .add(mapper_normal)
        .add(mapper_tumor)
        .add(mdup_normal)
        .add(mdup_tumor)
        .add(caller)
    )
                
    config = pipeline.build(workdir=workdir)

    pipeline_runner = PipelineRunner()
    pipeline_runner.run_pipeline(config)
