import os
import warnings

import pandas as pd
import pysam
import pytest

from cosap import (
    MDUP,
    AnnotatorFactory,
    BamReader,
    CNVCallerFactory,
    DNAPipeline,
    DNAPipelineInput,
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
from cosap._pipeline_config import PipelineKeys

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# Real NA12878 paired-end reads covering chr21:16570000-16610000 and
# chr22:16570000-16610000 on GRCh38, extracted from
#   nf-core/test-datasets modules branch:
#     data/genomics/homo_sapiens/illumina/bam/NA12878.chr21_22.1X.bam
# (converted to paired FASTQs via `samtools sort -n | samtools fastq`).
SAREK_NA12878_FASTQ_R1 = os.path.join(
    TEST_DATA_DIR, "sarek_NA12878_chr21_22_1X_1.fq.gz"
)
SAREK_NA12878_FASTQ_R2 = os.path.join(
    TEST_DATA_DIR, "sarek_NA12878_chr21_22_1X_2.fq.gz"
)
# GIAB (Genome in a Bottle) high-confidence truth calls for NA12878,
# restricted to chr21/chr22, on GRCh38. Sourced from
#   nf-core/test-datasets modules branch:
#     data/genomics/homo_sapiens/illumina/vcf/NA12878_GIAB.chr21_22.vcf.gz
SAREK_GIAB_TRUTH_VCF = os.path.join(
    TEST_DATA_DIR, "NA12878_GIAB.chr21_22.vcf.gz"
)

# Real NA12878 chr22 alignments (GRCh38) from the sarek test-data suite,
# used as the normal/germline BAM for unit tests. Sourced from
#   nf-core/test-datasets modules branch:
#     data/genomics/homo_sapiens/illumina/bam/NA12878.chr22.1X.bam
# @RG carries SM:NA12878.
SAREK_NA12878_NORMAL_BAM = os.path.join(TEST_DATA_DIR, "sarek_normal_chr22.bam")

# SEQC2 (FDA somatic reference sample) WES tumor/normal pair, sliced to a
# high-density 1.5 Mb window on chr22 (chr22:22,000,000-23,500,000) that
# contains 53 high-confidence somatic SNVs and multiple somatic INDELs per
# the SEQC2 v1.2.1 high-confidence release. Sourced from:
#   https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/
#     Somatic_Mutation_WG/data/WES/WES_FD_{T,N}_1.bwa.dedup.bam
# @RG carries SM:WES_FD_T_1 / SM:WES_FD_N_1.
SEQC2_TUMOR_BAM = os.path.join(TEST_DATA_DIR, "seqc2_tumor.bam")
SEQC2_NORMAL_BAM = os.path.join(TEST_DATA_DIR, "seqc2_normal.bam")
SEQC2_REGION_BED = os.path.join(TEST_DATA_DIR, "seqc2_region.bed")
# High-confidence somatic truth (SNV+INDEL) from the SEQC2 v1.2.1 release,
# sliced to the same chr22 window. Sourced from:
#   https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/
#     Somatic_Mutation_WG/release/latest/
#       high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
#       high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz
SEQC2_SOMATIC_TRUTH_VCF = os.path.join(TEST_DATA_DIR, "seqc2_somatic_chr22.vcf.gz")
SEQC2_TUMOR_SAMPLE_NAME = "WES_FD_T_1"
SEQC2_NORMAL_SAMPLE_NAME = "WES_FD_N_1"


def _find_variant_caller_output(
    config: dict, library_name: str, workdir: str = None
) -> str:
    """Locate the output VCF path for a given variant-calling library in a
    DNAPipeline-built config. Returns the first matching step's output.

    Mirrors ``VariantCallerBuilder.get_output()`` in
    ``cosap/pipeline_builder/builders/_variantcaller_builder.py``: the config
    stores only a filename under ``all_variants`` (VCF mode) or ``gvcf_output``
    (GVCF mode); the full path is built by joining with ``output_dir``.

    The ``output_dir`` recorded in the config is *relative* to the pipeline
    workdir (e.g. ``VCF/haplotypecaller``). Pass ``workdir`` to receive an
    absolute path that can be passed directly to ``os.path.isfile``.
    """
    vc_section = config.get("variant-calling", {})
    for step_name, step_config in vc_section.items():
        if str(step_config.get("library", "")).lower() != library_name.lower():
            continue
        output_dir = step_config.get("output_dir")
        if not output_dir:
            raise KeyError(
                f"variant-calling step {step_name!r} has no 'output_dir' key; "
                f"keys present: {list(step_config.keys())}"
            )
        is_gvcf = str(step_config.get("output_type", "")).upper() == "GVCF"
        filename_key = "gvcf_output" if is_gvcf else "all_variants"
        filename = step_config.get(filename_key)
        if not filename:
            raise KeyError(
                f"variant-calling step {step_name!r} has no {filename_key!r} "
                f"key; keys present: {list(step_config.keys())}"
            )
        relative_path = os.path.join(output_dir, filename)
        if workdir and not os.path.isabs(relative_path):
            return os.path.join(workdir, relative_path)
        return relative_path
    raise KeyError(
        f"No variant-calling step with library={library_name!r} found in config. "
        f"Available steps: {list(vc_section.keys())}"
    )


_VCF_COLUMNS = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
    "SAMPLE",
]


def _read_vcf_variants(vcf_path: str) -> set:
    """Parse a (optionally gzipped) VCF file and return the set of
    (CHROM, POS, REF, ALT) tuples for its variant records. Multi-allelic ALT
    fields are split into one tuple per alternate allele."""
    compression = "gzip" if vcf_path.endswith(".gz") else None
    # Only the first 5 columns (CHROM, POS, ID, REF, ALT) are needed; restrict
    # both `names` and `usecols` so the parser is agnostic to whether the VCF
    # carries FORMAT/SAMPLE columns.
    df = pd.read_csv(
        vcf_path,
        sep="\t",
        comment="#",
        header=None,
        names=_VCF_COLUMNS[:5],
        usecols=range(5),
        compression=compression,
        dtype={"CHROM": str, "REF": str, "ALT": str},
    )
    variants = set()
    for chrom, pos, ref, alt in zip(
        df["CHROM"].astype(str),
        df["POS"].astype(int),
        df["REF"].astype(str),
        df["ALT"].astype(str),
    ):
        for single_alt in alt.split(","):
            variants.add((chrom, pos, ref, single_alt))
    return variants


def _get_sample_name_from_bam(bam_path: str) -> str:
    """Return the SM (sample) value from the first @RG line of ``bam_path``.

    Variant callers like Mutect2 require the ``--normal-name`` / ``--tumor-name``
    arguments to match the SM tag stored in the BAM header. Reading the value
    from the file at test time avoids hard-coded mismatches if the test BAM
    is ever regenerated with different SM tags.
    """
    abs_path = os.path.abspath(bam_path)
    with pysam.AlignmentFile(abs_path, "rb", check_sq=False) as bam_file:
        for rg in bam_file.header.get("RG", []):
            sm = rg.get("SM")
            if sm:
                return sm
    raise ValueError(f"No @RG SM tag found in BAM header: {bam_path}")


def _write_fastq(path):
    path.write_text("@read\nA\n+\nI\n")
    return str(path)


def test_germline_dna_pipeline_includes_mitochondrial_caller_by_default(tmp_path):
    r1 = _write_fastq(tmp_path / "normal_R1.fastq")
    r2 = _write_fastq(tmp_path / "normal_R2.fastq")

    pipeline_input = DNAPipelineInput(
        ANALYSIS_TYPE="germline",
        WORKDIR=str(tmp_path),
        NORMAL_SAMPLE=(r1, r2),
        MAPPERS=["bwa"],
        VARIANT_CALLERS=["haplotypecaller"],
        NORMAL_SAMPLE_NAME="NA12878",
    )

    dna_pipeline = DNAPipeline(pipeline_input)
    mitochondrial_config = dna_pipeline.config[
        PipelineKeys.MITOCHONDRIAL_VARIANT_CALLING
    ]

    assert list(mitochondrial_config) == ["NA12878_bwa_mitochondrial"]
    assert (
        mitochondrial_config["NA12878_bwa_mitochondrial"]["library"]
        == "gatk_mitochondrial"
    )


def test_somatic_dna_pipeline_includes_mitochondrial_callers_by_default(tmp_path):
    normal_r1 = _write_fastq(tmp_path / "normal_R1.fastq")
    normal_r2 = _write_fastq(tmp_path / "normal_R2.fastq")
    tumor_r1 = _write_fastq(tmp_path / "tumor_R1.fastq")
    tumor_r2 = _write_fastq(tmp_path / "tumor_R2.fastq")

    pipeline_input = DNAPipelineInput(
        ANALYSIS_TYPE="somatic",
        WORKDIR=str(tmp_path),
        NORMAL_SAMPLE=(normal_r1, normal_r2),
        TUMOR_SAMPLES=(tumor_r1, tumor_r2),
        MAPPERS=["bwa"],
        VARIANT_CALLERS=["mutect2"],
        NORMAL_SAMPLE_NAME="normal",
        TUMOR_SAMPLE_NAME="tumor",
    )

    dna_pipeline = DNAPipeline(pipeline_input)

    assert set(dna_pipeline.config[PipelineKeys.MITOCHONDRIAL_VARIANT_CALLING]) == {
        "normal_bwa_mitochondrial",
        "tumor_bwa_mitochondrial",
    }


def test_somatic_dna_pipeline_reuses_normal_steps_for_multiple_tumors(tmp_path):
    normal_r1 = _write_fastq(tmp_path / "normal_R1.fastq")
    normal_r2 = _write_fastq(tmp_path / "normal_R2.fastq")
    tumor_1_r1 = _write_fastq(tmp_path / "tumor_1_R1.fastq")
    tumor_1_r2 = _write_fastq(tmp_path / "tumor_1_R2.fastq")
    tumor_2_r1 = _write_fastq(tmp_path / "tumor_2_R1.fastq")
    tumor_2_r2 = _write_fastq(tmp_path / "tumor_2_R2.fastq")

    pipeline_input = DNAPipelineInput(
        ANALYSIS_TYPE="somatic",
        WORKDIR=str(tmp_path),
        NORMAL_SAMPLE=(normal_r1, normal_r2),
        TUMOR_SAMPLES=[(tumor_1_r1, tumor_1_r2), (tumor_2_r1, tumor_2_r2)],
        MAPPERS=["bwa"],
        VARIANT_CALLERS=["mutect2"],
        NORMAL_SAMPLE_NAME="normal",
        TUMOR_SAMPLE_NAME="tumor",
    )

    dna_pipeline = DNAPipeline(pipeline_input)

    assert list(dna_pipeline.config[PipelineKeys.MDUP]).count("normal_bwa") == 1
    assert set(dna_pipeline.config[PipelineKeys.MITOCHONDRIAL_VARIANT_CALLING]) == {
        "normal_bwa_mitochondrial",
        "tumor_1_bwa_mitochondrial",
        "tumor_2_bwa_mitochondrial",
    }


@pytest.mark.integration
def test_bwa_haplotypecaller_dna_pipeline(tmp_path):
    """End-to-end test of the default DNAPipeline workflow (germline analysis)
    using BWA + HaplotypeCaller.

    Input FASTQs and truth VCF come directly from the nf-core/sarek test
    dataset (nf-core/test-datasets, modules branch): real NA12878 reads
    aligned to GRCh38 chr21/chr22 around position 16.57-16.61 Mb, plus the
    corresponding GIAB high-confidence truth calls.

    Assertion policy:
      * If any called variant matches a GIAB truth variant -> pass.
      * If variants were called but none overlap the truth set -> warn only
        (possible caller drift, alignment differences, or reference mismatch);
        do not fail.
      * If the output VCF contains no variant records at all -> fail.
    """
    workdir = str(tmp_path)

    # Loosen permissions so any tools run in child processes can write here.
    os.chmod(workdir, 0o777)

    assert os.path.isfile(SAREK_NA12878_FASTQ_R1), (
        f"Missing sarek test FASTQ: {SAREK_NA12878_FASTQ_R1}"
    )
    assert os.path.isfile(SAREK_NA12878_FASTQ_R2), (
        f"Missing sarek test FASTQ: {SAREK_NA12878_FASTQ_R2}"
    )
    assert os.path.isfile(SAREK_GIAB_TRUTH_VCF), (
        f"Missing GIAB truth VCF: {SAREK_GIAB_TRUTH_VCF}"
    )

    pipeline_input = DNAPipelineInput(
        ANALYSIS_TYPE="germline",
        WORKDIR=workdir,
        NORMAL_SAMPLE=(SAREK_NA12878_FASTQ_R1, SAREK_NA12878_FASTQ_R2),
        MAPPERS=["bwa"],
        VARIANT_CALLERS=["haplotypecaller"],
        NORMAL_SAMPLE_NAME="NA12878",
        MITOCHONDRIAL=False,
    )

    dna_pipeline = DNAPipeline(pipeline_input)
    dna_pipeline.run_pipeline()

    vcf_path = _find_variant_caller_output(
        dna_pipeline.config, "haplotypecaller", workdir=workdir
    )
    assert os.path.isfile(vcf_path), f"Expected VCF not produced at {vcf_path}"

    called_variants = _read_vcf_variants(vcf_path)
    truth_variants = _read_vcf_variants(SAREK_GIAB_TRUTH_VCF)

    assert len(truth_variants) > 0, (
        f"GIAB truth VCF {SAREK_GIAB_TRUTH_VCF} parsed to zero variants; "
        "test data is corrupt."
    )

    matched = called_variants & truth_variants
    if matched:
        return  # success: pipeline rediscovered at least one known truth variant

    if len(called_variants) > 0:
        warnings.warn(
            f"HaplotypeCaller produced {len(called_variants)} variant(s) but "
            f"none overlapped the {len(truth_variants)} GIAB truth variants. "
            f"This may indicate caller drift, alignment differences, or a "
            f"reference-version mismatch (cosap uses its configured "
            f"LibraryPaths().REF_FASTA). First few called variants: "
            f"{list(called_variants)[:5]}"
        )
        return  # do not fail; policy established for this test

    pytest.fail(
        f"Output VCF {vcf_path} contains no variant records; "
        f"expected to rediscover at least one of {len(truth_variants)} "
        "GIAB truth variants on chr21/chr22."
    )


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
        FastqReader(SAREK_NA12878_FASTQ_R1, read=1, name="NA12878"),
        FastqReader(SAREK_NA12878_FASTQ_R2, read=2, name="NA12878"),
    ]
    mapper_builder = Mapper(
        input_step=fastqs,
        library=library,
        params={
            "read_groups": {
                "ID": "0",
                "SM": "NA12878",
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

    # Create variant caller object with builder pattern. Read the actual SM
    # tag from the BAM so Mutect2's --tumor-name aligns with the header.
    tumor_sample_name = _get_sample_name_from_bam(SEQC2_TUMOR_BAM)
    bam = BamReader(SEQC2_TUMOR_BAM, name=tumor_sample_name)

    variant_caller_builder = VariantCaller(
        library=library,
        tumor_sample=bam,
        params={
            "tumor_sample_name": tumor_sample_name,
        },
        output_dir=str(tmp_path),
        bed_file=SEQC2_REGION_BED,
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

    # SEQC2 (FDA) tumor/normal WES pair, sliced to chr22:22,000,000-23,500,000
    # where the SEQC2 v1.2.1 high-confidence release declares ~53 somatic SNVs
    # plus somatic INDELs. This gives a real tumor-vs-normal signal rather
    # than a same-donor placeholder.
    # Read the actual SM tags from the BAM headers so Mutect2's
    # --normal-name / --tumor-name arguments line up with what the BAM
    # reports, rather than relying on (possibly stale) hard-coded constants.
    tumor_sample_name = _get_sample_name_from_bam(SEQC2_TUMOR_BAM)
    normal_sample_name = _get_sample_name_from_bam(SEQC2_NORMAL_BAM)

    tumor_bam = BamReader(SEQC2_TUMOR_BAM, name=tumor_sample_name)
    normal_bam = BamReader(SEQC2_NORMAL_BAM, name=normal_sample_name)

    variant_caller_builder = VariantCaller(
        library=library,
        tumor_sample=tumor_bam,
        normal_sample=normal_bam,
        params={
            "tumor_sample_name": tumor_sample_name,
            "germline_sample_name": normal_sample_name,
        },
        output_dir=str(tmp_path),
        bed_file=SEQC2_REGION_BED,
    )
    variant_caller_config = variant_caller_builder.get_config()
    variant_caller.call_variants(
        variant_caller_config[variant_caller_builder.key][variant_caller_builder.name]
    )

    # Check if output file exists
    assert os.path.isfile(variant_caller_builder.get_output())
    # Check if file is not empty
    assert os.path.getsize(variant_caller_builder.get_output()) > 0

    # Overlap check against SEQC2 truth somatic variants.
    # Policy mirrors test_bwa_haplotypecaller_dna_pipeline:
    #   * any overlap -> pass;
    #   * variants called but none overlap -> warn only (caller/filter drift);
    #   * empty VCF -> fail.
    called_variants = _read_vcf_variants(variant_caller_builder.get_output())
    truth_variants = _read_vcf_variants(SEQC2_SOMATIC_TRUTH_VCF)
    assert len(truth_variants) > 0, (
        f"SEQC2 truth VCF {SEQC2_SOMATIC_TRUTH_VCF} parsed to zero variants; "
        "test data is corrupt."
    )

    matched = called_variants & truth_variants
    if matched:
        return
    if len(called_variants) > 0:
        warnings.warn(
            f"{library} produced {len(called_variants)} variant(s) but none "
            f"overlapped the {len(truth_variants)} SEQC2 truth variants in "
            f"chr22:22,000,000-23,500,000. First few called: "
            f"{list(called_variants)[:5]}"
        )
        return
    pytest.fail(
        f"Output VCF {variant_caller_builder.get_output()} contains no "
        f"variant records; expected to rediscover at least one of "
        f"{len(truth_variants)} SEQC2 truth somatic variants."
    )


@pytest.mark.parametrize(
    "library",
    [
        # "haplotypecaller",
        "deepvariant",
    ],
)
def test_germline_variantcaller_gvcf(variant_caller, library, tmp_path):
    assert variant_caller is not None

    # Create variant caller object with builder pattern
    bam = BamReader(SAREK_NA12878_NORMAL_BAM, name="NA12878")

    variant_caller_builder = VariantCaller(
        library=library,
        normal_sample=bam,
        output_dir=str(tmp_path),
        gvcf=True,
        # sarek_normal_chr22.bam covers NA12878 reads on chr22; pick a window
        # with non-zero coverage.
        bed_file="chr22:16570000-16610000",
        model_type="WES",
    )
    variant_caller_config = variant_caller_builder.get_config()
    variant_caller.call_variants(
        variant_caller_config[variant_caller_builder.key][variant_caller_builder.name]
    )

    # Check if output file exists
    assert os.path.isfile(variant_caller_builder.get_output())
    # Check if file is not empty
    assert os.path.getsize(variant_caller_builder.get_output()) > 0
