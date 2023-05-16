import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class DeepVariantVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def create_make_examples_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ):
        
        input_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        sample_name = caller_config[VariantCallingKeys.GERMLINE_SAMPLE_NAME]
        tmpdir = os.path.dirname(input_bam)

        command = [
            "dv_make_examples.py ",
            "--cores",
            str(AppConfig.MAX_CORES_PER_JOB),
            "--ref",
            library_paths.REF_FASTA,
            "--reads",
            input_bam,
            "--sample",
            sample_name,
            "--examples",
            tmpdir,
            "--logdir",
            tmpdir,
            "--gvcf",
            tmpdir
        ]
        return command

    @classmethod
    def create_call_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ):
        input_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        sample_name = caller_config[VariantCallingKeys.GERMLINE_SAMPLE_NAME]
        tmpdir = os.path.dirname(input_bam)
        outfile = os.path.join(tmpdir, f"{sample_name}.tmp")

        command = [
            "dv_call_variants.py",
            "--cores",
            str(AppConfig.MAX_CORES_PER_JOB),
            "--outfile",
            outfile,
            "--sample",
            sample_name,
            "--examples",
            tmpdir,
            "--model",
            "wgs"
        ]
        return command

    @classmethod
    def create_postprocess_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ):

        input_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        sample_name = caller_config[VariantCallingKeys.GERMLINE_SAMPLE_NAME]
        tmpdir = os.path.dirname(input_bam)
        infile = os.path.join(tmpdir, f"{sample_name}.tmp")
        outfile = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
        gvcf_infile = os.path.join(tmpdir, f"{sample_name}.gvcf.tfrecord@{AppConfig.MAX_THREADS_PER_JOB}.gz")
        gvcf_outfile = caller_config[VariantCallingKeys.GVCF_OUTPUT]

        command = [
            "dv_postprocess_variants.py",
            "--ref",
            library_paths.REF_FASTA,
            "--infile",
            infile,
            "--outfile",
            outfile,
            "--gvcf_infile",
            gvcf_infile,
            "--gvcf_outfile",
            gvcf_outfile,
            ]
    
        return command


    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()

        make_examples_command = cls.create_make_examples_command(
            caller_config, library_paths
        )
        call_variants_command = cls.create_call_variants_command(
            caller_config, library_paths
        )
        postprocess_variants_command = cls.create_postprocess_variants_command(
            caller_config, library_paths
        )
        run(make_examples_command)
        run(call_variants_command)
        run(postprocess_variants_command)
