from abc import ABC, abstractmethod
from subprocess import run
import os
from ..._pipeline_config import AnnotatorKeys


class _Annotator(ABC):
    @abstractmethod
    def annotate(self):
        pass


class _Annotatable:
    @classmethod
    def chr_filter_vcf(cls, vcf_path: str) -> str:
        """
        Only get chromosomes 1-22 from given vcf and returns path of filtered vcf.
        """

        chromosomes = ",".join([f"chr{i}" for i in list(range(1, 23))])
        input_without_ext = "".join(vcf_path.split(".")[:-1])
        input_extension = vcf_path.split(".")[-1]
        output_file = f"{input_without_ext}_chr1_22.{input_extension}"

        command = [
            "bcftools",
            "view",
            vcf_path,
            "--targets",
            chromosomes,
            "-o",
            output_file,
        ]
        run(command)
        return output_file

    @classmethod
    def create_output_dir(cls, annotator_config: dict, workdir:str = None) -> str:
        """
        Creates output directory for annotator.
        """
        output_dir = os.path.dirname(annotator_config[AnnotatorKeys.OUTPUT])
        if workdir:
            output_dir = os.path.join(workdir, output_dir)
        else:
            output_dir = os.path.abspath(output_dir)

        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)
        return output_dir
