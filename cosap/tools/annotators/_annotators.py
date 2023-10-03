from abc import ABC, abstractmethod
from subprocess import run


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
        output_vcf = f"{input_without_ext}_chr1_22.vcf"
        command = [
            "bcftools",
            "view",
            vcf_path,
            "--targets",
            chromosomes,
            "-o",
            output_vcf,
        ]
        run(command)
        return output_vcf
