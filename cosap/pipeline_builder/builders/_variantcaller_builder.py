from copy import copy
from dataclasses import dataclass, field
from subprocess import PIPE, STDOUT, Popen
from typing import Dict

from ..._formats import FileFormats
from ..._pipeline_config import (DefaultValues, MappingKeys, PipelineKeys,
                                 VariantCallingKeys)
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class VariantCaller(_IPipelineStep, _PipelineStep):
    PY2_PACKAGES = ["strelka"]

    library: str
    params: dict = field(default_factory=dict)
    germline: str = None
    tumor: str = None
    name: str = None
    gvcf: bool = False
    key: str = PipelineKeys.VARIANT_CALLING

    def __post_init__(self):
        if self.name is None:

            name_temp = []
            if self.germline:
                name_temp.append(self.germline.name)
            if self.tumor:
                name_temp.append(self.tumor.name)
            name_temp.append(self.library)

            self.name = "_".join(name_temp)

        if VariantCallingKeys.GERMLINE_SAMPLE_NAME not in self.params and self.germline:
            self.params[
                VariantCallingKeys.GERMLINE_SAMPLE_NAME
            ] = self._get_sample_name_from_bam(self.germline.get_output())

    def _get_sample_name_from_bam(self, bam) -> str:
        cmd = f"samtools view -H {bam} | grep '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/g' | uniq"
        ps = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
        sample_name = str(ps.communicate()[0])
        return str(sample_name)

    def get_output(self):
        config = self.get_config()
        if self.gvcf:
            return config[self.key][self.name][VariantCallingKeys.GVCF_OUTPUT] 
        return config[self.key][self.name][VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]

    def get_config(self) -> Dict:
        unfiltered_variants_output_filename = FileFormats.GATK_UNFILTERED_OUTPUT.format(
            identification=self.name
        )
        filtered_variants_output_filename = FileFormats.GATK_FILTERED_OUTPUT.format(
            identification=self.name
        )
        snp_output_filename = FileFormats.GATK_SNP_OUTPUT.format(
            identification=self.name
        )
        indel_output_filename = FileFormats.GATK_INDEL_OUTPUT.format(
            identification=self.name
        )
        other_variants_output_filename = FileFormats.GATK_OTHER_VARIANTS_OUTPUT.format(
            identification=self.name
        )
        gvcf_output_filename = FileFormats.GATK_GVCF_OUTPUT.format(
            identification=self.name
        )

        vc_config = {
            self.name: {
                VariantCallingKeys.LIBRARY: self.library,
                VariantCallingKeys.PARAMS: self.params,
                VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT: unfiltered_variants_output_filename,
                VariantCallingKeys.FILTERED_VARIANTS_OUTPUT: filtered_variants_output_filename,
                VariantCallingKeys.SNP_OUTPUT: snp_output_filename,
                VariantCallingKeys.INDEL_OUTPUT: indel_output_filename,
                VariantCallingKeys.OTHER_VARIANTS_OUTPUT: other_variants_output_filename,
                VariantCallingKeys.GVCF_OUTPUT: gvcf_output_filename
            },
        }

        if self.tumor is not None:
            vc_config[self.name][
                VariantCallingKeys.TUMOR_INPUT
            ] = self.tumor.get_output()

        if self.germline is not None:
            vc_config[self.name][
                VariantCallingKeys.GERMLINE_INPUT
            ] = self.germline.get_output()

        config = {self.key: vc_config}
        return config
