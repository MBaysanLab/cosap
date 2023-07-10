import os
from copy import copy
from dataclasses import dataclass, field
from subprocess import PIPE, STDOUT, Popen
from typing import Dict

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import PipelineKeys, VariantCallingKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class VariantCaller(_IPipelineStep, _PipelineStep):
    PY2_PACKAGES = ["strelka", "manta"]

    library: str
    params: dict = field(default_factory=dict)
    germline: str = None
    tumor: _PipelineStep = None
    name: _PipelineStep = None
    gvcf: bool = False
    bed_file: str = None
    key: str = PipelineKeys.VARIANT_CALLING
    next_step: _PipelineStep = None

    def __post_init__(self):
        if self.name is None:

            name_temp = []
            if self.germline:
                name_temp.append(self.germline.name)
            if self.tumor:
                name_temp.append(self.tumor.name)
            name_temp.append(self.library)

            self.name = "_".join(name_temp)

        # TODO: Read sample names from bam.
        if VariantCallingKeys.GERMLINE_SAMPLE_NAME not in self.params and self.germline:
            self.params[VariantCallingKeys.GERMLINE_SAMPLE_NAME] = "normal"

        if VariantCallingKeys.TUMOR_SAMPLE_NAME not in self.params and self.tumor:
            self.params[VariantCallingKeys.TUMOR_SAMPLE_NAME] = "tumor"

        if self.germline:
            self.germline.next_step = self
        if self.tumor:
            self.tumor.next_step = self

    def _get_sample_name_from_bam(self, bam) -> str:
        abs_bam_path = os.path.abspath(os.path.normpath(bam))
        cmd = f"samtools view -H {abs_bam_path} | grep '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/g' | uniq"
        ps = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
        sample_name = str(ps.communicate()[0])
        return str(sample_name)

    def get_output(self):
        config = self.get_config()
        if self.gvcf:
            return config[self.key][self.name][VariantCallingKeys.GVCF_OUTPUT]
        return config[self.key][self.name][VariantCallingKeys.ALL_VARIANTS_OUTPUT]

    def get_config(self) -> Dict:
        unfiltered_variants_output_filename = FileFormats.GATK_UNFILTERED_OUTPUT.format(
            identification=self.name
        )
        all_variants_output_filename = FileFormats.ALL_VARIANTS_OUTPUT.format(
            identification=self.name
        )
        snp_output_filename = FileFormats.SNP_OUTPUT.format(identification=self.name)
        indel_output_filename = FileFormats.INDEL_OUTPUT.format(
            identification=self.name
        )
        other_variants_output_filename = FileFormats.OTHER_VARIANTS_OUTPUT.format(
            identification=self.name
        )
        gvcf_output_filename = FileFormats.GVCF_OUTPUT.format(identification=self.name)

        vc_config = {
            self.name: {
                VariantCallingKeys.LIBRARY: self.library,
                VariantCallingKeys.PARAMS: self.params,
                VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING,
                    self.library,
                    unfiltered_variants_output_filename,
                ),
                VariantCallingKeys.ALL_VARIANTS_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING,
                    self.library,
                    all_variants_output_filename,
                ),
                VariantCallingKeys.SNP_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING, self.library, snp_output_filename
                ),
                VariantCallingKeys.INDEL_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING, self.library, indel_output_filename
                ),
                VariantCallingKeys.GVCF_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING, self.library, gvcf_output_filename
                ),
                VariantCallingKeys.OTHER_VARIANTS_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING,
                    self.library,
                    other_variants_output_filename,
                ),
                VariantCallingKeys.OUTPUT_DIR: join_paths(
                    OutputFolders.VARIANT_CALLING, self.library
                ),
                VariantCallingKeys.OUTPUT_TYPE: "VCF" if not self.gvcf else "GVCF",
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

        if self.bed_file is not None:
            vc_config[self.name][VariantCallingKeys.BED_FILE] = self.bed_file

        config = {self.key: vc_config}
        return config
