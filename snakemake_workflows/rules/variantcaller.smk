from typing import List, Dict
import os
from cosap.variant_callers._variant_factory import VariantCallerFactory
from cosap._formats import FileFormats
from cosap._pipeline_config import VariantCallingKeys,PipelineKeys

rule variant_caller:
    input:
        germline_bam = lambda wildcards: FileFormats.CALIBRATION_OUTPUT.format(identification=wildcards.germline_identification),
        tumor_bam = lambda wildcards: FileFormats.CALIBRATION_OUTPUT.format(identification=wildcards.tumor_identification)
    output:
        vcf = FileFormats.GATK_SNP_OUTPUT.format(germline_identification="{germline_identification}", tumor_identification="{tumor_identification}")
    run:
        variant_caller = VariantCallerFactory(
            caller_type=config[PipelineKeys.VARIANT_CALLING][f"{wildcards.germline_identification}_{wildcards.tumor_identification}"][VariantCallingKeys.LIBRARY]
            )
        variant_caller.call_variants(config[PipelineKeys.VARIANT_CALLING][f"{wildcards.germline_identification}_{wildcards.tumor_identification}"])