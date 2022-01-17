from cosap.variant_callers._variant_factory import VariantCallerFactory
from cosap._pipeline_config import VariantCallingKeys, PipelineKeys

config = snakemake.config
wildcards = snakemake.wildcards

variant_caller = VariantCallerFactory.create(
    caller_type=config[PipelineKeys.VARIANT_CALLING][wildcards.identification][
        VariantCallingKeys.LIBRARY]
        )
variant_caller.call_variants(
    config[PipelineKeys.VARIANT_CALLING][wildcards.identification]
)