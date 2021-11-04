from ._mutect2_variantcaller import Mutect2VariantCaller
from ._somaticsniper_variantcaller import SomaticSniperVariantCaller
from ._variantcallers import _VariantCaller
from ._varscan_variantcaller import VarScanVariantCaller
from ._octopus_variantcaller import OctopusVariantCaller


class VariantCallerFactory:
    MUTECT2_CALLER = "mutect"
    SOMATICSNIPER_CALLER = "somaticsniper"
    STRELKA2_CALLER = "strelka"
    VARSCAN_CALLER = "varscan"
    OCTOPUS_CALLER = "octopus"

    @classmethod
    def create(cls, caller_type: str) -> _VariantCaller:
        caller_type = str(caller_type).lower()

        if caller_type == cls.MUTECT2_CALLER:
            caller = Mutect2VariantCaller
        elif caller_type == cls.SOMATICSNIPER_CALLER:
            caller = SomaticSniperVariantCaller
        elif caller_type == cls.VARSCAN_CALLER:
            caller = VarScanVariantCaller
        elif caller_type == cls.OCTOPUS_CALLER:
            caller = OctopusVariantCaller
        else:
            raise Exception(f"Unknown caller type: {caller_type}")

        return caller
