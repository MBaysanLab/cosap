from ._mutect2_variantcaller import Mutect2VariantCaller
from ._somaticsniper_variantcaller import SomaticSniperVariantCaller
from ._strelka_variantcaller import Strelka2VariantCaller
from ._varscan_variantcaller import VarScanVariantCaller
from ._variantcallers import _VariantCallable


class VariantCallerFactory:
    MUTECT2_CALLER = "mutect2" 
    SOMATICSNIPER_CALLER = "somaticsniper"
    STRELKA2_CALLER = "strelka2"
    VARSCAN_CALLER = "varscan"

    @classmethod
    def create(cls, caller_type: str) -> _VariantCallable:
        caller_type = str(caller_type).lower()

        if caller_type == cls.MUTECT2_CALLER:
            caller = Mutect2VariantCaller
        elif caller_type == cls.SOMATICSNIPER_CALLER:
            caller = SomaticSniperVariantCaller
        elif caller_type == cls.STRELKA2_CALLER:
            caller = Strelka2VariantCaller
        elif caller_type == cls.VARSCAN_CALLER:
            caller = VarScanVariantCaller
        else:
            raise Exception(f"Unknown caller type: {caller_type}")

        return caller
