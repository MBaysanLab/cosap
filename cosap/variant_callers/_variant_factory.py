from subprocess import call
from ._mutect2_variantcaller import Mutect2VariantCaller
from ._octopus_variantcaller import OctopusVariantCaller
from ._somaticsniper_variantcaller import SomaticSniperVariantCaller
from ._variantcallers import _VariantCaller
from ._varscan_variantcaller import VarScanVariantCaller
from ._muse_variantcaller import MuseVariantCaller
from ._vardict_variantcaller import VarDictVariantCaller


class VariantCallerFactory:
    MUTECT2_CALLER = "mutect"
    SOMATICSNIPER_CALLER = "somaticsniper"
    STRELKA2_CALLER = "strelka"
    VARSCAN_CALLER = "varscan"
    OCTOPUS_CALLER = "octopus"
    MUSE_CALLER = "muse"
    VARDICT_CALLER = "vardict"

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
        elif caller_type == cls.MUSE_CALLER:
            caller = MuseVariantCaller
        elif caller_type == cls.VARDICT_CALLER:
            caller = VarDictVariantCaller
        else:
            raise Exception(f"Unknown caller type: {caller_type}")

        return caller
