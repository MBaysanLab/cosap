from subprocess import call

from ._deepvariant_variantcaller import DeepVariantVariantCaller
from ._haplotypecaller_variantcaller import HaplotypeCallerVariantCaller
from ._manta_variantcaller import MantaVariantCaller
from ._muse_variantcaller import MuseVariantCaller
from ._mutect2_variantcaller import Mutect2VariantCaller
from ._octopus_variantcaller import OctopusVariantCaller
from ._somaticsniper_variantcaller import SomaticSniperVariantCaller
from ._strelka_variantcaller import StrelkaVariantCaller
from ._vardict_variantcaller import VarDictVariantCaller
from ._variantcallers import _VariantCaller
from ._varnet_variantcaller import VarNetVariantCaller
from ._varscan_germline_variantcaller import VarScanGermlineVariantCaller
from ._varscan_variantcaller import VarScanVariantCaller


class VariantCallerFactory:
    MUTECT2_CALLER = "mutect"
    SOMATICSNIPER_CALLER = "somaticsniper"
    STRELKA2_CALLER = "strelka"
    VARSCAN_CALLER = "varscan"
    OCTOPUS_CALLER = "octopus"
    MUSE_CALLER = "muse"
    VARDICT_CALLER = "vardict"
    HAPLOTYPECALLER = "haplotypecaller"
    VARSCAN_GERMLINE_CALLER = "varscan_germline"
    VARNET = "varnet"
    MANTA = "manta"
    DEEPVARIANT = "deepvariant"

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
        elif caller_type == cls.STRELKA2_CALLER:
            caller = StrelkaVariantCaller
        elif caller_type == cls.HAPLOTYPECALLER:
            caller = HaplotypeCallerVariantCaller
        elif caller_type == cls.VARSCAN_GERMLINE_CALLER:
            caller = VarScanGermlineVariantCaller
        elif caller_type == cls.VARNET:
            caller = VarNetVariantCaller
        elif caller_type == cls.MANTA:
            caller = MantaVariantCaller
        elif caller_type == cls.DEEPVARIANT:
            caller = DeepVariantVariantCaller
        else:
            raise Exception(f"Unknown caller type: {caller_type}")

        return caller
