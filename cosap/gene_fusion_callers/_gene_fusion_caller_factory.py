from ._genefuse_fusion_caller import GeneFuse
from ._gene_fusion_callers  import _GeneFusionCaller


class GeneFusionCallerFactory:
    GENEFUSE: str = "genefuse"

    @classmethod
    def create(cls, caller_type: str) -> _GeneFusionCaller:
        caller_type = str(caller_type).lower()

        if caller_type == cls.GENEFUSE:
            caller = GeneFuse
        else:
            raise Exception(f"Unknown mapper type: {caller_type}")

        return caller
