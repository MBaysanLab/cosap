from ._gatk_mitochondrial_caller import GATKMitochondrialCaller


class MitochondrialCallerFactory:
    GATK_MITOCHONDRIAL = "gatk_mitochondrial"

    @classmethod
    def create(cls, caller_type: str):
        if str(caller_type).lower() == cls.GATK_MITOCHONDRIAL:
            return GATKMitochondrialCaller
        raise ValueError(f"Unknown mitochondrial caller type: {caller_type}")
