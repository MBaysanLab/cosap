from ._cnv_callers import _CNVCaller
from ._cnvkit_cnv_caller import CNVKit
from ._gatk_germline_cnv_caller import GATKGermlineCNVCaller


class CNVCallerFactory:
    CNVKIT: str = "cnvkit"
    GATK_GERMLINE_CNV: str = "gatk_germline_cnv"

    @classmethod
    def create(cls, caller_type: str) -> _CNVCaller:
        caller_type = str(caller_type).lower()

        if caller_type == cls.CNVKIT:
            caller = CNVKit
        elif caller_type == cls.GATK_GERMLINE_CNV:
            caller = GATKGermlineCNVCaller
        else:
            raise Exception(f"Unknown copy number variation caller type: {caller_type}")

        return caller
