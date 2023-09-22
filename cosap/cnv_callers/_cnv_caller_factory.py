from ._cnv_callers import _CNVCaller
from ._cnvkit_cnv_caller import CNVKit


class CNVCallerFactory:
    CNVKIT: str = "cnvkit"

    @classmethod
    def create(cls, caller_type: str) -> _CNVCaller:
        caller_type = str(caller_type).lower()

        if caller_type == cls.CNVKIT:
            caller = CNVKit
        else:
            raise Exception(f"Unknown copy number variation caller type: {caller_type}")

        return caller
