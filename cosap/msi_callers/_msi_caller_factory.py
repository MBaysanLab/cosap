from ._msisensorpro_msicaller import MSISensorPro
from ._msi_callers  import _MSICaller


class MSICallerFactory:
    GENEFUSE: str = "msisensor"

    @classmethod
    def create(cls, caller_type: str) -> _MSICaller:
        caller_type = str(caller_type).lower()

        if caller_type == cls.GENEFUSE:
            caller = MSISensorPro
        else:
            raise Exception(f"Unknown microsatellite instability caller type: {caller_type}")

        return caller
