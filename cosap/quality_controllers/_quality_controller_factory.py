from re import I

from ._qualimap import Qualimap
from ._quality_controllers import _QualityController


class QualityContollerFactory:
    QUALIMAP = "qualimap"

    @classmethod
    def create(cls, quality_controller_type=str) -> _QualityController:
        quality_controller_type = str(quality_controller_type).lower()

        if quality_controller_type == cls.QUALIMAP:
            return Qualimap
        else:
            raise Exception("Unknown quality controller")
