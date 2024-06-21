from ._mosdepth import Mosdepth
from ._qualimap import Qualimap
from ._quality_controllers import _QualityController


class QualityControllerFactory:
    QUALIMAP = "qualimap"
    MOSDEPTH = "mosdepth"

    @classmethod
    def create(cls, quality_controller_type=str) -> _QualityController:
        quality_controller_type = str(quality_controller_type).lower()

        if quality_controller_type == cls.QUALIMAP:
            return Qualimap
        elif quality_controller_type == cls.MOSDEPTH:
            return Mosdepth
        else:
            raise Exception("Unknown quality controller")
