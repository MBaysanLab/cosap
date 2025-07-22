from ._annotators import _Annotator
from ._annotsv_annotator import AnnotSVAnnotator
from ._annovar_annotator import AnnovarAnnotator
from ._cancervar_annotatator import CancervarAnnotator
from ._ensembl_vep_annotator import VepAnnotator
from ._intervar_annotator import IntervarAnnotator
from ._pharmcat_annotator import PharmcatAnnotator


class AnnotatorFactory:
    ANNOVAR_ANNOTATOR = "annovar"
    INTERVAR_ANNOTATOR = "intervar"
    CANCERVAR_ANNOTATOR = "cancervar"
    VEP_ANNOTATOR = "vep"
    PHARMGKB_ANNOTATOR = "pharmcat"
    ANNOTSV = "annotsv"

    @classmethod
    def create(cls, annotator_type: str) -> _Annotator:
        annotator_type = str(annotator_type).lower()

        if annotator_type == cls.ANNOVAR_ANNOTATOR:
            annotator = AnnovarAnnotator
        elif annotator_type == cls.INTERVAR_ANNOTATOR:
            annotator = IntervarAnnotator
        elif annotator_type == cls.VEP_ANNOTATOR:
            annotator = VepAnnotator
        elif annotator_type == cls.PHARMGKB_ANNOTATOR:
            annotator = PharmcatAnnotator
        elif annotator_type == cls.CANCERVAR_ANNOTATOR:
            annotator = CancervarAnnotator
        elif annotator_type == cls.ANNOTSV:
            annotator = AnnotSVAnnotator
        else:
            raise Exception(f"Unknown annotator type: {annotator_type}")

        return annotator
