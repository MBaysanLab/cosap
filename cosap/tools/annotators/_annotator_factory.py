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
    PHARMGKB_ANNOTATOR = "pharmgkb"
    ANNOTSV = "annotsv"

    @classmethod
    def create(cls, annotator_tpye: str) -> _Annotator:
        annotator_tpye = str(annotator_tpye).lower()

        if annotator_tpye == cls.ANNOVAR_ANNOTATOR:
            annotator = AnnovarAnnotator
        elif annotator_tpye == cls.INTERVAR_ANNOTATOR:
            annotator = IntervarAnnotator
        elif annotator_tpye == cls.VEP_ANNOTATOR:
            annotator = VepAnnotator
        elif annotator_tpye == cls.PHARMGKB_ANNOTATOR:
            annotator = PharmcatAnnotator
        elif annotator_tpye == cls.CANCERVAR_ANNOTATOR:
            annotator = CancervarAnnotator
        elif annotator_tpye == cls.ANNOTSV:
            annotator = AnnotSVAnnotator
        else:
            raise Exception(f"Unknown annotator type: {annotator_tpye}")

        return annotator
