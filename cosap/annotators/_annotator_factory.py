from ._annotators import _Annotator
from ._annovar_annotator import AnnovarAnnotator
from ._ensembl_vep_annotator import VepAnnotator
from ._intervar_annotator import IntervarAnnotator
from ._pharmcat_annotator import PharmcatAnnotator


class AnnotatorFactory:
    ANNOVAR_ANNOTATOR = "annovar"
    INTERVAR_ANNOTATOR = "intervar"
    VEP_ANNOTATOR = "vep"
    PHARMGKB_ANNOTATOR = "pharmgkb"

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
        else:
            raise Exception(f"Unknown annotator type: {annotator_tpye}")

        return annotator
