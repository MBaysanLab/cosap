import pandas as pd

from .._utils import convert_list_to_annovar_input, convert_list_to_ensembl_vep_input
from ..pipeline_builder.builders import Annotator
from ..tools.annotators import AnnotatorFactory


class VariantMultipleAnnotator:
    def __init__(self, variants: list):
        self.variants = variants

    def annotate(self):
        vep_input_file = convert_list_to_ensembl_vep_input(self.variants)
        annovar_input_file = convert_list_to_annovar_input(self.variants)

        vep_config = Annotator(library="vep", input_file=vep_input_file)
        intervar_config = Annotator(library="intervar", input_file=annovar_input_file)
        cancervar_config = Annotator(library="cancervar", input_file=annovar_input_file)

        vep_annotator = AnnotatorFactory.create_annotator("vep")
        intervar_annotator = AnnotatorFactory.create_annotator("intervar")
        cancervar_annotator = AnnotatorFactory.create_annotator("cancervar")

        vep_annotator.annotate(vep_config)
        intervar_annotator.annotate(intervar_config)
        cancervar_annotator.annotate(cancervar_config)

        vep_results = self._parse_ensembl_vep_results(vep_config.get_output())
        intervar_results = self._parse_annovar_results(
            intervar_config.get_output(), subtype="intervar"
        )
        cancervar_results = self._parse_annovar_results(
            cancervar_config.get_output(), subtype="cancervar"
        )

        results = pd.merge(
            vep_results,
            intervar_results,
            on="variant_id",
            how="outer",
        )

        return results.to_dict(orient="records")

    def _parse_ensembl_vep_results(self, path: str) -> pd.DataFrame:
        """
        Parses the output of ensembl vep and returns a list of variants.
        """

        vep_column_mapping = {
            "variant_id": "#Uploaded_variation",
            "ref": "",
            "alt": "Allele",
            "gene_symbol": "Gene",
            "feature": "Feature",
            "feature_type": "Feature_type",
            "consequence": "Consequence",
            "cdna_position": "cDNA_position",
            "rs_id": "Existing_variation",
            "impact": "IMPACT",
            "hgvs_c": "HGVSc",
            "gnomad_af": "AF",
            "aa_change": "Amino_acids",
            "sift_score": "SIFT",
            "polyphen_score": "PolyPhen",
            "clinical_significane": "CLIN_SIG",
        }

        df = pd.read_csv(path, sep="\t", dtype={"#Uploaded_variation": str})
        df = df.rename(columns={v: k for k, v in vep_column_mapping.items()})
        df = df[[*vep_column_mapping.keys()]]
        return df

    def _parse_annovar_results(self, path: str, subtype: str) -> pd.DataFrame:
        """
        Parses the output of annovar and returns a list of variants.
        """

        subtypes = ["cancervar", "intervar", "annovar"]
        if subtype.lower() not in subtypes:
            raise ValueError(f"subtype must be one of {subtypes}")

        annovar_column_mapping = {
            "gene": "Gene.ensGene",
            f"{subtype}_classification": "Classifiction",
            "evidence": "Evidence",
            "interpro_domain": "Interpro_domain",
            "omim": "OMIM",
            "orphanet_info": "Orpha",
            "cosmic": "cosmic91",
        }

        df = pd.read_csv(path, sep="\t", dtype={"#Uploaded_variation": str})
        df["variant_id"] = (
            "chr"
            + df["Chr"].astype(str)
            + "_"
            + df["Start"].astype(str)
            + "_"
            + df["Ref"].astype(str)
            + "/"
            + df["Alt"].astype(str)
        )
        df = df.rename(columns={v: k for k, v in annovar_column_mapping.items()})

        df = df[
            [annovar_column_mapping["variant_id"], *annovar_column_mapping.values()]
        ]

        return df
