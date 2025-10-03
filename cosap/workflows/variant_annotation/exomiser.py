import os

import pandas as pd

from ..._pipeline_config import PipelineKeys
from ...pipeline_builder.builders import Annotator, VCFReader
from ...tools.annotators import AnnotatorFactory


class ExomiserPrioritizer:
    def __init__(self, vcf_path: str, hpo_list: list, workdir: str = None):
        self.vcf_path = vcf_path
        self.hpo_list = hpo_list
        self.workdir = workdir
        self.output_file = None

    def prioritize(self):
        """
        Runs Exomiser to prioritize variants based on the provided VCF and HPO terms.

        Returns:
            str: Path to the Exomiser output file.
        """
        vcf = VCFReader(self.vcf_path)
        sample_name = vcf.get_sample_name()

        exomiser_builder = Annotator(
            library="exomiser",
            input_step=vcf,
            sample_name=sample_name,
            name=sample_name,
            input_type="vcf",
            output_dir=self.workdir,
            params={
                "hpo_terms": self.hpo_list,
            },
        )
        exomiser_config = exomiser_builder.get_config()[PipelineKeys.ANNOTATION][
            sample_name
        ]

        exomiser_annotator = AnnotatorFactory.create("exomiser")
        exomiser_annotator.annotate(exomiser_config)

        self.output_file = os.path.join(
            self.workdir, f"annotated_{sample_name}.variants.tsv"
        )

        variants_json = self._parse_exomiser_results(self.output_file)

        return variants_json

    def _parse_exomiser_results(self, path: str) -> pd.DataFrame:
        """
        Parses the output of Exomiser and returns a dataframe.
        """
        df = pd.read_csv(path, sep="\t")
        df.rename(columns={"#RANK": "RANK", "P-VALUE": "P_VALUE"}, inplace=True)
        df["variant_id"] = (
            df[["CONTIG", "START", "REF", "ALT"]].astype(str).agg("_".join, axis=1)
        )

        # Lowercase all columns
        df.columns = [c.lower() for c in df.columns]

        # Only need RANK, P_VALUE and VARIANT_ID
        df = df[["rank", "p_value", "variant_id"]]

        return df.to_json(orient="records")
