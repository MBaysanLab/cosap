import io
import os

import pandas as pd

from .._pipeline_config import PipelineKeys
from .._utils import (
    convert_list_to_annovar_input,
    convert_list_to_ensembl_vep_input,
    join_paths,
)
from ..pipeline_builder.builders import Annotator, VCFReader
from ..tools.annotators import AnnotatorFactory


class VariantMultipleAnnotator:
    def __init__(self, variants: list, workdir: str):
        self.variants = variants
        self.workdir = workdir

    def annotate(self):
        """
        Annotates the variants using VEP, Intervar, and Cancervar.

        Returns:
            dict: A dictionary containing the annotated variant records.
        """
        vep_input_file = convert_list_to_ensembl_vep_input(self.variants, self.workdir)
        annovar_input_file = convert_list_to_annovar_input(self.variants, self.workdir)

        vep_input = VCFReader(filename=vep_input_file)
        annovar_input = VCFReader(filename=annovar_input_file)

        vep_config = Annotator(
            library="vep",
            input_step=vep_input,
            sample_name="multiple_anotator",
            name="multiple_anotator",
            input_type="ensembl",
            output_dir=self.workdir,
        )
        intervar_config = Annotator(
            library="intervar",
            input_step=annovar_input,
            sample_name="multiple_anotator",
            input_type="AVinput",
            name="multiple_anotator",
            output_dir=self.workdir,
        )
        cancervar_config = Annotator(
            library="cancervar",
            input_step=annovar_input,
            sample_name="multiple_anotator",
            input_type="AVinput",
            name="multiple_anotator",
            output_dir=self.workdir,
        )

        vep_annotator = AnnotatorFactory.create("vep")
        intervar_annotator = AnnotatorFactory.create("intervar")
        cancervar_annotator = AnnotatorFactory.create("cancervar")

        if not os.path.exists(join_paths(self.workdir, vep_config.get_output())):
            vep_annotator.annotate(
                self._get_annotator_config(vep_config), workdir=self.workdir
            )

        if not os.path.exists(join_paths(self.workdir, intervar_config.get_output())):
            intervar_annotator.annotate(
                self._get_annotator_config(intervar_config), workdir=self.workdir
            )

        if not os.path.exists(join_paths(self.workdir, cancervar_config.get_output())):
            cancervar_annotator.annotate(
                self._get_annotator_config(cancervar_config), workdir=self.workdir
            )

        vep_results = self._parse_ensembl_vep_results(
            join_paths(self.workdir, vep_config.get_output())
        )
        intervar_results = self._parse_annovar_results(
            join_paths(self.workdir, intervar_config.get_output()), subtype="intervar"
        )
        cancervar_results = self._parse_annovar_results(
            join_paths(self.workdir, cancervar_config.get_output()), subtype="cancervar"
        )

        results = pd.merge(
            vep_results,
            intervar_results,
            on="variant_id",
            how="outer",
        ).merge(
            cancervar_results,
            on="variant_id",
            how="outer",
            suffixes=("_drop", ""),
        )

        results = results.drop(
            columns=[c for c in results.columns if c.endswith("_drop")]
        )

        # Make all columns lowercase
        results.columns = [c.lower() for c in results.columns]
        results = self._convert_null_values_to_nan(results)

        # Remove temp input files
        os.remove(vep_input_file)
        os.remove(annovar_input_file)

        return results.to_dict(orient="records")

    def _parse_ensembl_vep_results(self, path: str) -> pd.DataFrame:
        """
        Parses the output of ensembl vep and returns a list of variants.
        """

        vep_column_mapping = {
            "variant_id": "#Uploaded_variation",
            "ref": "REF_ALLELE",
            "alt": "Allele",
            "gene_id": "Gene",
            "gene_symbol": "SYMBOL",
            "feature": "Feature",
            "feature_type": "Feature_type",
            "consequence": "Consequence",
            "cdna_position": "cDNA_position",
            "rs_id": "Existing_variation",
            "impact": "IMPACT",
            "hgvs_g": "HGVSg",
            "hgvs_c": "HGVSc",
            "hgvs_p": "HGVSp",
            "mane_select": "MANE_SELECT",
            "gnomad_af": "gnomADg_AF",
            "aa_change": "Amino_acids",
            "sift_score": "SIFT",
            "polyphen_score": "PolyPhen",
            "clinical_significane": "CLIN_SIG",
        }

        # Skip rows that start with ##
        df = self._read_ensembl_vcf_to_dataframe(path)
        df = df.rename(columns={v: k for k, v in vep_column_mapping.items()})
        # df = df[[*vep_column_mapping.keys()]]ü

        # Convert chr_pos_ref/alt to chr_pos_ref_alt
        df["variant_id"] = df["variant_id"].str.replace("/", "_")

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
            f"{subtype}_classification": "Classification",
            "evidence": "Evidence",
            "interpro_domain": "Interpro_domain",
            "omim": "OMIM",
            "orphanet_info": "Orpha",
        }
        if subtype == "cancervar":
            annovar_column_mapping["cosmic"] = "cosmic91"

        df = pd.read_csv(path, sep="\t", dtype={"#Uploaded_variation": str})
        df["variant_id"] = (
            "chr"
            + df["#Chr"].astype(str)
            + "_"
            + df["Start"].astype(str)
            + "_"
            + df["Ref"].astype(str)
            + "_"
            + df["Alt"].astype(str)
        )
        df = df.rename(columns={v: k for k, v in annovar_column_mapping.items()})

        df = df[["variant_id", *annovar_column_mapping.keys()]]

        return df

    def _get_annotator_config(self, Annotator) -> dict:
        annotator_name = Annotator.name
        annotator_config = Annotator.get_config()

        return annotator_config[PipelineKeys.ANNOTATION][annotator_name]

    def _read_ensembl_tabs_to_dataframe(self, path: str) -> pd.DataFrame:
        """
        Reads the output of ensembl vep and returns a dataframe.
        """
        with open(path, "r") as f:
            lines = [l for l in f if not l.startswith("##")]
            df = pd.read_csv(
                io.StringIO("".join(lines)),
                sep="\t",
            )

        df.reset_index(inplace=True)

        return df

    def _read_ensembl_vcf_to_dataframe(self, path: str) -> pd.DataFrame:
        """
        Reads the output of ensembl vep and returns a dataframe.
        """
        with open(path, "r") as f:
            lines = list(f)
            header_lines = [l for l in lines if l.startswith("##")]
            data_lines = [l for l in lines if not l.startswith("##")]
            df = pd.read_csv(
                io.StringIO("".join(data_lines)),
                sep="\t",
            )

            # Get vep fields from vcf header and add them to the dataframe
            vep_info_field = [
                l for l in header_lines if l.startswith("##INFO=<ID=CSQ")
            ][0]
            vep_fields = vep_info_field.split("Format: ")[1].split("|")

            # The vep data is the last value of the INFO field
            df_vep = (
                df["INFO"].str.split("CSQ=", expand=True)[1].str.split("|", expand=True)
            )
            df_vep.columns = vep_fields

            df_vep["variant_id"] = (
                df["#CHROM"]
                + "_"
                + df["POS"].astype(str)
                + "_"
                + df["REF"]
                + "_"
                + df["ALT"]
            )
            df_vep["ref"] = df["REF"]
            df_vep["alt"] = df["ALT"]

        df_vep.reset_index(inplace=True)

        return df_vep

    def _convert_null_values_to_nan(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Converts all None values to NaN in a dataframe.

        None values can be "-", ".", or any other string that represents a null value.
        """
        return df.replace("-", pd.NA).replace(".", pd.NA).replace("", pd.NA)
