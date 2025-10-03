import io
import os
from typing import Optional, Dict, List

import pandas as pd
import numpy as np

from .._pipeline_config import PipelineKeys
from .._utils import (
    convert_list_to_annovar_input,
    convert_list_to_vcf_input,
    join_paths,
)
from ..pipeline_builder.builders import Annotator, VCFReader
from ..tools.annotators import AnnotatorFactory
from logging import getLogger

logger = getLogger(__name__)


class VariantAnnotationConfig:
    """Configuration class for variant annotation mappings and thresholds."""

    VEP_COLUMN_MAPPING = {
        "variant_id": "Uploaded_variation",
        "chr": "#CHROM",
        "pos": "POS",
        "ref": "REF",
        "alt": "ALT",
        "gene_id": "Gene",
        "gene_symbol": "SYMBOL",
        "feature": "Feature",
        "feature_type": "Feature_type",
        "consequence": "Consequence",
        "cdna_position": "cDNA_position",
        "cds_position": "CDS_position",
        "rs_id": "Existing_variation",
        "impact": "IMPACT",
        "hgvsg": "HGVSg",
        "hgvsc": "HGVSc",
        "hgvsp": "HGVSp",
        "vep_pick": "PICK",
        "mane": "MANE",
        "canonical": "CANONICAL",
        "gnomadg_af": "gnomADg_AF",
        "gnomadg_afr_af": "gnomADg_AFR_AF",
        "gnomadg_ami_af": "gnomADg_AMI_AF",
        "gnomadg_amr_af": "gnomADg_AMR_AF",
        "gnomadg_asj_af": "gnomADg_ASJ_AF",
        "gnomadg_eas_af": "gnomADg_EAS_AF",
        "gnomadg_fin_af": "gnomADg_FIN_AF",
        "gnomadg_mid_af": "gnomADg_MID_AF",
        "gnomadg_nfe_af": "gnomADg_NFE_AF",
        "gnomadg_oth_af": "gnomADg_OTH_AF",
        "gnomadg_sas_af": "gnomADg_SAS_AF",
        "aminoacid_change": "Amino_acids",
        "sift_score": "SIFT",
        "polyphen_score": "PolyPhen",
        "clinical_significance": "CLIN_SIG",
        "alphamissense_class": "am_class",
        "alphamissense_pathogenicity": "am_pathogenicity",
        "revel_score": "REVEL",
        "spliceai_pred_dp_ag": "SpliceAI_pred_DP_AG",
        "spliceai_pred_dp_al": "SpliceAI_pred_DP_AL",
        "spliceai_pred_dp_dg": "SpliceAI_pred_DP_DG",
        "spliceai_pred_dp_dl": "SpliceAI_pred_DP_DL",
        "spliceai_pred_ds_ag": "SpliceAI_pred_DS_AG",
        "spliceai_pred_ds_al": "SpliceAI_pred_DS_AL",
        "spliceai_pred_ds_dg": "SpliceAI_pred_DS_DG",
        "spliceai_pred_ds_dl": "SpliceAI_pred_DS_DL",
        "spliceai_pred_symbol": "SpliceAI_pred_SYMBOL",
        "pubmed": "PUBMED",
    }

    SPLICEAI_SCORE_COLUMNS = [
        "spliceai_pred_ds_ag",
        "spliceai_pred_ds_al",
        "spliceai_pred_ds_dg",
        "spliceai_pred_ds_dl",
    ]

    SPLICEAI_THRESHOLDS = {
        "high_impact": 0.8,
        "moderate_impact": 0.5,
        "low_impact": 0.2,
    }

    REVEL_THRESHOLDS = {"pathogenic": 0.75, "likely_pathogenic": 0.5, "uncertain": 0.25}


class VariantMultipleAnnotator:
    def __init__(self, variants: list, workdir: str):
        self.variants = variants
        self.workdir = workdir
        self.vep_input_file = None
        self.vep_output_file = None
        self.annovar_input_file = None
        self.intervar_output_file = None
        self.cancervar_output_file = None

        # Initialize processors
        self.config = VariantAnnotationConfig()
        self.vep_processor = VEPResultsProcessor(self.config)

    def annotate(self):
        """
        Annotates the variants using VEP, Intervar, and Cancervar.

        Returns:
            dict: A dictionary containing the annotated variant records.
        """
        vep_input_file = convert_list_to_vcf_input(self.variants, self.workdir)
        annovar_input_file = convert_list_to_annovar_input(self.variants, self.workdir)

        self.vep_input_file = vep_input_file
        self.annovar_input_file = annovar_input_file

        vep_input = VCFReader(filename=vep_input_file)
        annovar_input = VCFReader(filename=annovar_input_file)

        vep_config = Annotator(
            library="vep",
            input_step=vep_input,
            sample_name="multiple_anotator",
            name="multiple_anotator_vep",
            input_type="vcf",
            output_dir=self.workdir,
        )
        intervar_config = Annotator(
            library="intervar",
            input_step=annovar_input,
            sample_name="multiple_anotator",
            input_type="AVinput",
            name="multiple_anotator_intervar",
            output_dir=self.workdir,
        )
        cancervar_config = Annotator(
            library="cancervar",
            input_step=annovar_input,
            sample_name="multiple_anotator",
            input_type="AVinput",
            name="multiple_anotator_cancervar",
            output_dir=self.workdir,
        )

        vep_annotator = AnnotatorFactory.create("vep")
        intervar_annotator = AnnotatorFactory.create("intervar")
        cancervar_annotator = AnnotatorFactory.create("cancervar")

        if not os.path.exists(join_paths(self.workdir, vep_config.get_output())):
            logger.info(
                f"Annotating {len(self.variants)} variants with VEP using input file: {self.vep_input_file}"
            )
            vep_annotator.annotate(self._get_annotator_config(vep_config))

        if not os.path.exists(join_paths(self.workdir, intervar_config.get_output())):
            logger.info(
                f"Annotating {len(self.variants)} variants with Intervar using input file: {self.annovar_input_file}"
            )
            intervar_annotator.annotate(self._get_annotator_config(intervar_config))

        if not os.path.exists(join_paths(self.workdir, cancervar_config.get_output())):
            logger.info(
                f"Annotating {len(self.variants)} variants with Cancervar using input file: {self.annovar_input_file}"
            )
            cancervar_annotator.annotate(self._get_annotator_config(cancervar_config))

        self.vep_output_file = join_paths(self.workdir, vep_config.get_output())
        self.intervar_output_file = join_paths(
            self.workdir, intervar_config.get_output()
        )
        self.cancervar_output_file = join_paths(
            self.workdir, cancervar_config.get_output()
        )

    def get_results(self):
        vep_results = self._parse_ensembl_vep_results(self.vep_output_file)
        intervar_results = self._parse_annovar_results(
            self.intervar_output_file, subtype="intervar"
        )
        cancervar_results = self._parse_annovar_results(
            self.cancervar_output_file, subtype="cancervar"
        )
        if any([vep_results.empty, intervar_results.empty, cancervar_results.empty]):
            logger.warning(
                "One or more annotation results are empty. Returning empty DataFrame."
            )
            return pd.DataFrame(
                columns=vep_results.columns
                | intervar_results.columns
                | cancervar_results.columns
            ).to_json(orient="records")

        results = pd.merge(
            intervar_results,
            cancervar_results,
            on="variant_id",
            how="inner",
            suffixes=("_intervar", "_cancervar"),
        ).merge(
            vep_results,
            on="variant_id",
            how="right",
            suffixes=("", "_vep"),
        )

        # Merge orpha, omim and interpro columns from intervar and cancervar
        results["orpha_info"] = results["orpha_info_intervar"].combine_first(
            results["orpha_info_cancervar"]
        )
        results["omim"] = results["omim_intervar"].combine_first(
            results["omim_cancervar"]
        )
        results["interpro_domain"] = results["interpro_domain_intervar"].combine_first(
            results["interpro_domain_cancervar"]
        )

        # Drop the columns that are not needed
        results.drop(
            columns=[
                "orpha_info_intervar",
                "orpha_info_cancervar",
                "omim_intervar",
                "omim_cancervar",
                "interpro_domain_intervar",
                "interpro_domain_cancervar",
            ],
            inplace=True,
        )

        # Make all columns lowercase
        results.columns = [c.lower() for c in results.columns]
        results = self._convert_null_values_to_na(results)

        # # Remove temp input files
        # os.remove(self.vep_input_file)
        # os.remove(self.annovar_input_file)

        return results.to_json(orient="records")

    def _parse_ensembl_vep_results(self, path: str) -> pd.DataFrame:
        """Parses the output of ensembl vep and returns a list of variants."""
        return self.vep_processor.process(path)

    def _parse_annovar_results(self, path: str, subtype: str) -> pd.DataFrame:
        """
        Parses the output of annovar and returns a list of variants.
        """

        subtypes = ["cancervar", "intervar", "annovar"]
        if subtype.lower() not in subtypes:
            raise ValueError(f"subtype must be one of {subtypes}")

        annovar_column_mapping = {
            f"{subtype}_classification": "Classification",
            "evidence": "Evidence",
            "interpro_domain": "Interpro_domain",
            "omim": "OMIM",
            "orpha_info": "Orpha",
        }
        if subtype == "cancervar":
            annovar_column_mapping["cosmic"] = "cosmic91"

        try:
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

        except Exception as e:
            logger.error(f"Error parsing {subtype} results: {e}")

            # Return an empty DataFrame if parsing fails
            return pd.DataFrame(columns=["variant_id", *annovar_column_mapping.keys()])

        return df

    def _get_annotator_config(self, Annotator) -> dict:
        annotator_name = Annotator.name
        annotator_config = Annotator.get_config()

        return annotator_config[PipelineKeys.ANNOTATION][annotator_name]

    def _convert_null_values_to_na(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Converts all None values to NaN in a dataframe.

        None values can be "-", ".", or any other string that represents a null value.
        """
        return (
            df.replace([
                "None", "none", "NA", "N/A", "na", ".", "-", "", np.nan
            ], pd.NA)).fillna(pd.NA)


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


class VEPResultsProcessor:
    """Dedicated class for processing VEP results."""

    def __init__(self, config: VariantAnnotationConfig):
        self.config = config

    def process(self, path: str) -> pd.DataFrame:
        """Main processing method with comprehensive error handling."""
        try:
            if not os.path.exists(path):
                raise FileNotFoundError(f"VEP results file not found: {path}")

            df = self._load_and_validate_data(path)
            if df.empty:
                logger.warning(f"No data found in VEP results file: {path}")
                return self._create_empty_dataframe()

            return self._apply_transformations(df)

        except Exception as e:
            logger.error(f"Error processing VEP results from {path}: {e}")
            return self._create_empty_dataframe()

    def _load_and_validate_data(self, path: str) -> pd.DataFrame:
        """Load and validate VEP data."""
        df = self._read_ensembl_vcf_to_dataframe(path)

        # Validate required columns exist
        required_cols = ["#CHROM", "POS", "REF", "ALT"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in VEP output: {missing_cols}")

        return df

    def _apply_transformations(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply all data transformations."""
        df = self._load_and_map_vep_data(df)
        df = self._process_variant_identifiers(df)
        df = self._process_boolean_columns(df)
        df = self._process_consequence_column(df)
        df = self._process_prediction_scores(df)
        return df

    def _load_and_map_vep_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply column mapping to VEP data."""
        df = df.rename(
            columns={v: k for k, v in self.config.VEP_COLUMN_MAPPING.items()}
        )
        return df[
            [col for col in self.config.VEP_COLUMN_MAPPING.keys() if col in df.columns]
        ]

    def _process_variant_identifiers(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process variant ID formatting."""
        df["variant_id"] = df["variant_id"].str.replace("/", "_")
        df["variant_id"] = np.where(
            df["variant_id"].str.startswith("chr"),
            df["variant_id"],
            "chr" + df["variant_id"],
        )
        return df

    def _process_boolean_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert string columns to boolean."""
        if "canonical" in df.columns:
            df["canonical"] = df["canonical"].str.lower() == "yes"
        if "vep_pick" in df.columns:
            df["vep_pick"] = df["vep_pick"] == "1"
        return df

    def _process_consequence_column(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process consequence column to get most severe consequence."""
        if "consequence" in df.columns:
            df["consequence"] = df["consequence"].str.split("&").str[0]
        return df

    def _process_prediction_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process all prediction scores (SpliceAI, REVEL, etc.)."""
        df = self._process_spliceai_scores(df)
        df = self._process_revel_scores(df)
        return df

    def _process_spliceai_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process SpliceAI prediction scores."""
        score_columns = self.config.SPLICEAI_SCORE_COLUMNS

        # Convert to numeric
        for col in score_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        # Calculate max score and interpretation
        available_cols = [col for col in score_columns if col in df.columns]
        if available_cols:
            df["spliceai_max_score"] = df[available_cols].max(axis=1)
            df["spliceai_interpretation"] = self._interpret_spliceai_score(
                df["spliceai_max_score"]
            )

        return df

    def _process_revel_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process REVEL prediction scores."""
        if "revel_score" in df.columns:
            df["revel_score"] = pd.to_numeric(df["revel_score"], errors="coerce")
            df["revel_interpretation"] = self._interpret_revel_score(df["revel_score"])
        return df

    def _interpret_spliceai_score(self, scores: pd.Series) -> pd.Series:
        """Interpret SpliceAI scores based on thresholds."""
        thresholds = self.config.SPLICEAI_THRESHOLDS
        return np.where(
            scores >= thresholds["high_impact"],
            "high_impact",
            np.where(
                scores >= thresholds["moderate_impact"],
                "moderate_impact",
                np.where(
                    scores >= thresholds["low_impact"], "low_impact", "minimal_impact"
                ),
            ),
        )

    def _interpret_revel_score(self, scores: pd.Series) -> pd.Series:
        """Interpret REVEL scores based on thresholds."""
        thresholds = self.config.REVEL_THRESHOLDS
        return np.where(
            scores >= thresholds["pathogenic"],
            "pathogenic",
            np.where(
                scores >= thresholds["likely_pathogenic"],
                "likely_pathogenic",
                np.where(scores >= thresholds["uncertain"], "uncertain", "benign"),
            ),
        )

    def _create_empty_dataframe(self) -> pd.DataFrame:
        """Create empty DataFrame with expected columns."""
        return pd.DataFrame(columns=self.config.VEP_COLUMN_MAPPING.keys())

    def _read_ensembl_vcf_to_dataframe(self, path: str) -> pd.DataFrame:
        """
        Reads the output of ensembl vep and returns a dataframe.
        """
        with open(path, "r") as f:
            lines = list(f)

        header_lines = [l for l in lines if l.startswith("##")]
        data_lines = [l for l in lines if not l.startswith("##")]

        if not data_lines:
            raise ValueError("No data lines found in the vep VCF file.")

        df = pd.read_csv(
            io.StringIO("".join(data_lines)),
            sep="\t",
        )

        # Get vep fields from vcf header and add them to the dataframe
        vep_info_field = [l for l in header_lines if l.startswith("##INFO=<ID=CSQ")][0]
        vep_fields = vep_info_field.split("Format: ")[1].rstrip("'>\n").split("|")

        # The vep data is the last value of the INFO field
        df["INFO"] = df["INFO"].str.extract(r"CSQ=([^;]+)")[0].str.split(",")
        df = df.explode("INFO")
        vep_columns = df["INFO"].str.split("|", expand=True)
        vep_columns.columns = vep_fields
        df_vep = pd.concat([df, vep_columns], axis=1)

        df_vep["variant_id"] = (
            df["#CHROM"]
            + "_"
            + df["POS"].astype(str)
            + "_"
            + df["REF"]
            + "_"
            + df["ALT"]
        )
        df_vep.reset_index(inplace=True, drop=True)

        return df_vep

    def _read_ensemb_default_output_to_dataframe(self, vep_file: str) -> pd.DataFrame:
        """
        Parses the default Ensembl VEP output into a Pandas DataFrame.

        Parameters:
        vep_file (str): Path to the VEP output file.

        Returns:
        pd.DataFrame: A DataFrame containing the parsed VEP output.
        """
        # Read the VEP output file into a DataFrame
        df = pd.read_csv(vep_file, sep="\t", comment="#", header=None)

        # Define the column names based on the VEP output format
        columns = [
            "Uploaded_variation",
            "Location",
            "Allele",
            "Gene",
            "Feature",
            "Feature_type",
            "Consequence",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "Amino_acids",
            "Codons",
            "Existing_variation",
            "Extra",
        ]

        # Assign the column names to the DataFrame
        df.columns = columns

        # Parse the 'Extra' column into separate columns
        extra_data = df["Extra"].str.split(";", expand=True)

        # Create a dictionary to store the extra information
        extra_dict = {}
        for col in extra_data.columns:
            # Split each key-value pair in the column
            key_value = extra_data[col].str.split("=", expand=True)
            if not key_value.empty:
                # Extract the key and values
                keys = key_value[0].dropna().unique()
                for key in keys:
                    if key:  # Ensure the key is not empty
                        # Extract values for the current key
                        values = key_value[key_value[0] == key][1]
                        extra_dict[key] = values

        # Add the extra information to the DataFrame
        for key, values in extra_dict.items():
            df[key] = values

        # Drop the original 'Extra' column
        df = df.drop(columns=["Extra"])

        return df
