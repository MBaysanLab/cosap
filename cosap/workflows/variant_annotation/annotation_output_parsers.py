import os
from logging import getLogger

import numpy as np
import pandas as pd
import io

from ._config import VariantAnnotationConfig

logger = getLogger(__name__)


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


class AnnovarResultsProcessor:
    """Dedicated class for processing Annovar results."""

    def __init__(self, config: VariantAnnotationConfig):
        self.config = config

    def process(self, path: str, subtype: str) -> pd.DataFrame:
        """Main processing method with comprehensive error handling."""
        try:
            self._validate_subtype(subtype)

            if not os.path.exists(path):
                raise FileNotFoundError(f"Annovar results file not found: {path}")

            df = self._load_and_validate_data(path)
            if df.empty:
                logger.warning(f"No data found in Annovar results file: {path}")
                return self._create_empty_dataframe(subtype)

            return self._apply_transformations(df, subtype)

        except Exception as e:
            logger.error(f"Error processing Annovar results from {path}: {e}")
            return self._create_empty_dataframe(subtype)

    def _validate_subtype(self, subtype: str) -> None:
        """Validate that the subtype is supported."""
        valid_subtypes = ["cancervar", "intervar", "annovar"]
        if subtype.lower() not in valid_subtypes:
            raise ValueError(f"subtype must be one of {valid_subtypes}")

    def _load_and_validate_data(self, path: str) -> pd.DataFrame:
        """Load and validate Annovar data."""
        df = pd.read_csv(path, sep="\t", dtype={"#Uploaded_variation": str})

        # Validate required columns exist
        required_cols = ["#Chr", "Start", "Ref", "Alt"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(
                f"Missing required columns in Annovar output: {missing_cols}"
            )

        return df

    def _apply_transformations(self, df: pd.DataFrame, subtype: str) -> pd.DataFrame:
        """Apply all data transformations."""
        df = self._process_variant_identifiers(df)
        df = self._apply_column_mapping(df, subtype)
        df = self._select_output_columns(df, subtype)
        return df

    def _process_variant_identifiers(self, df: pd.DataFrame) -> pd.DataFrame:
        """Create variant ID from chromosome, position, ref, and alt."""
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
        return df

    def _apply_column_mapping(self, df: pd.DataFrame, subtype: str) -> pd.DataFrame:
        """Apply column mapping based on subtype."""
        column_mapping = self._get_column_mapping(subtype)
        return df.rename(columns={v: k for k, v in column_mapping.items()})

    def _get_column_mapping(self, subtype: str) -> dict:
        """Get column mapping configuration for the specified subtype."""
        column_mapping = {
            f"{subtype}_classification": "Classification",
            "evidence": "Evidence",
            "interpro_domain": "Interpro_domain",
            "omim": "OMIM",
            "orpha_info": "Orpha",
        }

        if subtype == "cancervar":
            column_mapping["cosmic"] = "cosmic91"

        return column_mapping

    def _select_output_columns(self, df: pd.DataFrame, subtype: str) -> pd.DataFrame:
        """Select only the relevant output columns."""
        column_mapping = self._get_column_mapping(subtype)
        output_columns = ["variant_id", *column_mapping.keys()]

        # Only select columns that exist in the dataframe
        available_columns = [col for col in output_columns if col in df.columns]
        return df[available_columns]

    def _create_empty_dataframe(self, subtype: str) -> pd.DataFrame:
        """Create empty DataFrame with expected columns."""
        column_mapping = self._get_column_mapping(subtype)
        columns = ["variant_id", *column_mapping.keys()]
        return pd.DataFrame(columns=columns)
