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
        df = self._process_clinical_significance(df)
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

    def _read_ensembl_vcf_from_string(self, vcf_string: str) -> pd.DataFrame:
        """
        Reads VCF content from a string and returns a dataframe.
        Reuses the logic from _read_ensembl_vcf_to_dataframe.

        Args:
            vcf_string: VCF content as string

        Returns:
            pd.DataFrame: Parsed VCF data
        """
        lines = vcf_string.split("\n")

        header_lines = [l for l in lines if l.startswith("##")]
        data_lines = [l for l in lines if not l.startswith("##") and l.strip()]

        if not data_lines:
            logger.warning("No data lines found in VCF content.")
            return pd.DataFrame()

        df = pd.read_csv(
            io.StringIO("\n".join(data_lines)),
            sep="\t",
        )

        # Get vep fields from vcf header and add them to the dataframe
        csq_lines = [l for l in header_lines if l.startswith("##INFO=<ID=CSQ")]
        if not csq_lines:
            logger.warning("No CSQ field found in VCF header.")
            return pd.DataFrame()

        vep_info_field = csq_lines[0]
        vep_fields = vep_info_field.split("Format: ")[1].rstrip('">\n').split("|")

        # The vep data is the last value of the INFO field
        df["INFO"] = df["INFO"].str.extract(r"CSQ=([^;]+)")[0].str.split(",")
        df = df.explode("INFO")

        if df["INFO"].isna().all():
            logger.warning("No CSQ data found in INFO field.")
            return pd.DataFrame()

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

    def _read_ensembl_vcf_to_dataframe(self, path: str) -> pd.DataFrame:
        """
        Reads the output of ensembl vep and returns a dataframe.
        """
        with open(path, "r", encoding="utf-8") as f:
            vcf_content = f.read()

        return self._read_ensembl_vcf_from_string(vcf_content)

    def process_record(self, record, header_lines: list = None) -> dict:
        """
        Process a single pysam.VariantFile record or a list of records with the same transformations as process().

        Args:
            record: pysam.VariantRecord object or list of pysam.VariantRecord objects
            header_lines: List of VCF header lines (optional, for efficiency)

        Returns:
            dict or list of dict: Processed record data, or None if processing fails
                                 If input is a list, returns a list of dicts
        """
        # Handle list of records
        if isinstance(record, list):
            return self._process_multiple_records(record, header_lines)

        # Handle single record
        return self._process_single_record(record, header_lines)

    def _process_single_record(self, record, header_lines: list = None) -> dict:
        """Process a single pysam.VariantFile record."""
        try:
            import io

            # If header_lines not provided, extract from record
            if header_lines is None:
                header_lines = []
                for header_rec in record.header.records:
                    header_lines.append(str(header_rec))

            # Create a minimal VCF content with just this record
            vcf_content = []

            # Add header lines
            vcf_content.extend(header_lines)

            # Add column header line
            vcf_content.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            # Add the record line
            alt = record.alts[0] if record.alts else "."
            info_str = ";".join(
                [
                    f"{key}={','.join(map(str, val)) if isinstance(val, (list, tuple)) else val}"
                    for key, val in record.info.items()
                ]
            )

            record_line = f"{record.chrom}\t{record.pos}\t.\t{record.ref}\t{alt}\t.\t.\t{info_str}\n"
            vcf_content.append(record_line)

            # Create a StringIO object
            vcf_string = "".join(vcf_content)

            # Use existing _read_ensembl_vcf_to_dataframe logic
            df = self._read_ensembl_vcf_from_string(vcf_string)

            if df.empty:
                return None

            # Apply transformations using existing pipeline
            df = self._apply_transformations(df)

            # Return the first (and only) record as a dict
            return df.iloc[0].to_dict()

        except Exception as e:
            logger.error(f"Error processing record: {e}")
            return None

    def _process_multiple_records(self, records: list, header_lines: list = None) -> list:
        """
        Process multiple pysam.VariantFile records at once for better efficiency.

        Args:
            records: List of pysam.VariantRecord objects
            header_lines: List of VCF header lines (optional, for efficiency)

        Returns:
            list of dict: List of processed record data
        """
        try:
            import io

            if not records:
                return []

            # If header_lines not provided, extract from first record
            if header_lines is None:
                header_lines = []
                for header_rec in records[0].header.records:
                    header_lines.append(str(header_rec))

            # Create VCF content with all records
            vcf_content = []

            # Add header lines
            vcf_content.extend(header_lines)

            # Add column header line
            vcf_content.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            # Add all record lines
            for record in records:
                alt = record.alts[0] if record.alts else "."
                info_str = ";".join(
                    [
                        f"{key}={','.join(map(str, val)) if isinstance(val, (list, tuple)) else val}"
                        for key, val in record.info.items()
                    ]
                )

                record_line = f"{record.chrom}\t{record.pos}\t.\t{record.ref}\t{alt}\t.\t.\t{info_str}\n"
                vcf_content.append(record_line)

            # Create a StringIO object
            vcf_string = "".join(vcf_content)

            # Use existing _read_ensembl_vcf_to_dataframe logic
            df = self._read_ensembl_vcf_from_string(vcf_string)

            if df.empty:
                return []

            # Apply transformations using existing pipeline
            df = self._apply_transformations(df)
            # Return all records as list of dicts
            return df.to_dict(orient="records")

        except Exception as e:
            logger.error(f"Error processing multiple records: {e}")
            return []

    def _create_empty_dataframe(self) -> pd.DataFrame:
        """Create empty DataFrame with expected columns."""
        return pd.DataFrame(columns=self.config.VEP_COLUMN_MAPPING.keys())
    
    def _process_clinical_significance(self, df: pd.DataFrame) -> pd.DataFrame:
        """Split clinical_significance into array-compatible values."""
        if "clinical_significance" in df.columns:
            df["clinical_significance"] = (
                df["clinical_significance"]
                .astype(str)
                .str.strip()
                .str.split(r"\s*,\s*", regex=True)
            )
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

    def process_dataframe(self, df: pd.DataFrame, subtype: str) -> pd.DataFrame:
        """
        Process a DataFrame of Annovar records.

        Args:
            df: DataFrame with Annovar data
            subtype: Type of annotation ("annovar", "intervar", or "cancervar")

        Returns:
            pd.DataFrame: Processed data
        """
        try:
            self._validate_subtype(subtype)

            if df.empty:
                logger.warning("Empty DataFrame provided")
                return self._create_empty_dataframe(subtype)

            return self._apply_transformations(df, subtype)

        except Exception as e:
            logger.error(f"Error processing DataFrame: {e}")
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
        df = self._normalize_intervar_classification(df, subtype)
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

    def _normalize_intervar_classification(self, df: pd.DataFrame, subtype: str) -> pd.DataFrame:
        """Normalize InterVar classification values."""
        if subtype == "intervar" and "intervar_classification" in df.columns:
            df["intervar_classification"] = (
                df["intervar_classification"]
                .astype(str)
                .str.strip()
                .str.lower()
                .str.replace(r"\s+", "_", regex=True)
            )
        return df

    def _create_empty_dataframe(self, subtype: str) -> pd.DataFrame:
        """Create empty DataFrame with expected columns."""
        column_mapping = self._get_column_mapping(subtype)
        columns = ["variant_id", *column_mapping.keys()]
        return pd.DataFrame(columns=columns)


class MetaRnnResultsProcessor:
    """Dedicated class for processing MetaRNN results."""

    def process(self, path: str) -> pd.DataFrame:
        """Process MetaRNN results from the specified file path."""
        try:
            if not os.path.exists(path):
                raise FileNotFoundError(f"MetaRNN results file not found: {path}")

            df = pd.read_csv(path, sep="\t")

            if df.empty:
                logger.warning(f"No data found in MetaRNN results file: {path}")
                return self._create_empty_dataframe()

            df = df.assign(
                TRANSCRIPT=df["TRANSCRIPT"].str.split(";"),
                MetaRNN_score=df["MetaRNN_score"].str.split(";")
            ).explode(
                ["TRANSCRIPT", "MetaRNN_score"]
            )
            return self._apply_transformations(df)

        except Exception as e:
            logger.error(f"Error processing MetaRNN results from {path}: {e}")
            return self._create_empty_dataframe()

    def process_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Process a DataFrame of MetaRNN records.

        Args:
            df: DataFrame with MetaRNN data

        Returns:
            pd.DataFrame: Processed data
        """
        try:
            if df.empty:
                logger.warning("Empty DataFrame provided")
                return self._create_empty_dataframe()

            # Handle transcript and score expansion if needed
            if "TRANSCRIPT" in df.columns and "MetaRNN_score" in df.columns:
                if df["TRANSCRIPT"].dtype == "object" and df["TRANSCRIPT"].str.contains(";", na=False).any():
                    df = df.assign(
                        TRANSCRIPT=df["TRANSCRIPT"].str.split(";"),
                        MetaRNN_score=df["MetaRNN_score"].str.split(";")
                    ).explode(
                        ["TRANSCRIPT", "MetaRNN_score"]
                    )

            return self._apply_transformations(df)

        except Exception as e:
            logger.error(f"Error processing MetaRNN DataFrame: {e}")
            return self._create_empty_dataframe()

    def _apply_transformations(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply all data transformations."""
        df = self._process_variant_identifiers(df)
        df = self._process_metarnn_scores(df)

        df = df.drop(
            columns=["#CHROM", "POS", "REF", "ALT"],
        )

        # Make colum names lowercase
        df.columns = [col.lower() for col in df.columns]

        return df

    def _process_variant_identifiers(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process variant ID formatting."""
        df["variant_id"] = df[["#CHROM", "POS", "REF", "ALT"]].astype(str).agg("_".join, axis=1)
        df["variant_id"] = np.where(
            df["variant_id"].str.startswith("chr"),
            df["variant_id"],
            "chr" + df["variant_id"],
        )
        return df

    def _process_metarnn_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process MetaRNN prediction scores."""
        if "MetaRNN_score" in df.columns:
            df["MetaRNN_score"] = pd.to_numeric(df["MetaRNN_score"], errors="coerce")
            df["MetaRNN_interpretation"] = np.where(
                df["MetaRNN_score"] >= 0.5, "pathogenic", "benign"
            )
        return df

    def _create_empty_dataframe(self) -> pd.DataFrame:
        """Create empty DataFrame with expected columns."""
        return pd.DataFrame(columns=["variant_id", "metarnn_score", "metarnn_interpretation"])