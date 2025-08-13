import io
import os

import pandas as pd

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


class VariantMultipleAnnotator:
    def __init__(self, variants: list, workdir: str):
        self.variants = variants
        self.workdir = workdir
        self.vep_input_file = None
        self.vep_output_file = None
        self.annovar_input_file = None
        self.intervar_output_file = None
        self.cancervar_output_file = None

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
        results = self._convert_null_values_to_nan(results)

        # # Remove temp input files
        # os.remove(self.vep_input_file)
        # os.remove(self.annovar_input_file)

        return results.to_json(orient="records")

    def _parse_ensembl_vep_results(self, path: str) -> pd.DataFrame:
        """
        Parses the output of ensembl vep and returns a list of variants.
        """

        vep_column_mapping = {
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
            "pubmed": "PUBMED",
        }

        # Skip rows that start with ##
        df = self._read_ensembl_vcf_to_dataframe(path)
        df = df.rename(columns={v: k for k, v in vep_column_mapping.items()})

        df = df[[col for col in vep_column_mapping.keys() if col in df.columns]]

        # Convert chr_pos_ref/alt to chr_pos_ref_alt
        df["variant_id"] = df["variant_id"].str.replace("/", "_")

        # Convert cannonical and pick to boolean
        df["canonical"] = df["canonical"].str.lower() == "yes"
        df["vep_pick"] = df["vep_pick"] == "1"

        # Consequence column is a list of consequences separated by &, \
        # get first consequence since it is the most severe
        df["consequence"] = df["consequence"].str.split("&").str[0]
        return df

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

    def _convert_null_values_to_nan(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Converts all None values to NaN in a dataframe.

        None values can be "-", ".", or any other string that represents a null value.
        """
        return df.replace("-", pd.NA).replace(".", pd.NA).replace("", pd.NA)
