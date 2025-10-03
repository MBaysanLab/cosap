import os
from logging import getLogger
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Tuple

import numpy as np
import pandas as pd

from ..._pipeline_config import PipelineKeys
from ..._utils import (
    convert_list_to_annovar_input,
    convert_list_to_vcf_input,
    join_paths,
)
from ...pipeline_builder.builders import Annotator, VCFReader
from ...tools.annotators import AnnotatorFactory
from ._config import VariantAnnotationConfig
from .annotation_output_parsers import VEPResultsProcessor, AnnovarResultsProcessor

logger = getLogger(__name__)


class VariantMultipleAnnotator:
    def __init__(self, variants: list, workdir: str, max_workers: int = 3):
        self.variants = variants
        self.workdir = workdir
        self.max_workers = max_workers
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
        Annotates the variants using VEP, Intervar, and Cancervar in parallel.

        Returns:
            dict: A dictionary containing the annotated variant records.
        """
        vep_input_file = convert_list_to_vcf_input(self.variants, self.workdir)
        annovar_input_file = convert_list_to_annovar_input(self.variants, self.workdir)

        self.vep_input_file = vep_input_file
        self.annovar_input_file = annovar_input_file

        # Prepare annotation tasks
        annotation_tasks = self._prepare_annotation_tasks(
            vep_input_file, annovar_input_file
        )

        # Execute annotations in parallel
        self._execute_parallel_annotations(annotation_tasks)

    def _prepare_annotation_tasks(
        self, vep_input_file: str, annovar_input_file: str
    ) -> list:
        """Prepare annotation task configurations."""
        vep_input = VCFReader(filename=vep_input_file)
        annovar_input = VCFReader(filename=annovar_input_file)

        tasks = []

        # VEP task
        vep_config = Annotator(
            library="vep",
            input_step=vep_input,
            sample_name="multiple_anotator",
            name="multiple_anotator_vep",
            input_type="vcf",
            output_dir=self.workdir,
        )
        tasks.append(
            {
                "name": "vep",
                "config": vep_config,
                "annotator": AnnotatorFactory.create("vep"),
                "input_file": self.vep_input_file,
            }
        )

        # Intervar task
        intervar_config = Annotator(
            library="intervar",
            input_step=annovar_input,
            sample_name="multiple_anotator",
            input_type="AVinput",
            name="multiple_anotator_intervar",
            output_dir=self.workdir,
        )
        tasks.append(
            {
                "name": "intervar",
                "config": intervar_config,
                "annotator": AnnotatorFactory.create("intervar"),
                "input_file": self.annovar_input_file,
            }
        )

        # Cancervar task
        cancervar_config = Annotator(
            library="cancervar",
            input_step=annovar_input,
            sample_name="multiple_anotator",
            input_type="AVinput",
            name="multiple_anotator_cancervar",
            output_dir=self.workdir,
        )
        tasks.append(
            {
                "name": "cancervar",
                "config": cancervar_config,
                "annotator": AnnotatorFactory.create("cancervar"),
                "input_file": self.annovar_input_file,
            }
        )

        return tasks

    def _execute_parallel_annotations(self, tasks: list):
        """Execute annotation tasks in parallel using ThreadPoolExecutor."""
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_task = {
                executor.submit(self._run_annotation_task, task): task for task in tasks
            }

            # Process completed tasks
            for future in as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    output_file = future.result()
                    logger.info(f"{task['name']} annotation completed: {output_file}")

                    # Store output file paths
                    self._store_output_path(task["name"], output_file)

                except Exception as e:
                    logger.error(f"Error during {task['name']} annotation: {e}")
                    raise

    def _run_annotation_task(self, task: dict) -> str:
        """Run a single annotation task."""
        name = task["name"]
        config = task["config"]
        annotator = task["annotator"]
        input_file = task["input_file"]

        output_path = join_paths(self.workdir, config.get_output())

        if not os.path.exists(output_path):
            logger.info(
                f"Annotating {len(self.variants)} variants with {name.upper()} "
                f"using input file: {input_file}"
            )
            annotator.annotate(self._get_annotator_config(config))
        else:
            logger.info(f"{name.upper()} output already exists, skipping annotation")

        return output_path

    def _store_output_path(self, task_name: str, output_path: str):
        """Store the output path for each annotation task."""
        if task_name == "vep":
            self.vep_output_file = output_path
        elif task_name == "intervar":
            self.intervar_output_file = output_path
        elif task_name == "cancervar":
            self.cancervar_output_file = output_path

    def get_results(self):
        """Process and merge annotation results."""
        # Parse results in parallel
        vep_results, intervar_results, cancervar_results = (
            self._parse_results_parallel()
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

        # Merge results
        results = self._merge_annotation_results(
            vep_results, intervar_results, cancervar_results
        )

        # Post-process results
        results = self._post_process_results(results)

        return results.to_json(orient="records")

    def _parse_results_parallel(
        self,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Parse annotation results in parallel."""
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit parsing tasks
            vep_future = executor.submit(
                self._parse_ensembl_vep_results, self.vep_output_file
            )
            intervar_future = executor.submit(
                self._parse_annovar_results, self.intervar_output_file, "intervar"
            )
            cancervar_future = executor.submit(
                self._parse_annovar_results, self.cancervar_output_file, "cancervar"
            )

            # Wait for all to complete
            vep_results = vep_future.result()
            intervar_results = intervar_future.result()
            cancervar_results = cancervar_future.result()

        return vep_results, intervar_results, cancervar_results

    def _merge_annotation_results(
        self,
        vep_results: pd.DataFrame,
        intervar_results: pd.DataFrame,
        cancervar_results: pd.DataFrame,
    ) -> pd.DataFrame:
        """Merge annotation results from different sources."""
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

        return results

    def _post_process_results(self, results: pd.DataFrame) -> pd.DataFrame:
        """Apply post-processing to merged results."""
        # Make all columns lowercase
        results.columns = [c.lower() for c in results.columns]
        results = self._convert_null_values_to_na(results)

        return results

    def _parse_ensembl_vep_results(self, path: str) -> pd.DataFrame:
        """Parses the output of ensembl vep and returns a list of variants."""
        return self.vep_processor.process(path)

    def _parse_annovar_results(self, path: str, subtype: str) -> pd.DataFrame:
        """Parses the output of annovar and returns a list of variants."""
        annovar_processor = AnnovarResultsProcessor(self.config)
        return annovar_processor.process(path, subtype)

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
            df.replace(["None", "none", "NA", "N/A", "na", ".", "-", "", np.nan], pd.NA)
        ).fillna(pd.NA)
