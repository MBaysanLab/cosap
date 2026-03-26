import gzip
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from logging import getLogger
from typing import Tuple

import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

from ..._pipeline_config import PipelineKeys
from ..._utils import (
    convert_list_to_vcf_input,
    join_paths,
)
from ...pipeline_builder.builders import Annotator, VCFReader
from ...tools.annotators import AnnotatorFactory
from ._config import VariantAnnotationConfig
from .annotation_output_parsers import (
    MetaRnnResultsProcessor,
    VEPResultsProcessor,
)

logger = getLogger(__name__)


class VariantMultipleAnnotator:
    def __init__(self, variants: list, workdir: str, max_workers: int = 3):
        self.variants = variants
        self.workdir = (
            workdir if os.path.isabs(workdir) else join_paths(os.getcwd(), workdir)
        )
        self.max_workers = max_workers
        self.vep_input_file = None
        self.vep_output_file = None
        self.metarnn_output_file = None

        # Initialize processors
        self.config = VariantAnnotationConfig()
        self.vep_processor = VEPResultsProcessor(self.config)

    def _is_bgzipped(self, file_path: str) -> bool:
        """Check if a file is bgzipped."""
        try:
            with gzip.open(file_path, "rb") as f:
                f.read(1)
            return True
        except (gzip.BadGzipFile, OSError):
            return False

    def _bgzip_file(self, file_path: str) -> str:
        """Bgzip a file and return the path to the bgzipped file."""
        bgzipped_path = file_path + ".gz"

        # Check if already bgzipped
        if os.path.exists(bgzipped_path):
            logger.info(f"File already bgzipped: {bgzipped_path}")
            return bgzipped_path

        logger.info(f"Bgzipping file: {file_path}")
        try:
            # Use bgzip command
            subprocess.run(
                ["bgzip", "-c", file_path], stdout=open(bgzipped_path, "wb"), check=True
            )
            logger.info(f"Successfully bgzipped: {bgzipped_path}")
            return bgzipped_path
        except subprocess.CalledProcessError as e:
            logger.error(f"Error bgzipping file {file_path}: {e}")
            raise
        except FileNotFoundError:
            logger.error("bgzip command not found. Please install htslib/tabix.")
            raise

    def _tabix_index_file(
        self,
        file_path: str,
        sequence_col: int = 1,
        begin_col: int = 2,
        end_col: int = 2,
    ):
        """Create a tabix index for a bgzipped file."""
        index_path = file_path + ".tbi"

        # Check if index already exists
        if os.path.exists(index_path):
            logger.info(f"Tabix index already exists: {index_path}")
            return

        logger.info(f"Creating tabix index for: {file_path}")
        try:
            subprocess.run(
                [
                    "tabix",
                    "-s",
                    str(sequence_col),
                    "-b",
                    str(begin_col),
                    "-e",
                    str(end_col),
                    file_path,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            logger.info(f"Successfully created tabix index: {index_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating tabix index for {file_path}: {e.stderr}")
            raise
        except FileNotFoundError:
            logger.error("tabix command not found. Please install htslib/tabix.")
            raise

    def _prepare_tabix_indexed_file(
        self, file_path: str, file_type: str = "annovar"
    ) -> str:
        """
        Prepare a file for tabix indexing by bgzipping if needed and creating index.

        Args:
            file_path: Path to the file to prepare
            file_type: Type of file - "annovar" or "vcf"

        Returns:
            str: Path to the bgzipped and indexed file
        """
        if not os.path.exists(file_path):
            logger.error(f"File not found: {file_path}")
            raise FileNotFoundError(f"File not found: {file_path}")

        # Check if file is already bgzipped
        if file_path.endswith(".gz"):
            bgzipped_path = file_path
            if not self._is_bgzipped(file_path):
                logger.warning(
                    f"File has .gz extension but is not properly bgzipped: {file_path}"
                )
                # Remove .gz extension and try again
                original_path = file_path[:-3]
                if os.path.exists(original_path):
                    bgzipped_path = self._bgzip_file(original_path)
                else:
                    raise ValueError(f"Cannot process file: {file_path}")
        else:
            # Check if .gz version already exists
            gz_path = file_path + ".gz"
            if os.path.exists(gz_path) and self._is_bgzipped(gz_path):
                logger.info(f"Bgzipped version already exists: {gz_path}")
                bgzipped_path = gz_path
            else:
                # Bgzip the file
                bgzipped_path = self._bgzip_file(file_path)

        # Check if tabix index already exists
        index_path = bgzipped_path + ".tbi"
        if os.path.exists(index_path):
            logger.info(f"Tabix index already exists: {index_path}")
            return bgzipped_path

        # Create tabix index
        logger.info(f"Creating tabix index for {bgzipped_path}")
        if file_type == "vcf":
            # For VCF files, use -p vcf option
            try:
                subprocess.run(
                    ["tabix", "-p", "vcf", bgzipped_path],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                logger.info(f"Successfully created VCF tabix index: {index_path}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error creating VCF tabix index: {e.stderr}")
                raise
        else:
            # For ANNOVAR files (txt files with chrom, start, end columns)
            # -s 1: chromosome column (1-based)
            # -b 2: begin position column
            # -e 2: end position column
            self._tabix_index_file(
                bgzipped_path, sequence_col=1, begin_col=2, end_col=2
            )

        return bgzipped_path

    def _prepare_annotation_outputs(self):
        """Prepare all annotation output files for tabix indexing."""
        logger.info("Preparing annotation output files for tabix indexing...")

        # Prepare VEP VCF file
        if self.vep_output_file and os.path.exists(self.vep_output_file):
            logger.info("Preparing VEP output file...")
            self.vep_output_file = self._prepare_tabix_indexed_file(
                self.vep_output_file, file_type="vcf"
            )
        else:
            logger.warning(f"VEP output file not found: {self.vep_output_file}")

        # Prepare MetaRNN output file
        if self.metarnn_output_file and os.path.exists(self.metarnn_output_file):
            logger.info("Preparing MetaRNN output file...")
            self.metarnn_output_file = self._prepare_tabix_indexed_file(
                self.metarnn_output_file, file_type="annovar"
            )

        logger.info("All annotation output files prepared for tabix indexing")

    def annotate(self):
        """
        Annotates the variants using VEP and MetaRNN.

        Returns:
            dict: A dictionary containing the annotated variant records.
        """
        vep_input_file = convert_list_to_vcf_input(self.variants, self.workdir)

        self.vep_input_file = vep_input_file

        # Prepare annotation tasks
        annotation_tasks = self._prepare_annotation_tasks(vep_input_file)

        # Execute annotations
        self._execute_annotations(annotation_tasks)

        # Prepare output files for tabix indexing
        self._prepare_annotation_outputs()

    def _prepare_annotation_tasks(self, vep_input_file: str) -> list:
        """Prepare annotation task configurations."""
        vep_input = VCFReader(filename=vep_input_file)

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

        # MetaRNN task
        metarnn_config = Annotator(
            library="metarnn",
            input_step=vep_input,
            sample_name="multiple_anotator",
            input_type="vcf",
            name="multiple_anotator_metarnn",
            output_dir=self.workdir,
        )

        tasks.append(
            {
                "name": "metarnn",
                "config": metarnn_config,
                "annotator": AnnotatorFactory.create("metarnn"),
                "input_file": self.vep_input_file,
            }
        )

        return tasks

    def _execute_annotations(self, tasks: list):
        """Execute annotation tasks in parallel using ThreadPoolExecutor."""
        logger.info(
            f"Running {len(tasks)} annotation tasks with {self.max_workers} workers"
        )

        results = {}

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_task = {
                executor.submit(self._run_annotation_task, task): task for task in tasks
            }

            # Create progress bar
            pbar = tqdm(total=len(tasks), desc="Annotation progress", unit="task")

            # Process completed tasks
            for future in as_completed(future_to_task):
                task = future_to_task[future]
                task_name = task["name"]

                try:
                    output_file = future.result()
                    results[task_name] = output_file
                    logger.info(f"{task_name} annotation completed: {output_file}")

                    # Store output file paths
                    self._store_output_path(task_name, output_file)

                except Exception as e:
                    logger.error(f"Error during {task_name} annotation: {e}")
                    pbar.close()
                    raise
                finally:
                    pbar.update(1)

            pbar.close()

        logger.info(f"All {len(tasks)} annotation tasks completed successfully")
        return results

    def _run_annotation_task(self, task: dict) -> str:
        """Run a single annotation task."""
        name = task["name"]
        config = task["config"]
        annotator = task["annotator"]
        input_file = task["input_file"]

        output_path = join_paths(self.workdir, config.get_output())
        output_path_gz = output_path + ".gz"

        # Check if output already exists (either regular or bgzipped version)
        if os.path.exists(output_path):
            logger.info(
                f"{name.upper()} output already exists: {output_path}, skipping annotation"
            )
            return output_path
        elif os.path.exists(output_path_gz):
            logger.info(
                f"{name.upper()} output already exists (bgzipped): {output_path_gz}, skipping annotation"
            )
            return output_path_gz
        else:
            logger.info(
                f"Annotating {len(self.variants)} variants with {name.upper()} "
                f"using input file: {input_file}"
            )
            annotator.annotate(self._get_annotator_config(config))

            # Return the actual output path (check again which version was created)
            if os.path.exists(output_path_gz):
                return output_path_gz
            elif os.path.exists(output_path):
                return output_path
            else:
                logger.warning(
                    f"Expected output file not found after annotation: {output_path}"
                )
                return output_path

    def _store_output_path(self, task_name: str, output_path: str):
        """Store the output path for each annotation task."""
        if task_name == "vep":
            self.vep_output_file = output_path
        elif task_name == "metarnn":
            self.metarnn_output_file = output_path

    def _fetch_tabix_region(
        self, tbx: pysam.TabixFile, chrom: str, start: int = None, end: int = None
    ) -> list:
        """
        Fetch records from tabix file, handling chromosome naming inconsistencies.

        Args:
            tbx: TabixFile object
            chrom: Chromosome name
            start: Start position (optional)
            end: End position (optional)

        Returns:
            list: List of fetched records
        """
        records = []
        chrom_normalized, chrom_alt = self._normalize_chromosome(chrom)

        # Get available contigs from tabix file
        try:
            contigs = list(tbx.contigs)
        except:
            contigs = []

        # Determine which chromosome naming convention to use
        chrom_to_use = None
        if chrom_normalized in contigs:
            chrom_to_use = chrom_normalized
        elif chrom_alt in contigs:
            chrom_to_use = chrom_alt
        else:
            logger.warning(
                f"Chromosome {chrom} (or {chrom_alt}) not found in tabix file contigs: {contigs[:5]}..."
            )
            return records

        # Fetch records
        try:
            if start is not None and end is not None:
                for row in tbx.fetch(chrom_to_use, start, end):
                    records.append(row.split("\t"))
            elif start is not None:
                for row in tbx.fetch(chrom_to_use, start):
                    records.append(row.split("\t"))
            else:
                for row in tbx.fetch(chrom_to_use):
                    records.append(row.split("\t"))
        except Exception as e:
            logger.warning(f"Could not fetch from {chrom_to_use}: {e}")

        return records

    def _fetch_all_regions_in_batch(
        self,
        tbx: pysam.TabixFile,
        record_batch: list,
    ) -> list:
        """
        Fetch all records from tabix file that match any variant in the batch.
        Handles batches that span multiple chromosomes.

        Args:
            tbx: TabixFile object
            record_batch: List of variant records

        Returns:
            list: List of all fetched records from all chromosomes in the batch
        """
        records = []

        # Group records by chromosome
        chrom_ranges = {}
        for record in record_batch:
            chrom = record.chrom
            pos = record.pos

            if chrom not in chrom_ranges:
                chrom_ranges[chrom] = {"min": pos, "max": pos}
            else:
                chrom_ranges[chrom]["min"] = min(chrom_ranges[chrom]["min"], pos)
                chrom_ranges[chrom]["max"] = max(chrom_ranges[chrom]["max"], pos)

        # Fetch records for each chromosome range
        for chrom, pos_range in chrom_ranges.items():
            chrom_records = self._fetch_tabix_region(
                tbx, chrom, pos_range["min"], pos_range["max"] + 1
            )
            records.extend(chrom_records)

        return records

    def _parse_metarnn_tabix_by_batch(
        self,
        tabix_file: str,
        record_batch: list,
    ) -> pd.DataFrame:
        """
        Parse MetaRNN tabix-indexed tab-separated file for all variants in a batch.
        Handles batches that span multiple chromosomes.

        Args:
            tabix_file: Path to tabix-indexed MetaRNN file
            record_batch: List of pysam.VariantRecord objects

        Returns:
            DataFrame with parsed MetaRNN results
        """
        # Check if file exists and is tabix-indexed
        if not tabix_file or not os.path.exists(tabix_file):
            logger.warning(f"MetaRNN output file not found: {tabix_file}")
            return pd.DataFrame()

        if not os.path.exists(tabix_file + ".tbi"):
            logger.warning(
                f"MetaRNN output file is not tabix-indexed. Preparing file..."
            )
            tabix_file = self._prepare_tabix_indexed_file(
                tabix_file, file_type="annovar"
            )
            self.metarnn_output_file = tabix_file

        try:
            tbx = pysam.TabixFile(tabix_file, encoding="utf-8")

            # Get header from the file
            header = None
            with gzip.open(tabix_file, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        header = line.strip().split("\t")
                        break
                    else:
                        break

            if not header:
                logger.warning(f"No header found in {tabix_file}")
                tbx.close()
                return pd.DataFrame()

            # Fetch records for all chromosomes in the batch
            records = self._fetch_all_regions_in_batch(tbx, record_batch)

            tbx.close()

            if not records:
                # Get chromosome info from batch for logging
                chroms = set(r.chrom for r in record_batch)
                logger.warning(
                    f"No MetaRNN records found for batch spanning chromosomes: {chroms}"
                )
                return pd.DataFrame()

            # Create DataFrame from records
            df = pd.DataFrame(records, columns=header)

            # Use MetaRnnResultsProcessor to process the dataframe
            metarnn_processor = MetaRnnResultsProcessor()
            return metarnn_processor.process_dataframe(df)

        except Exception as e:
            logger.error(f"Error parsing MetaRNN tabix file: {e}")
            return pd.DataFrame()

    def _parse_batch(
        self,
        record_batch: list,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Parse annotation results for a batch of variants.
        Handles batches that span multiple chromosomes.

        Args:
            record_batch: List of pysam.VariantRecord objects from VEP

        Returns:
            Tuple of (vep_df, metarnn_df)
        """
        # Process VEP records
        header_lines = (
            [str(rec) for rec in record_batch[0].header.records] if record_batch else []
        )
        parsed_vep = self.vep_processor.process_record(record_batch, header_lines)
        vep_results = pd.DataFrame(parsed_vep) if parsed_vep else pd.DataFrame()
        vep_results = vep_results[vep_results["vep_pick"] == True]

        # Parse MetaRNN for this batch (handles multiple chromosomes)
        metarnn_results = self._parse_metarnn_tabix_by_batch(
            self.metarnn_output_file,
            record_batch,
        )

        return vep_results, metarnn_results

    def _parse_results_in_batches(self, batch_size: int = 1000):
        """
        Generator that yields batches of parsed annotation results.
        Handles batches that span multiple chromosomes.

        Args:
            batch_size: Number of VEP records to process in each batch

        Yields:
            Tuple of (vep_df, metarnn_df) for each batch
        """
        # Check if VEP file is tabix-indexed
        if not os.path.exists(self.vep_output_file + ".tbi"):
            logger.warning(f"VEP output file is not tabix-indexed. Preparing file...")
            self.vep_output_file = self._prepare_tabix_indexed_file(
                self.vep_output_file, file_type="vcf"
            )

        # Open VEP VCF file
        vcf = pysam.VariantFile(self.vep_output_file)

        # Create an iterator over all records
        record_batch = []

        for record in vcf.fetch():
            record_batch.append(record)

            # Process batch when it reaches batch_size
            if len(record_batch) >= batch_size:
                # Get chromosome range info for logging
                chroms = set(r.chrom for r in record_batch)
                first_record = record_batch[0]
                last_record = record_batch[-1]

                if len(chroms) == 1:
                    logger.info(
                        f"Processing batch: {first_record.chrom}:{first_record.pos}-{last_record.pos} "
                        f"({len(record_batch)} variants)"
                    )
                else:
                    logger.info(
                        f"Processing batch spanning {len(chroms)} chromosomes: "
                        f"{first_record.chrom}:{first_record.pos} to {last_record.chrom}:{last_record.pos} "
                        f"({len(record_batch)} variants)"
                    )

                # Parse this batch
                vep_results, metarnn_results = self._parse_batch(record_batch)

                yield vep_results, metarnn_results

                # Clear memory
                record_batch = []

        # Process remaining records
        if record_batch:
            chroms = set(r.chrom for r in record_batch)
            first_record = record_batch[0]
            last_record = record_batch[-1]

            if len(chroms) == 1:
                logger.info(
                    f"Processing final batch: {first_record.chrom}:{first_record.pos}-{last_record.pos} "
                    f"({len(record_batch)} variants)"
                )
            else:
                logger.info(
                    f"Processing final batch spanning {len(chroms)} chromosomes: "
                    f"{first_record.chrom}:{first_record.pos} to {last_record.chrom}:{last_record.pos} "
                    f"({len(record_batch)} variants)"
                )

            vep_results, metarnn_results = self._parse_batch(record_batch)

            yield vep_results, metarnn_results

        vcf.close()

    def get_results(self, output_file: str = None, batch_size: int = 1000):
        """
        Process and merge annotation results in batches.

        Args:
            output_file: Optional path to save results to file (writes immediately, memory-efficient)
            batch_size: Number of variants to process in each batch

        Returns:
            JSON string of all results (only if output_file is None)
        """
        if output_file:
            # Write directly to file for memory efficiency
            return self.save_results(output_file, batch_size)

        # Legacy behavior: keep in memory and return JSON string
        logger.info(f"Processing annotation results in batches of {batch_size}")
        logger.warning(
            "Loading all results in memory. Consider using output_file parameter for large datasets."
        )

        all_results = []
        batch_count = 0
        total_variants = 0

        for vep_results, metarnn_results in self._parse_results_in_batches(batch_size):
            batch_count += 1

            if vep_results.empty:
                logger.warning(f"Batch {batch_count}: No results found, skipping")
                continue

            # Merge results for this batch
            batch_results = self._merge_annotation_results(
                vep_results,
                metarnn_results,
            )

            # Post-process results
            batch_results = self._post_process_results(batch_results)

            all_results.append(batch_results)
            total_variants += len(batch_results)

            logger.info(
                f"Batch {batch_count} completed: {len(batch_results)} variants "
                f"(total: {total_variants})"
            )

            # Clear memory
            del (
                vep_results,
                metarnn_results,
                batch_results,
            )

        if not all_results:
            logger.warning("No results found in any batch")
            return pd.DataFrame().to_json(orient="records")

        # Concatenate all batch results
        final_results = pd.concat(all_results, ignore_index=True)
        logger.info(f"Total variants processed: {len(final_results)}")
        return final_results.to_json(orient="records")

    def _save_results_as_json(
        self, output_file: str, batch_size: int, batch_count: int, total_variants: int
    ):
        """Save results as gzipped newline-delimited JSON format, writing each batch immediately."""
        # Ensure output file has .gz extension
        if not output_file.endswith(".gz"):
            output_file = output_file + ".gz"

        with gzip.open(output_file, "wt", encoding="utf-8") as f:
            # Create progress bar
            pbar = tqdm(desc="Writing JSON batches", unit="batch", total=None)

            try:
                for vep_results, metarnn_results in self._parse_results_in_batches(
                    batch_size
                ):
                    batch_count += 1

                    if vep_results.empty:
                        logger.error(f"Batch {batch_count}: No results found, skipping")
                        raise ValueError("No results found in batch")

                    # Merge results for this batch
                    batch_results = self._merge_annotation_results(
                        vep_results,
                        metarnn_results,
                    )

                    # Post-process results
                    batch_results = self._post_process_results(batch_results)

                    # Convert to records and write each as a separate JSON line
                    records = batch_results.to_dict(orient="records")

                    for record in records:
                        json.dump(record, f)
                        f.write("\n")

                    variant_count = len(records)
                    total_variants += variant_count

                    # Update progress bar
                    pbar.update(1)
                    pbar.set_postfix(
                        {"variants": total_variants, "batch_size": variant_count}
                    )

                    logger.info(
                        f"Batch {batch_count} saved: {variant_count} variants "
                        f"(total: {total_variants})"
                    )

                    # Flush to disk immediately
                    f.flush()

                    # Clear memory
                    del (
                        vep_results,
                        metarnn_results,
                        batch_results,
                        records,
                    )
            finally:
                pbar.close()

        logger.info(
            f"Annotation results saved to {output_file}: {total_variants} total variants"
        )

    def _save_results_as_ipc(
        self, output_file: str, batch_size: int, batch_count: int, total_variants: int
    ):
        """Save results as Apache Arrow IPC format, writing each batch immediately."""
        import pyarrow as pa
        import pyarrow.ipc as ipc

        # Ensure output file has correct extension
        if not output_file.endswith(".arrow") and not output_file.endswith(".feather"):
            output_file = output_file + ".arrow"

        writer = None
        schema = None

        # Create progress bar
        pbar = tqdm(desc="Writing Arrow IPC batches", unit="batch", total=None)

        try:
            for vep_results, metarnn_results in self._parse_results_in_batches(
                batch_size
            ):
                batch_count += 1

                if vep_results.empty:
                    logger.error(f"Batch {batch_count}: No results found, skipping")
                    raise ValueError("No results found in batch")

                # Merge results for this batch
                batch_results = self._merge_annotation_results(
                    vep_results, metarnn_results
                )

                # Post-process results
                batch_results = self._post_process_results(batch_results)

                # Convert pandas DataFrame to Arrow Table
                arrow_table = pa.Table.from_pandas(batch_results)

                # Initialize writer on first batch
                if writer is None:
                    schema = arrow_table.schema
                    writer = ipc.new_file(output_file, schema)

                # Write batch
                writer.write_table(arrow_table)

                variant_count = len(batch_results)
                total_variants += variant_count

                # Update progress bar
                pbar.update(1)
                pbar.set_postfix(
                    {"variants": total_variants, "batch_size": variant_count}
                )

                logger.info(
                    f"Batch {batch_count} saved: {variant_count} variants "
                    f"(total: {total_variants})"
                )

                # Clear memory
                del (
                    vep_results,
                    metarnn_results,
                    batch_results,
                    arrow_table,
                )
        finally:
            if writer is not None:
                writer.close()
            pbar.close()

        logger.info(
            f"Annotation results saved to {output_file}: {total_variants} total variants"
        )

    def save_results(
        self, output_file: str, batch_size: int = 10_000, format: str = "json"
    ) -> str:
        """
        Save annotation results to a file by processing in batches.
        Reads tabix-indexed files, merges results, and writes immediately.
        This is memory-efficient for large datasets.

        Args:
            output_file: Path to the output file
            batch_size: Number of variants to process in each batch
            format: Output format - "json", "tsv", or "ipc" (Arrow IPC format)

        Returns:
            str: Path to the saved file
        """
        # Ensure the output directory exists
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            logger.info(f"Created output directory: {output_dir}")

        logger.info(
            f"Saving annotation results in batches of {batch_size} to {output_file}"
        )

        total_variants = 0
        batch_count = 0

        format_lower = format.lower()
        if format_lower == "json":
            self._save_results_as_json(
                output_file, batch_size, batch_count, total_variants
            )
        elif format_lower in ["ipc", "arrow", "feather"]:
            self._save_results_as_ipc(
                output_file, batch_size, batch_count, total_variants
            )
        else:
            raise ValueError(
                f"Unsupported format: {format}. Use 'json', 'ipc', or 'arrow'."
            )
        return output_file

    def _merge_annotation_results(
        self,
        vep_results: pd.DataFrame,
        metarnn_results: pd.DataFrame = None,
    ) -> pd.DataFrame:
        """Merge annotation results from VEP and optional MetaRNN."""
        results = vep_results.copy()

        # Merge MetaRNN results if available
        if metarnn_results is not None and not metarnn_results.empty:
            results = results.merge(
                metarnn_results,
                left_on=["variant_id", "feature"],
                right_on=["variant_id", "transcript"],
                how="left",
                suffixes=("", "_metarnn"),
            )

        # Drop the columns that are not needed
        if "transcript" in results.columns:
            results.drop(columns=["transcript"], inplace=True)

        return results

    def _post_process_results(self, results: pd.DataFrame) -> pd.DataFrame:
        """Apply post-processing to merged results."""
        # Make all columns lowercase
        results.columns = [c.lower() for c in results.columns]
        results = self._convert_null_values_to_na(results)

        return results

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

    def _normalize_chromosome(self, chrom: str) -> Tuple[str, str]:
        """
        Normalize chromosome name and return both normalized and alternative versions.

        Args:
            chrom: Chromosome name (e.g., "chr1", "1", "chrX", "X")

        Returns:
            Tuple of (normalized_chrom, alternative_chrom)
            - normalized_chrom: chromosome with 'chr' prefix
            - alternative_chrom: chromosome without 'chr' prefix
        """
        chrom_str = str(chrom)

        if chrom_str.startswith("chr"):
            # If it starts with 'chr', normalized is as-is, alternative is without 'chr'
            normalized = chrom_str
            alternative = chrom_str[3:]
        else:
            # If it doesn't start with 'chr', normalized has 'chr', alternative is as-is
            normalized = f"chr{chrom_str}"
            alternative = chrom_str

        return normalized, alternative
