import os
from datetime import datetime
from pathlib import Path
from subprocess import run
from typing import Dict, List

import pysam
import yaml

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import AnnotatorKeys
from ._annotators import _Annotatable, _Annotator


class Exomiser(_Annotatable, _Annotator):
    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, app_config: AppConfig, annotator_config: Dict
    ) -> List:
        input_vcf = annotator_config[AnnotatorKeys.INPUT]
        output_file = annotator_config[AnnotatorKeys.OUTPUT]
        hpo_terms = annotator_config[AnnotatorKeys.PARAMS].get("hpo_terms", [])
        output_dir = annotator_config[AnnotatorKeys.OUTPUT_DIR]
        job_file = cls._create_job_file(input_vcf, hpo_terms, output_file, output_dir)

        spring_config_location = (
            os.path.dirname(library_paths.EXOMISER_JAR)
            + "/docker/application.properties"
            if os.path.exists("/.dockerenv")
            else os.path.dirname(library_paths.EXOMISER_JAR) + "/application.properties"
        )

        command = [
            "java",
            "-Xmx4G",
            "-jar",
            library_paths.EXOMISER_JAR,
            "--job",
            job_file,
            f"--spring.config.location={spring_config_location}",
        ]

        return command

    @classmethod
    def _create_job_file(
        cls, input_vcf: str, hpo_terms: list, output_file: str, output_dir: str
    ) -> str:
        sample_name = pysam.VariantFile(input_vcf).header.samples[0]
        phenotypic_features = [{"type": {"id": hpo_id}} for hpo_id in hpo_terms]
        current_time = datetime.utcnow().isoformat() + "Z"

        config = {
            "phenopacket": {
                "id": sample_name.split("_")[0] if "_" in sample_name else sample_name,
                "subject": {"id": sample_name},
                "phenotypicFeatures": phenotypic_features,
                "htsFiles": [
                    {"uri": input_vcf, "htsFormat": "VCF", "genomeAssembly": "hg38"}
                ],
                "metaData": {
                    "created": current_time,
                    "createdBy": "COSAP",
                    "resources": [
                        {
                            "id": "hp",
                            "name": "human phenotype ontology",
                            "url": "http://purl.obolibrary.org/obo/hp.owl",
                            "version": "hp/releases/2025-09-01",
                            "namespacePrefix": "HP",
                            "iriPrefix": "http://purl.obolibrary.org/obo/HP_",
                        }
                    ],
                    "phenopacketSchemaVersion": "1.0",
                },
            },
            "preset": "EXOME",
            "outputOptions": {
                "outputContributingVariantsOnly": False,
                "numGenes": 0,
                "outputDirectory": output_dir,
                "outputFileName": str(Path(output_file).stem),
                "outputFormats": ["JSON", "TSV_VARIANT"],
            },
        }

        job_config_file = str(Path(output_dir) / f"exomiser_job_{sample_name}.yaml")

        try:
            with open(job_config_file, "w") as file:
                yaml.dump(
                    config,
                    file,
                    default_flow_style=False,
                    allow_unicode=True,
                    indent=2,
                    sort_keys=False,
                    width=120,
                )
        except Exception as e:
            raise Exception(f"Error writing Exomiser job file: {e}")

        return job_config_file

    @classmethod
    def annotate(cls, annotator_config: Dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        command = cls.create_command(
            library_paths=library_paths,
            app_config=app_config,
            annotator_config=annotator_config,
        )
        run(command, cwd=annotator_config[AnnotatorKeys.OUTPUT_DIR], check=True)
