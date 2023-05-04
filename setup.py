#!/usr/bin/env python

from distutils.core import setup

setup(
    name="cosap",
    version="0.1",
    description="Cosap",
    packages=[
        "cosap",
        "cosap/mappers",
        "cosap/parsers",
        "cosap/pipeline_builder",
        "cosap/pipeline_runner",
        "cosap/preprocessors",
        "cosap/variant_callers",
        "cosap/pipeline_builder/builders",
        "cosap/scatter_gather",
        "cosap/annotators",
        "cosap/default_pipelines",
        "cosap/quality_controllers",
        "cosap/memory_handler",
        "cosap/genefusion_callers",
        "cosap/celery",
        "cosap/msi_callers",
        "cosap/cnv_callers",
    ],
    license="MIT license",
    entry_points="""
        [console_scripts]
        cosap=cosap._cosap:cosap_cli
    """,
)
