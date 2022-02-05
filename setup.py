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
    ],
)
