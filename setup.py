#!/usr/bin/env python

from setuptools import find_namespace_packages, setup

from cosap import __version__

setup(
    name="cosap",
    version=__version__,
    description="COSAP python library for NGS data analysis",
    packages=find_namespace_packages(),
    package_dir={"cosap": "cosap"},
    package_data={
        "cosap": ["snakemake_workflows/*"],
        "cosap.snakemake_workflows": ["rules/*"],
    },
    license="MIT license",
    entry_points={
        "console_scripts": [
            "cosap = cosap._cosap:cosap_cli",
        ],
    },
)
