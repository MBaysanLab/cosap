#!/usr/bin/env python

from setuptools import find_packages, setup
from cosap._version import version

setup(
    name="cosap",
    version=version,
    description="COSAP python library for NGS data analysis",
    packages=find_packages(),
    package_data={"cosap": ["snakemake_workflows/*"]},
    license="MIT license",
    entry_points={
        "console_scripts": [
            "cosap = cosap._cosap:cosap_cli",
        ],
    },
)
