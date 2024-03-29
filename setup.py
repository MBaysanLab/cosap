#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name="cosap",
    version="0.1.0",
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
