#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name="cosap",
    version="0.0.2",
    description="COSAP python library for NGS data analysis",
    packages=find_packages(),
    license="MIT license",
    entry_points={
        "console_scripts": [
            "cosap = cosap._cosap:cosap_cli",
        ],
    },
)
