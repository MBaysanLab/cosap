"""Shared test fixtures and configuration for the cosap test suite."""

import os

import pytest

# Root of the tests/ directory — all data paths are relative to this.
TESTS_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(TESTS_DIR, "data")
PROJECT_ROOT = os.path.abspath(os.path.join(TESTS_DIR, os.pardir))


@pytest.fixture
def data_dir():
    """Return the path to the shared test data directory."""
    return DATA_DIR


@pytest.fixture
def project_root():
    """Return the path to the project root directory."""
    return PROJECT_ROOT
