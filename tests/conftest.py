"""Shared fixtures for seqproc paper analysis tests."""

import sys
import os
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent.parent


@pytest.fixture(scope="session", autouse=True)
def scripts_on_path():
    """Add scripts/ to sys.path so test modules can import them."""
    scripts_dir = str(PROJECT_ROOT / "scripts")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)


@pytest.fixture(scope="session")
def project_root():
    return PROJECT_ROOT
