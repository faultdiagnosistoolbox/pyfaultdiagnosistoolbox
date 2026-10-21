import shutil
from pathlib import Path

import pytest

GENERATED_DIR = Path(__file__).resolve().parents[2] / "generated"


@pytest.fixture(scope="session", autouse=True)
def cleanup_generated_dir():
    """Fixture to clean up the generated directory before and after tests."""
    yield
    shutil.rmtree(GENERATED_DIR, ignore_errors=True)
