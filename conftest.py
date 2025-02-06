# conftest.py
import pytest

def pytest_configure(config):
    """Register custom marks"""
    config.addinivalue_line(
        "markers",
        "integration: mark test as an integration test"
    )