# conftest.py

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run_tb", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "tb: mark test as involving M. tuberculosis and therefore slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run_tb"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_tb = pytest.mark.skip(reason="need --run_tb option to run")
    for item in items:
        if "tb" in item.keywords:
            item.add_marker(skip_tb)
